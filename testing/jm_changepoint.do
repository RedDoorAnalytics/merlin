//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 300
gen id1 = _n
gen u1 = rnormal()
gen u2 = exp(rnormal(.69314718,1))
gen trt = runiform()>0.5

survsim stime died, hazard(0.1:*1.2:*#t:^0.2 :* 					///
							exp(									///
								0.2 :* (u1 :+ 0.1:*#t :+ (-0.3):*(#t:>2):*#t) 		///
								)) cov(trt -0.5) maxt(5)

expand 10
bys id1 : gen time = _n-1
//drop if time>stime

gen xb = u1 + 1*time - 2 * (time > 2) * time
gen y = rnormal(xb,0.5)

gen start = time
bys id: gen stop = start[_n+1]
gen event = 0
bys id : replace event = died if _n==_N
bys id: replace stop = stime if _n==_N
stset stop, enter(start) id(id) f(event)

bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1
drop if time>10

//change point at 2 years
mata:
real colvector tfun(gml,t)
{
	return((t:>2):*t)
}
end

//estimate changepoint
mata:
real matrix tfun2(gml,t)
{
	cp = merlin_util_ap(gml,1)
	return((t:>cp):*t)
}
end

// merlin	(y fp(time, pow(1)) mf(tfun) M1[id1]@1, family(gaussian) timevar(time)) 	///
// 		(stime EV[y] , family(weibull, failure(died)))								///
// 		, 
mat init = J(1,5,1)

merlin	/*(stime EV[y] , family(weibull, failure(died)) timevar(stime))*/										///
		(y fp(time, pow(1)) mf(tfun2), family(gaussian, nap(1)) timevar(time)) 	///
		, apstartvalues(1) from(init)
		
	
