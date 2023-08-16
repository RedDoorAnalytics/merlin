//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
adopath ++ "./jm"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 200
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal()
gen u3 = rnormal()
gen trt = runiform()>0.5
survsim stime died, hazard(0.1:*1.2:*{t}:^0.2 :* 	///
	exp(0.2:*(u1 :+ (0.1:+u2):*{t} :+ 0.5:*trt) :+ ///
		0.1:*(u1:*{t} :+ 0.5:*(0.1:+u2):*{t}:^2 :+ 0.5:*trt:*{t})	///
		:- 0.2:*(u1 :+ (0.1:+u2):*{t} :+ 0.5:*trt) :*(u3 :+ (0.2):*{t} ) ///
	)) cov(trt -0.5) maxt(5)

expand 5
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u1 + (0.1+u2)*time + 0.5*trt
gen y = rnormal(xb,0.5)
gen xb2 = u3 + (0.2)*time 
gen y2 = rnormal(xb2,0.5)

gen start = time
bys id: gen stop = start[_n+1]
gen event = 0
bys id : replace event = died if _n==_N
bys id: replace stop = stime if _n==_N
stset stop, enter(start) id(id) f(event)
timer clear
timer on 1
//stjm y trt, panel(id) survm(weib) rfp(1) gh(7) survcov(trt) cov(indep) 
timer off 1
cap est store m1
bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1

timer on 2

merlin 	(y trt fp(time, pow(1)) fp(time, pow(1))#M1[id1]@1 M2[id1]@1,family(gaussian) timevar(time)) ///
		(y2 fp(time, pow(1)) M3[id1]@1,family(gaussian) timevar(time)) ///
		(stime trt EV[y] iEV[y] EV[y]#EV[y2], family(weib, failure(died)) timevar(stime)) ///
		, intmethod(mvagh) intpoints(7) adaptopts(log)
		
timer off 2
timer list
	
	
	
