//local drive Z:\
local drive /Users/Michael/Documents

cd "`drive'/stdev/gml"
adopath ++ "./gml"
//adopath ++ "Z:\stdev\gml\ssc\version_0_0_4_temp"
clear all
do ./build/buildmlib.do

set seed 7254
pr drop _all
clear
set obs 500
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal()
gen u3 = rnormal()
gen trt = runiform()>0.5
survsim stime died, hazard(0.1:*1.2:*#t:^0.2 :* 	///
	exp(0.2:*(u1 :+ (0.1:+u2):*#t :+ 0.1:* #t:^2 :+ 0.5:*trt) ///
	)) cov(trt -0.5) maxt(5)

expand 5
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u1 + (0.1+u2)*time + 0.1*time^2+ 0.5*trt
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

//gsem (stime <- trt , family(weib, failure(died)) ) ///
//	(y <- trt time M2[id1],family(gaussian) ) ///
//	, intmethod(mvagh) intpoints(7) adaptopts(log) 

timer on 2
gml (stime trt EV[y]@alpha1, family(weib, failure(died)) timevar(stime)) ///
	(y trt {&t}@lam1 {&t}#M1[id1] M2[id1],family(gaussian) timevar(time)) ///
	, intmethod(mvagh) intpoints(7) adaptopts(log) //showinit //alpha1(trt)
timer off 2
timer list
	
	
	
