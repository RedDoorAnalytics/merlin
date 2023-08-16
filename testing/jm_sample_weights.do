//simulate a large joint model dataset
//sample 5% and then upweight everyone

//local drive Z:
local drive /Users/Michael/Documents
cd "`drive'/megenreg/Stata"
adopath ++ "./megenreg"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 1000
gen id1 = _n
gen u1 = rnormal(0,1)
gen u2 = rnormal(0,0.1)
gen trt = runiform()>0.5

survsim stime died, hazard(0.1:*1.2:*#t:^0.2 :* 					///
							exp(									///
								0.2 :* (u1 :+ (0.1:+u2):*#t)		///
								)) cov(trt -0.5) maxt(5)
expand 5
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u1 + (0.1+u2)*time
gen y = rnormal(xb,0.5)

gen start = time
bys id: gen stop = start[_n+1]
gen event = 0
bys id : replace event = died if _n==_N
bys id: replace stop = stime if _n==_N
stset stop, enter(start) id(id) f(event)

bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1

gen wts1 = 1
gen wts2 = 10 

megenreg	(stime trt EV[2]@a1 				, family(weibull, failure(died))) ///
			(y fp(1)@lam1 fp(1)#M2[id1] M1[id1] , family(gaussian) timevar(time)) ///
			, weights(wts1 wts2)
	
	
