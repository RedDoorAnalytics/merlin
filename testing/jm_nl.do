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
set obs 200
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal(0,0.1)
gen trt = runiform()>0.5
survsim stime died, hazard(0.1:*1.2:*#t:^0.2 :* 		///
					exp(								///
					0.2:*(u1 :+ (0.1:+u2):*#t)			///
					:+ 0.01 :* (u1 :+ (0.1:+u2):*#t):^2	///
					)) cov(trt -0.5) maxt(5)

expand 5
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u1 + (0.1+u2)*time 
gen y = rnormal(xb,0.5)

bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1

merlin 	(y 										///
				fp(time, pow(1)) 				///
				fp(time, pow(1))#M1[id1]@1 		///
				M2[id1]@1						///
				, family(gaussian) 				///
				timevar(time)) 					///
		(stime 									///
				trt 							///
				rcs(EV[y], df(3))				///
				, family(weib, failure(died)) 	///
				timevar(stime)) 
	

	
