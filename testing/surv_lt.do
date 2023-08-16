//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 72549

clear
set obs 500
gen id 	= _n
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5) maxt(10)

gen t0 = 0
replace t0 = runiform() * 5 if _n>200

drop if stime<t0

// merlin (stime trt , family(w, ltruncated(t0) failure(died))) 

stset stime, enter(t0) f(died)
stcox trt, nohr tvc(trt) texp(log(_t))

merlin (stime trt trt#rcs(stime, df(1) log), family(cox, ltruncated(t0) failure(died)) timevar(stime)) 
