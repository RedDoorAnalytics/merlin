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
set obs 100

gen trt = runiform()>0.5
gen age = rnormal(50,5)
survsim stime died, dist(weib) lambda(0.001) gamma(1.2) cov(trt -0.5 age 0.1) maxt(5) tde(trt 0.01)

// replace stime = stime[2] if _n==1

gen id = _n
stset stime, f(died) id(id)

// set seed 87876
// bsurvci, generate(scurve) at(trt 1) id(id) reps(200)

merlin (stime trt age trt#rcs(stime, df(1) log), 	///
			family(cox, failure(died)) timevar(stime))  

predict s1, surv at(trt 1)
predict s2, surv at(trt 1) standardise ci reps(100)
