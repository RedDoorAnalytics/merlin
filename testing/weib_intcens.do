
//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
adopath ++ "`drive'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 72549

pr drop _all
clear

set obs 100000
gen trt = runiform()>0.5


survsim stime event , dist(weib) lambda(0.1) gamma(1.2) covariates(trt -0.5) maxt(5)

local maxt 5

gen st1 = floor(stime) if stime < `maxt'
replace stime = st1 + 1 if st1 < `maxt'
replace event = 2 if event == 1

// replace st1 = 0 if st1>0 & st1!=.
timer clear
timer on 1
merlin 	(stime trt , family(rp, df(1) noorthog failure(event) linterval(st1))) , evaltype(gf0)
timer off 1
timer on 2
merlin 	(stime trt , family(rp, df(1) noorthog failure(event) linterval(st1))) , evaltype(gf2) 
timer off 2
// timer list

replace st1 = stime if event==0
replace stime = . if event==0
timer on 3
stintreg trt, dist(weib) interval(st1 stime) nohr
timer off 3
timer list

// gen died = event>0
// stset stime, f(died)
// replace st1 = _t if st1==.
// stpm trt, scale(h) df(3) left(st1)
