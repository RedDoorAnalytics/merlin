//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

pr drop _all
set seed 7254
clear
set obs 10000

gen trt = runiform()>0.5
gen t0 = 3*runiform()
survsim stime died, hazard(({t}:<1):*0.1 :+ ({t}:>=1):*0.2) cov(trt -0.5) maxt(10) ltruncated(t0)

qui stset stime, f(died)

timer clear
timer on 1
merlin (stime trt trt#fp(stime, pow(0)), family(pwe, knots(1 2 3) failure(died) ltruncated(t0))), ///
			evaltype(gf0)
timer off 1

timer on 2
merlin (stime trt trt#fp(stime, pow(0)), family(pwe, knots(1 2 3) failure(died) ltruncated(t0))),	///
		evaltype(gf1) 
timer off 2

timer on 3
merlin (stime trt , family(pwe, knots(2) failure(died) ltruncated(t0))),	///
		evaltype(gf2debug)
timer off 3

// est store m1 

// predict s1, survival zeros 
// predict h1, haz zeros 

// predictms, singleevent survival at1() models(m1)

timer list
