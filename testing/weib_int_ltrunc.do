//local drive Z:/
local drive /Users/Michael
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 72549

clear
set obs 100000
gen id 	= _n
gen trt = runiform()>0.5
gen u1 	= rnormal(0,1) 

// expand 10

survsim stime died, dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 u1 1) maxt(10)

gen t0 = 0
replace t0 = runiform() * 5 if _n<=50000

drop if stime<t0

merlin (stime trt M1[id]@1, family(weib, failure(died) 		///
					ltruncated(t0, marginal))) ,	///
	adaptopts(log) intpoints(11)



//need to handle mix of t0=0 and del entry

stset stime, f(died) id(id) enter(t0)

gsem (_t <- trt M1[id]@1 , family(weibull, failure(_d) ltruncated(_t0))) 


merlin (_t trt M1[id]@1, family(w, failure(_d)          ///
                                ltruncated(_t0)))       ///
        , //adaptopts(log)
