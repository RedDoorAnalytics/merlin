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
set obs 1000

gen trt = runiform()>0.5

gen id = _n
gen t0 = runiform()*5
survsim stime died, hazard(0.1:*1.2:*{t}:^0.2 :* exp(0.1:*{t}:*trt)) cov(trt -0.5) maxt(20) //ltruncated(t0)

qui stset stime, f(died) id(id)

// stpm2 trt, df(3) scale(h) //tvc(trt) dftvc(1) 
// staft trt, df(3) tvc(trt) dftvc(1)

stsplit _new, at(10)

merlin (stime 	trt , family(weib, failure(died)))  	///
		, 
