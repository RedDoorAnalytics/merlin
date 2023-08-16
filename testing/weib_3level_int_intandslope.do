//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/megenreg/Stata"
adopath ++ "./megenreg"
clear all

do ./build/buildmlib.do
mata mata clear


set seed 7254
clear all

set obs 300
gen id1 = _n
expand 7
bys id1: gen id2 = _n
expand 7	
gen id3 = _n
gen trt = runiform()>0.5
bys id1 (id2 id3) : gen u1 = rnormal() if _n==1
bys id1 (id2 id3) : replace u1 = u1[1]
bys id1 id2 (id3) : gen u2 = rnormal() if _n==1
bys id1 id2 (id3) : replace u2 = u2[1]
bys id1 id2 (id3) : gen u3 = rnormal() if _n==1
bys id1 id2 (id3) : replace u3 = u3[1]
gen trtui = (-0.5+u3) * trt


survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(u1 1 u2 1 trtui 1) maxt(5)
stset stime, f(dead)

timer clear
timer on 1
megenreg (stime trt M1[id1] M2[id1>id2] trt#M3[id1>id2], family(weibull, failure(dead))) ///
			, intmethod(gh) intpoints(7) 
timer off 1
timer list

