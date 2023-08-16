//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear


set seed 7254
clear
set obs 100
gen id1 = _n
expand 3
bys id1: gen id2 = _n
expand 4
bys id1 id2: gen id3 = _n
expand 5
gen id4 = _n

gen trt = runiform()>0.5
bys id1 (id2 id3 id4) : gen double u1 = rnormal() if _n==1
bys id1 (id2 id3 id4) : replace u1 = u1[1]
bys id1 id2 (id3 id4) : gen double u2 = rnormal() if _n==1
bys id1 id2 (id3 id4) : replace u2 = u2[1]
bys id1 id2 id3 (id4) : gen double u3 = rnormal() if _n==1
bys id1 id2 id3 (id4) : replace u3 = u3[1]

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 u1 1 u2 1 u3 1) maxt(5)
stset stime, f(dead)
timer clear
timer on 1
	merlin (stime trt M1[id1]@1 M2[id1>id2]@1 M3[id1>id2>id3]@1, family(weibull, failure(dead))) ///
		, //intmethod(gh)
timer off 1

timer on 2
// gsem (stime <- trt M1[id1] M2[id2<id1] M3[id3<id2<id1], family(weibull, failure(dead))) , ///
// 	adaptopts(log)
timer off 2
timer list
