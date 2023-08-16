local drive Z:\

cd "`drive'\stdev\gml"
adopath ++ ".\gml"
clear all

do ./build/buildmlib.do

set seed 7254
pr drop _all
clear
set obs 300
gen id1 = _n
expand 5
bys id1: gen id2 = _n
bys id1 (id2) : gen u1 = rnormal(0,0.5) if _n==1
bys id1 (id2) : replace u1 = u1[1]
bys id1 id2 : gen u2 = rnormal(0,0.5) if _n==1
bys id1 id2 : replace u2 = u2[1]
gen trt = runiform()>0.5

survsim stime died , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 u1 1 u2 1) maxt(5)

expand 5
bys id1 id2: gen id3 = _n
bys id1 id2 : gen time = _n-1
drop if time>stime
bys id1 id2 : replace stime = . if _n<_N

gen double y = u1 + u2 -0.5* trt + 0.2*time + rnormal(0,0.5)

timer clear
timer on 1
gml (stime trt M1[id1]@lam1 M2[id1>id2]@lam2, family(weib, failure(died)))	///
 	(y time trt M1[id1] M2[id1>id2], family(gaussian)),		///
	intmethod(gh) intpoints(7)
timer off 1
est store m1

timer on 2
gsem (stime <- trt M1[id1]@lam1 M2[id1>id2]@lam2, family(weib, failure(died)))	///
 	(y <- time trt M1[id1]@1 M2[id1>id2]@1, family(gaussian)),		///
	intmethod(gh) intpoints(7) //evaltype(gf0)
timer off 2
est store m2

timer list

/*
gsem (stime <- trt M1[id1]@lam1 M2[id1>id2]@lam2, family(weib, failure(died)))	///
 	(y <- time trt M1[id1]@1 M2[id1>id2]@1, family(gaussian)),		///
	intmethod(mvagh) intpoints(7)
//est store m2
*/

