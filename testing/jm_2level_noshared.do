
local drive Z:\

cd "`drive'\stdev\gml"
adopath ++ ".\gml"
//adopath ++ "Z:\stdev\gml\ssc\version_0_0_4_temp"
clear all
do ./build/buildmlib.do

set seed 7254
pr drop _all
clear
set obs 500
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal()
gen u3 = rnormal()
gen u4 = rnormal()

gen trt = runiform()>0.5
survsim stime died, dist(weib) lambda(0.21) gamma(1.2) cov(trt -0.5) maxt(5)

expand 5
bys id1 : gen time = _n-1

drop if time>stime

gen xb = u1 + time + time * u2 + 0.5*trt
gen xb2 = u3 + time + /*time * u4 +*/ 0.5*trt
gen y = rnormal(xb,0.5)
gen y2 = rnormal(xb2,0.5)

bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1


timer clear
timer on 1
gml (stime trt , family(weib, failure(died)))	///
 	(y time trt M1[id1] time*M2[id1], family(gaussian))		///
	/*(y2 time trt M3[id1] , family(gaussian))*/, 			///
	intmethod(mvagh) cov(unstr) intpoints(7)
timer off 1
est store m1
timer on 2
gsem (stime <- trt , family(weib, failure(died)))	///
 	(y <- time trt M1[id1]@1 c.time#M2[id1], family(gaussian))		///
	/*(y2 time trt M3[id1] , family(gaussian))*/, 			///
	intmethod(mvagh) intpoints(7) evaltype(gf0) trace

timer off 2
est store m2
timer list


est tab m1 m2
est stats m1 m2
