//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/forge"
adopath ++ "./forge"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 1000
gen id = _n
gen u1 = rnormal(0,0.5)
gen u2 = rnormal(0,0.1)
gen trt = runiform()>0.5

expand 5
bys id : gen time = _n-1

gen xb = u1 + 0.1*time
gen cost = rgamma(xb,1.5) + 1

gen qaly = rnormal(xb + 0.5 * cost,1)

//gsem (cost <- time M1[id], family(gamma))



forge 	(cost time M1[id]@1 			,family(gamma)) 			///
// 		(qaly time M2[id]@1 EV[cost] 	,family(gaussian))

	
