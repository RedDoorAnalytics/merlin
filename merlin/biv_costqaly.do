//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/megenreg/Stata"
adopath ++ "./megenreg"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 30					//number of clusters
gen id = _n
gen u1 = rnormal(0,1)
gen u2 = rnormal(0,1)

expand 20					//patients per cluster
gen trt = runiform()>0.5

//lin pred for cost
gen xb1 = 100 + u1 -20 * trt
//observed cost
gen cost = rnormal(xb1,1)


//lin pred for qaly
gen xb2 = u2 + 2 * trt + 1 * xb1
gen qaly = rnormal(xb2,1)

megenreg 	(cost trt M1[id] 					, family(gamma)) 	///
			(qaly trt M2[id] EV[cost]@theta 	, family(gaussian))		
	
