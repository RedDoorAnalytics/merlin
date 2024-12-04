//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"

clear all
tr:do ./build/buildmlib.do
mata mata clear

pr drop _all

set seed 7254
clear
set obs 1000
gen id1 = _n
gen u0 = rnormal(0,1)
gen u1 = rnormal(0,0.3)
gen trt = runiform()>0.5

survsim stime died, hazard(0.1:*1.2:*{t}:^0.2 :* /// baseline hazard
        exp(0.1 :* (0.5 :*trt :+ u0 :+ (0.1:+u1):*{t} :+ 0.2:*{t}:^2) ///
		:+ 0.3 :* (0.5 :*trt :+ u0))) 	///
        covariates(trt -0.5) 	///
        maxt(10)		//	admin. censoring

expand 10
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u0 + (0.1+u1) * time + 0.2 * time^2 + 0.5 * trt
gen y = rnormal(xb,0.5)

bys id1 (time) : replace stime 	= . if _n>1
bys id1 (time) : replace died 	= . if _n>1
	
merlin 	(y trt time fp(time, pow(2)) time#M2[id]@1 M1[id]@1, 	///
		family(gaussian) timevar(time))			///
	(stime trt EV[y, time(0)] XB[y], 			///
		family(w, failure(died)) timevar(stime))	///
		, cov(unstr)

merlin 	(y trt time fp(time, pow(2)) time#M2[id]@1 M1[id]@1, 	///
		family(gaussian) timevar(time))			///
	(stime trt EV[y, time(0)], 				///
		family(w, failure(died)) timevar(stime))	///
		, cov(unstr)

merlin 	(y trt time fp(time, pow(2)) time#M2[id]@1 M1[id]@1, 	///
		family(gaussian) timevar(time))			///
	(stime trt EV[y, time(0)], 				///
		family(w, failure(died)))			///
		, cov(unstr)
