//source paths
local drive /Users/michael/Library/CloudStorage
local drive `drive'/OneDrive-RedDoorAnalyticsAB/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"


//build mlib
clear all
do ./build/buildmlib.do
mata mata clear

pr drop _all

clear 
set seed 725
set obs 1000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen bmi = rnormal(30,3)

gen t0 = 0 
replace t0 = 3*runiform() //in 1/100

survsim stime died , 	dist(weib) lambda(0.1) gamma(1.2)       ///
			cov(trt -0.5 age 0.01 bmi -0.05)        ///
                        maxt(10) //ltruncated(t0)	

mata:
real matrix gw_logh(transmorphic gml, real matrix t)
{
	logalpha = merlin_util_xzb(gml)
	beta  = exp(merlin_util_ap(gml,1))
	return(logalpha :+ log(t) :- log(beta:^2 :+ t:^2))
}
real matrix gw_logh2(transmorphic gml, real matrix t)
{
	logalpha = merlin_util_xzb(gml)
	beta  = exp(merlin_util_xzb_mod(gml,2))
	return(logalpha :+ log(t) :- log(beta:^2 :+ t:^2))
}
end

clear
set obs 100000
gen trt = runiform()>0.5
local alpha 1
local lnbeta 2 + 0.5 * trt
survsim stime died, hazard(`alpha' :* {t} :/ (exp(`lnbeta'):^2 :+ {t}:^2))   ///
                        cov(trt -0.5) maxtime(20)
merlin (stime trt , family(user, failure(died) loghfunction(gw_logh2))) 		///
	(trt, family(null))
est store m1

exit

range tvar 0 20 100
predict h0, hazard zeros timevar(tvar)
predict s0, survival zeros timevar(tvar)
predict h1, hazard timevar(tvar) at(trt 1)
predict s1, survival timevar(tvar) at(trt 1)

predictms , singleevent model(m1) hazard survival timevar(tvar)
