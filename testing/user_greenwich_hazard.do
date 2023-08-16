//local drive Z:/
local drive /Users/Michael/Documents/merlin

cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

tr:do ./build/buildmlib.do
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
end


merlin (stime trt rcs(age,df(3)) bmi ,  ///
        family(user, failure(died) loghfunction(gw_logh) nap(1))) 	

clear
set obs 1000
gen trt = runiform()>0.5
local alpha = 1
local beta = 10
survsim stime died, hazard(`alpha' :* {t} :/ ((`beta'):^2 :+ {t}:^2))   ///
                        cov(trt -0.5) maxtime(20)
merlin (stime trt , family(user, failure(died) loghfunction(gw_logh) nap(1))) 	
est store m1

range tvar 0 20 100
predict h1, hazard zeros timevar(tvar)
predict s1, survival zeros timevar(tvar)

predictms , singleevent model(m1) hazard survival timevar(tvar)
