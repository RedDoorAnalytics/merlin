//local drive Z:/
local drive /Users/Michael/My Drive/software/merlin
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

survsim stime died , 	dist(weib) lambda(0.1) gamma(1.2) 				///
						cov(trt -0.5 age 0.01 bmi -0.05) maxt(10) ltruncated(t0)	

stset stime, f(died) enter(t0)

timer clear
timer on 1
// streg trt age bmi, dist(weib) nohr
// stpm2 trt age bmi , scale(h) df(3) tvc(trt) dftvc(1) 
// strcs trt age bmi, nohr df(3) //tvc(trt) dftvc(1) 
// stcox trt age bmi, nohr tvc(trt) texp(log(_t)) 
// stmerlin trt age bmi, dist(cox) 
timer off 1

timer on 2
merlin (_t 	trt rcs(age,df(3)) bmi 			///
		if _n<100, family(weib, failure(died))) 	///
                , eform
timer off 2
timer list


predict h1, hazard
predict s1, surv standardise at(trt 1)


// mata:
// real colvector predf(gml)
// {	
// 	t = merlin_util_timevar(gml)
// 	gam = merlin_util_dap(gml,1)
// 	return(exp(merlin_util_xzb(gml,t)):*gam:*t:^(gam:-1))
// }
// end

// predict h1, userfunction(predf) at(trt 1) zeros ci

// predict h2, hazard at(trt 1) zeros ci

