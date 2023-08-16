set seed 72549

clear
set obs 1000

//exp mortailty
local truebhaz = 0.1
gen bhaz = exp(rnormal(log(`truebhaz'),0.5))
gen lbhaz = log(bhaz)

//excess mortality
local lambda = 0.1
local gamma  = 1.2
local beta = -0.5
gen trt = rbinomial(1,0.5)


survsim stime dead , hazard(`truebhaz' :+ `lambda':*`gamma':*{t}:^(`gamma':-1):*exp(`beta':*trt)) maxt(5)
stset stime, f(dead)

//merlin (stime trt, family(weibull, failure(dead) bhazard(bhaz)))
mata mata clear

mata:
real matrix logl1(gml, | G, H)
{
	y1 		= merlin_util_depvar(gml)
	linpred = merlin_util_xzb(gml)
	gam 	= exp(merlin_util_ap(gml,1))
	bhaz 	= exp(merlin_util_xzb_mod(gml,2))
	return(y1[,2] :* log(bhaz :+ exp(linpred :+ log(gam) :+ (gam:-1):*log(y1[,1]))) :- exp(linpred):*y1[,1]:^gam)
}
end

merlin 	(stime trt, family(user, llfunction(logl1) failure(dead) nap(1))) ///
		(lbhaz , family(gaussian))
est store m1
do ./cert/predictions.do

merlin (stime trt, family(weibull, failure(dead) bhazard(bhaz)))
est store m2
do ./cert/predictions.do






