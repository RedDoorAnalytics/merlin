//local drive Z:/
local drive /Users/Michael
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
clear all

do ./build/buildmlib.do
mata mata clear

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


survsim stime dead , hazard(`truebhaz' :* `lambda':*`gamma':*{t}:^(`gamma':-1):*exp(`beta':*trt)) maxt(5)
stset stime, f(dead)

merlin (stime trt, family(loglogistic, failure(dead)))
merlin (stime trt, family(loglogistic, failure(dead) bhazard(bhaz)))

stpm2 trt, scale(odds) df(1) noorthog bhazard(bhaz)
predict s1, surv
stmerlin trt, dist(loglogistic) bhazard(bhaz)
predict s2, surv
streg trt, dist(loglog)





