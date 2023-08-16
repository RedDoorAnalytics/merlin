set seed 72549

clear
set obs 10000

//exp mortailty
local truebhaz = 0.01
gen bhaz = exp(rnormal(log(`truebhaz'),0.5))
gen lbhaz = log(bhaz)

//excess mortality
local lambda = 0.1
local gamma  = 1.2
local beta = -0.5
gen trt = rbinomial(1,0.5)


survsim stime dead , hazard(bhaz :+ `lambda':*`gamma':*{t}:^(`gamma':-1):*exp(`beta':*trt)) maxt(5)
stset stime, f(dead)

stpm2 trt, scale(haz) df(3)  bhazard(bhaz)
predict s1, surv
stmerlin trt, dist(rp) df(3) bhazard(bhaz)
predict s2, surv

assert reldif(s1,s2)<1e-05


stpm2 trt, scale(haz) df(3)  bhazard(bhaz) tvc(trt) dftvc(2)
predict s3, surv
stmerlin trt, dist(rp) df(3) bhazard(bhaz) tvc(trt) dftvc(2)
predict s4, surv

assert reldif(s3,s4)<1e-05








