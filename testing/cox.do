//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"


//build mlib
clear all
do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 300

gen trt = runiform()>0.5
gen age = rnormal(50,5)
gen lt0 = 0 //runiform()*2
survsim stime died, dist(weib) lambda(0.01) gamma(2) cov(trt -0.5 age 0.01) maxt(5) ltruncated(lt0)
// survsim stime died, hazard(0.1:*1.2:*#t:^0.2 :* exp(0.1:*log(#t):*trt)) cov(trt -0.5) maxt(5)
//survsim stime died, cumhazard(0.1:*#t:^1.2 :* exp(0.1:*#t:*trt)) cov(trt -0.5 u1 1) maxt(5)

replace stime = ceil(stime)
qui stset stime, f(died) enter(lt0)

// timer clear
// timer on 1
stcox trt age, nohr 
// predict h1, basehc
// predict ch1, basechazard
// timer off 1

// timer on 2
merlin (stime trt age if _n>5, family(cox, failure(died) ))  , //evaltype(gf0)
// merlin (stime trt age, family(cox, failure(died) ))  , evaltype(gf0)
// timer off 2

// predict h2, hazard
// predict ch2, chazard

predict f1, surv 
predict f2, surv at1(trt 1) at2(trt 0) ci standardise reps(10)
predict f3, surv at1(trt 1) at2(trt 0) ci 

// predict h1, basehazard 

// predict h2, hazard at(trt 1)
// predict h3, hazard zeros ci

// predict s1, surv ci

// timer on 3
// merlin (stime trt age, family(cox, failure(died)))  , evaltype(gf1)
// timer off 3

// timer on 4
// merlin (stime trt age, family(cox, failure(died)))  , evaltype(gf2)
// timer off 4


// timer list

// predict e1, eta ci
// replace e1 = exp(e1)
// predict hr1, hratio at1(trt 1) at2(trt 0) ci
		
// timer clear
// timer on 1

