//local drive Z:/
local drive /Users/Michael
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 98798
clear
set obs 1000
gen id 	= _n
gen trt = runiform()>0.5
gen sd1 = exp(log(0.1))
gen u1 	= rnormal(0,sd1)
gen age = rnormal()
expand 10
sort id 
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 age 0.1 u1 1) maxt(5)
stset stime, f(died)

timer clear
timer on 1
mestreg trt age || id:, dist(weib) //evaltype(gf0))
predict s1, surv marginal
predict s5, surv cond(ebmeans)
predict d1, density marginal
timer off 1
// timer list
// predict s1, surv cond(ebmeans)
// // predictnl s2 = predict(surv cond(ebmeans)), ci(s2_lci s2_uci)
// predict r1, reffects reses(ser1)

timer on 2
merlin (stime trt age M1[id]@1, family(weib, failure(died)))	///
		, 
timer off 2

predict d2, density marginal
predict s2, survival marginal
predict s3, survival fixedonly
predict s6, surv fitted
// predict r2, reffects
// timer on 3
// merlin (stime trt age M1[id]@1, family(weib, failure(died)))	///
// 		, evaltype(gf1debug) //grad 
// timer off 3
// timer list
