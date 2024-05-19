//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"

clear all
tr:do ./build/buildmlib.do
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
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) ///
	cov(trt -0.5 age 0.1 u1 1) maxt(5)
stset stime, f(died)

mestreg trt age || id:, dist(weib) 
predict r1, reffects reses(ser1)

merlin (stime trt age M1[id]@1, family(weib, failure(died)))
predict r2, reffects debug 
predict ser2, reses debug 


list r1 r2 if inlist(_n,1,11)
