//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"

clear all
tr:do ./build/buildmlib.do
mata mata clear

set seed 7254
clear 

set obs 100
gen id1 = _n
gen age = runiform()>0.5
gen sd1 = exp(log(1)+0*age)
gen u1 = rnormal(0,sd1)
expand 5
bys id1: gen id2 = _n
gen trt = runiform()>0.5
gen sd2 = exp(log(1)+0*trt)
gen u2 = rnormal(0,sd2)
expand 5
gen id3 = _n
sort id1 id2 id3

survsim stime1 dead1 , dist(weib) lambda(0.1) gamma(1.2) cov(age -0.5 u1 1 u2 1) maxt(5)
stset stime1, f(dead1)

replace id2 = _n
// gsem (stime1 <- age M2[id1>id2]@1 M1[id1]@1 , family(weib, failure(dead1))) , //intmethod(gh)
// predict gs1, survival fixedonly 
// predict gs2, survival marginal 
			
merlin 	(stime1 age M2[id1>id2]@1 M1[id1]@1 , family(weib, failure(dead1))) ///
			, intmethod(gh) intpoints(7)

merlin 	(stime1 age M2[id1>id2]@1 M1[id1]@1 , family(cox, failure(dead1))) ///
			, intmethod(gh) intpoints(7) devcode5(294820)
			
range tvar 0 10 100

pr drop _all
predict s1, survival fixedonly //ci //timevar(tvar)
predict s2, survival marginal ci //timevar(tvar)
predict s3, hazard fixedonly timevar(tvar)
predict s4, hazard marginal timevar(tvar)
predict s5, rmst fixedonly ci //timevar(tvar)
predict s6, rmst marginal ci //timevar(tvar)
predict s7, chazard fixedonly timevar(tvar)
predict s8, chazard marginal timevar(tvar)
predict s9, eta fixedonly timevar(tvar)
predict s10, eta marginal timevar(tvar)
predict s11, mu fixedonly timevar(tvar)
predict s12, mu marginal timevar(tvar)

predictnl s13 = predict(rmst marginal), ci(s13_lci s13_uci)
			
