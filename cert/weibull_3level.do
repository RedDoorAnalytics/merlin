
clear
set seed 498093
set obs 15
gen trial = _n
gen u1 = rnormal(0,1)

expand 100
sort trial
gen id = _n
gen u2 = rnormal(0,1)
gen trt = runiform()>0.5

expand 2
sort trial id

survsim stime died, dist(weibull) lambda(0.1) gamma(1.2) maxtime(5) 	///
					covariates(trt -0.5 u1 1 u2 1)
					
merlin (stime trt M1[trial]@1 M2[trial>id]@1, family(weibull, failure(died))) , //redist(normal t) df(3) intmethod(mc)
do ./cert/predictions.do

clear
set seed 498093
set obs 15
gen trial = _n
gen b11 = rnormal(0,1)
gen b12 = rnormal(0,1)
expand 100
sort trial
gen id = _n
gen b2 = rnormal(0,1)
gen trt = runiform()>0.5
expand 2
sort trial id
gen trtb12 = (-0.5 + b12) * trt
survsim stime died, dist(weibull) lambda(0.1) gamma(1.2) maxtime(5) 	///
			covariates(trtb12 1 b11 1 b2 1)
merlin (stime trt trt#M1[trial]@1 M2[trial]@1 M3[trial>id]@1, family(weibull, failure(died))) 
         
do ./cert/predictions.do
