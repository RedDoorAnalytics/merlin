
set seed 725498

clear
set obs 1000
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.01) gamma(1.2) cov(trt -0.5) maxt(5)
qui stset stime, f(died)

stcox trt, nohr 
mat b1 = e(b)
local b1 : display %9.5f b1[1,1]
merlin (stime trt , family(cox, failure(died)))  
mat b2 = e(b)
local b2 : display %9.5f b2[1,1]

assert `b1'==`b2'

predict s1, survival
predict s2, survival ci 
predict h1, hazard
predict s3, survival standardise

drop stime died s1 s2 h1 s3
gen lt0 = 3*runiform()
survsim stime died, dist(weib) lambda(0.01) gamma(1.2) cov(trt -0.5) maxt(10) ltruncated(lt0)
qui stset stime, f(died) enter(lt0)

stcox trt, nohr 
mat b1 = e(b)
local b1 : display %9.5f b1[1,1]
stmerlin trt , dist(cox)
mat b2 = e(b)
local b2 : display %9.5f b2[1,1]

assert `b1'==`b2'

predict s1, survival
// predict s2, survival ci 
predict h1, hazard
predict s3, survival standardise



//with ties

set seed 725498

clear
set obs 1000
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.01) gamma(1.2) cov(trt -0.5) maxt(5)
replace stime = ceil(stime)
qui stset stime, f(died)

stcox trt, nohr 
mat b1 = e(b)
local b1 : display %9.5f b1[1,1]
merlin (stime trt , family(cox, failure(died)))  
mat b2 = e(b)
local b2 : display %9.5f b2[1,1]

assert `b1'==`b2'

predict s1, survival
predict s2, survival ci 
predict h1, hazard
predict s3, survival standardise

drop stime died s1 s2 h1 s3
gen lt0 = 3*runiform()
survsim stime died, dist(weib) lambda(0.01) gamma(1.2) cov(trt -0.5) maxt(10) ltruncated(lt0)
replace stime = ceil(stime)
qui stset stime, f(died) enter(lt0)

stcox trt, nohr 
mat b1 = e(b)
local b1 : display %9.5f b1[1,1]
stmerlin trt , dist(cox)
mat b2 = e(b)
local b2 : display %9.5f b2[1,1]

assert `b1'==`b2'

predict s1, survival
// predict s2, survival ci 
predict h1, hazard
predict s3, survival standardise
