
set seed 725498

clear
set obs 1000
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.01) gamma(1.5) cov(trt -0.5) maxt(5)
qui stset stime, f(died)


	stpm2 trt, scale(h) df(3)
	mat b1 = e(b)
	local b1 : display %9.5f b1[1,1]
	local s1 : display %9.5f b1[1,2]
	
	//merlin
	merlin (stime trt , family(rp, failure(died) df(3)))  
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.5f b2[1,3]

	assert `s1'==`s2'
	
	do ./cert/predictions.do
	
	//stmerlin
	stmerlin trt , dist(rp) df(3)
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.5f b2[1,3]
	assert `b1'==`b2'
	assert `s1'==`s2'
	
	do ./cert/predictions.do


