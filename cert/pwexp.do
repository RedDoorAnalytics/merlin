//=============================================================================================================================//
// cert. script for merlin, family(pwexp)



set seed 725498

clear
set obs 1000
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.01) gamma(1.5) cov(trt -0.5) maxt(5)
gen id = _n
qui stset stime, f(died) id(id)

	preserve
	stsplit sp, at(1 2 3)
	streg trt ibn.sp, dist(exp) nohr noconst
	
	mat b1 = e(b)
	local b1 : display %9.5f b1[1,1]
	local s1 : display %9.3f b1[1,2]
	restore 
	
	//merlin
	merlin (stime trt , family(pwexp, knots(1 2 3) failure(died)))  
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.3f b2[1,2]
	assert `s1'==`s2'
	
	do ./cert/predictions.do
	
	//stmerlin
	stmerlin trt , dist(pwexp) knots(1 2 3)
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.3f b2[1,2]
	assert `b1'==`b2'
	assert `s1'==`s2'
	
	do ./cert/predictions.do
