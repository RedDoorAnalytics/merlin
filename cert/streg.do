
set seed 725498

clear
set obs 500
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.1) gamma(1.5) cov(trt -0.5) maxt(5)
qui stset stime, f(died)

foreach dist in exp weib gomp loglogistic lognormal ggamma {
	streg trt, dist(`dist')
	mat b1 = e(b)
	local b1 : display %9.5f b1[1,1]
	
	//merlin
	merlin (stime trt , family(`dist', failure(died)))  
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]

	assert `b1'==`b2'
	
	do ./cert/predictions.do
	di "stmerlin, dist(`dist')"
	//stmerlin
	stmerlin trt , dist(`dist')
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	
	assert `b1'==`b2'
	
	do ./cert/predictions.do
}


set seed 725498

clear
set obs 500
gen trt = runiform()>0.5
gen lt0 = runiform()*4
survsim stime died, dist(weib) lambda(0.1) gamma(1.5) cov(trt -0.5) maxt(10) ///
        ltruncated(lt0)
qui stset stime, f(died) enter(lt0)

foreach dist in exp weib gomp loglogistic lognormal ggamma {
	streg trt, dist(`dist')
	mat b1 = e(b)
	local b1 : display %9.5f b1[1,1]
	
	//merlin
	merlin (stime trt , family(`dist', failure(died) ltruncated(lt0)))  
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]

	assert `b1'==`b2'
	
	do ./cert/predictions.do
	di "stmerlin, dist(`dist')"
	//stmerlin
	stmerlin trt , dist(`dist')
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	
	assert `b1'==`b2'
	
	do ./cert/predictions.do
}

