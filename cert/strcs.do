
set seed 725498

clear
set obs 300
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.01) gamma(1.5) cov(trt -0.5) maxt(5)
qui stset stime, f(died)


	strcs trt, df(3) nohr
	mat b1 = e(b)
	local b1 : display %9.5f b1[1,1]
	local s1 : display %9.3f b1[1,2]
	
	//merlin
	merlin (stime trt rcs(stime, df(3) log orthog event) , family(loghazard, failure(died)) timevar(stime))  
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.3f b2[1,2]
	assert `s1'==`s2'
	
	do ./cert/predictions.do
	
	//stmerlin
	stmerlin trt , dist(rcs) df(3)
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.3f b2[1,2]
	assert `b1'==`b2'
	assert `s1'==`s2'
	
	do ./cert/predictions.do

//noorthog
	
	strcs trt, df(3) noorthog nohr
	mat b1 = e(b)
	local b1 : display %9.5f b1[1,1]
	local s1 : display %9.2f b1[1,2]
	
	//merlin
	merlin (stime trt rcs(stime, df(3) log event) , family(loghazard, failure(died)) timevar(stime))  
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.2f b2[1,2]
	assert `s1'==`s2'
	
	do ./cert/predictions.do
	
	//stmerlin
	stmerlin trt , dist(rcs) df(3) noorthog
	mat b2 = e(b)
	local b2 : display %9.5f b2[1,1]
	local s2 : display %9.2f b2[1,2]
	assert `b1'==`b2'
	assert `s1'==`s2'
	
	do ./cert/predictions.do
