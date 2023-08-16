
set seed 725498

clear
set obs 1000
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.1) gamma(1.5) cov(trt -0.5) maxt(5)
qui stset stime, f(died)

foreach dist in cox exp weib gomp /*loglogistic lognormal ggamma*/ {
	
	merlin (stime trt , family(`dist', failure(died)))  , evaltype(gf0)
	mat b1 = e(b)
	mat v1 = e(V)
	local b1 : display %9.5f b1[1,1]
	local v1 : display %9.5f v1[1,1]
	
	merlin (stime trt , family(`dist', failure(died)))  , evaltype(gf1)
	mat b2 = e(b)
	mat v2 = e(V)
	local b2 : display %9.5f b2[1,1]
	local v2 : display %9.5f v2[1,1]
	
	merlin (stime trt , family(`dist', failure(died)))  , evaltype(gf2)
	mat b3 = e(b)
	mat v3 = e(V)
	local b3 : display %9.5f b3[1,1]
	local v3 : display %9.5f v3[1,1]
	
	assert `b3'==`b2'
	assert `v3'==`v2'
	
	//tvc
	merlin (stime trt trt#fp(stime, power(0)), timevar(stime) family(`dist', failure(died)))  , evaltype(gf0)
	mat b1 = e(b)
	mat v1 = e(V)
	local b1 : display %9.5f b1[1,1]
	local v1 : display %9.5f v1[1,1]
	
	merlin (stime trt trt#fp(stime, power(0)), timevar(stime) family(`dist', failure(died)))  , evaltype(gf1)
	mat b2 = e(b)
	mat v2 = e(V)
	local b2 : display %9.5f b2[1,1]
	local v2 : display %9.5f v2[1,1]
	
	merlin (stime trt trt#fp(stime, power(0)), timevar(stime) family(`dist', failure(died)))  , evaltype(gf2)
	mat b3 = e(b)
	mat v3 = e(V)
	local b3 : display %9.5f b3[1,1]
	local v3 : display %9.5f v3[1,1]
	
	assert `b3'==`b2'
	assert `v3'==`v2'
	
}

