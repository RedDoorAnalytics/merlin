
//mestreg check

set seed 98798
clear
set obs 300
gen id 	= _n
gen trt = runiform()>0.5
gen sd1 = exp(log(0.5))
gen u1 	= rnormal(0,sd1) 
gen age = rnormal()
expand 10
sort id 
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) cov(trt -0.5 age 0.1 u1 1) maxt(5)
stset stime, f(died)

foreach dist in weibull {

	mestreg trt age || id:, dist(`dist')
	mat b1 = e(b)
	local b1 : display %9.4f b1[1,1]
	local b2 : display %9.4f b1[1,2]
	local v1 : display %9.4f b1[1,5]
	
	cap drop r1 ser1
	predict r1, reffects reses(ser1)

	cap drop s1
	predict s1, surv cond(ebmeans)
	
        cap drop s8
        predict s8, surv marginal
        
        cap drop s4
        predict s4, surv cond(fixedonly)
        
        cap drop d1
        predict d1, density marginal
        
        cap drop h1
        predict h1, haz marginal

	merlin (stime trt age M1[id]@1, family(`dist', failure(died)))
	mat b2 = e(b)
	local b3 : display %9.4f b2[1,1]
	local b4 : display %9.4f b2[1,2]
	local v2 =  exp(b2[1,6])^2
	local v2 : display %9.4f `v2'

	cap drop refs* serefs*
	predict refs*, reffects 
	predict serefs*, reses
	cap drop s2
	predict s2, surv fitted
        cap drop s9
        predict s9, surv marginal
        cap drop s5
        predict s5, surv fixedonly
        cap drop d2
        predict d2, density marginal
        cap drop h2
        predict h2, haz marginal
        
	//check
	assert `b1'==`b3'
	assert `b2'==`b4'
	assert `v1'==`v2'
	cap drop t1 t2 t3 t4 t5
	gen t1 = abs(r1-refs1)
	gen t2 = abs(ser1 - serefs1)
	gen t3 = abs(s1-s2)
        gen t4 = abs(s8-s9)
        gen t5 = abs(s4-s5)
        gen t6 = abs(d1-d2)
        gen t7 = abs(h1-h2)
	assert t1<=0.00001
	assert t2<=0.00001
	assert t3<=0.00001
        assert t4<=0.00001
        assert t5<=0.00001
        assert t6<=0.00001
        assert t7<=0.00001
	
}

// with ltrunc

set seed 98798
clear
set obs 300
gen id 	= _n
gen trt = runiform()>0.5
gen sd1 = exp(log(0.5))
gen u1 	= rnormal(0,sd1) 
gen age = rnormal()
expand 10
sort id 
gen t0 = runiform()
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
	cov(trt -0.5 age 0.1 u1 1) maxt(5) ltruncated(t0)
stset stime, f(died) enter(t0)

foreach dist in weibull {

	gsem _t <- trt age M1[id]@1 , family(`dist', failure(_d) ltruncated(_t0))
	mat b1 = e(b)
	local b1 : display %9.4f b1[1,1]
	local b2 : display %9.4f b1[1,2]
	local v1 : display %9.4f b1[1,5]
	
	cap drop s1
	predict s1, surv cond(ebmeans)
	
        cap drop s8
        predict s8, surv marginal
        
        cap drop s4
        predict s4, surv cond(fixedonly)
        
        cap drop d1
        predict d1, density marginal
        
        merlin (stime trt age M1[id]@1, family(`dist', failure(died) ltruncated(t0)))
	mat b2 = e(b)
	local b3 : display %9.4f b2[1,1]
	local b4 : display %9.4f b2[1,2]
	local v2 =  exp(b2[1,6])^2
	local v2 : display %9.4f `v2'

	predict s2, surv fitted
        cap drop s9
        predict s9, surv marginal
        cap drop s5
        predict s5, surv fixedonly
        cap drop d2
        predict d2, density marginal
        
	//check
	assert `b1'==`b3'
	assert `b2'==`b4'
	cap drop t1 t2 t3 t4 t5
	gen t3 = abs(s1-s2)
        gen t4 = abs(s8-s9)
        gen t5 = abs(s4-s5)
        gen t6 = abs(d1-d2)
        assert t3<=0.00001
        assert t4<=0.00001
        assert t5<=0.00001
        assert t6<=0.00001
        
	
}


// with marginal ltrunc

set seed 98798
clear
set obs 3000
gen id 	= _n
gen trt = runiform()>0.5
gen u1 	= rnormal(0,0.5)
gen age = rnormal()
sort id 
gen t0 = 2*runiform()
survsim stime died , dist(weib) lambda(0.1) gamma(1.2) 	///
	cov(trt -0.5 age 0.1 u1 1) maxt(10) ltruncated(t0)
stset stime, f(died) enter(t0)

     
merlin (stime trt age M1[id]@1, family(weibull, failure(died) ///
	ltruncated(t0, marginal))), intpoints(11)

assert reldif( e(ll)          , -6149.969100965593) <  1E-8
	
qui {
mat T_b = J(1,6,0)
mat T_b[1,1] =  -.578985451169165
mat T_b[1,2] =  .0801909076614118
mat T_b[1,3] =                  1
mat T_b[1,4] = -2.209600641461436
mat T_b[1,5] =  .1747601315470237
mat T_b[1,6] = -.5554496156848963
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-6
