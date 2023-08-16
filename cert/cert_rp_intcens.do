set seed 72549

clear

set obs 1000
gen trt = runiform()>0.5
survsim stime event , dist(weib) lambda(0.1) gamma(1.2) covariates(trt -0.5) maxt(5)

gen st1 = floor(stime) if stime < 5
replace stime = st1 + 1 if st1 < 5

// merlin identifies interval censored oberservations through the `failure()` indicator being coded as a 2. 
replace event = 2 if event == 1

foreach dist in rp {
	di "model: `dist'"
	forvalues df = 1/3 {
		di "df = `df'"
		
		merlin (stime trt `baseline', family(`dist', df(`df') failure(event) linterval(st1)))  , evaltype(gf0)
		mat b1 = e(b)
		mat v1 = e(V)
		local b1 : display %9.5f b1[1,1]
		local v1 : display %9.5f v1[1,1]
		
		merlin (stime trt `baseline' , family(`dist', df(`df') failure(event) linterval(st1)))  , evaltype(gf1)
		mat b2 = e(b)
		mat v2 = e(V)
		local b2 : display %9.5f b2[1,1]
		local v2 : display %9.5f v2[1,1]
		
		merlin (stime trt `baseline' , family(`dist', df(`df') failure(event) linterval(st1))) 
		assert "`e(ml_method)'"=="gf2"
		mat b3 = e(b)
		mat v3 = e(V)
		local b3 : display %9.5f b3[1,1]
		local v3 : display %9.5f v3[1,1]
		
		assert `b3'==`b2'
		assert `v3'==`v2'
		
		//stpm bench
		capture {
			replace st1 = stime if st1==.
			gen died = event>0
		}
		stset stime, f(died)
		timer on 2
		stpm trt, scale(h) df(`df') left(st1)
		mat b4 = e(b)
		mat v4 = e(V)
		
		local b4 : display %9.5f b4[1,`=`df'+1']
		local v4 : display %9.5f v4[`=`df'+1',`=`df'+1']
		
		assert `b4'==`b3'
		assert `v4'==`v3'
		
		
		
	}
	
}

