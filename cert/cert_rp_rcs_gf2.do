
set seed 725498

clear
set obs 1000
gen trt = runiform()>0.5

survsim stime died, dist(weib) lambda(0.1) gamma(2) cov(trt -0.5) maxt(5)
qui stset stime, f(died)

gen bhaz = 0.1*runiform()

foreach dist in hazard loghazard rp {

	di "model: `dist'"
	
	forvalues df = 1/3 {
		di "df = `df'"
		local df df(`df')
		local baseline
		local tvar
		if "`dist'"!="rp" {
			local baseline rcs(stime, `df' log orthog)
			local tvar timevar(stime)
			local df
		}
		
		merlin (stime trt `baseline', family(`dist', `df' failure(died)) `tvar')  , evaltype(gf0)
		mat b1 = e(b)
		mat v1 = e(V)
		local b1 : display %9.5f b1[1,1]
		local v1 : display %9.5f v1[1,1]
		
		merlin (stime trt `baseline' , family(`dist', `df' failure(died)) `tvar')  , evaltype(gf1)
		mat b2 = e(b)
		mat v2 = e(V)
		local b2 : display %9.5f b2[1,1]
		local v2 : display %9.5f v2[1,1]
		
		merlin (stime trt `baseline' , family(`dist', `df' failure(died)) `tvar')  , evaltype(gf2)
		mat b3 = e(b)
		mat v3 = e(V)
		local b3 : display %9.5f b3[1,1]
		local v3 : display %9.5f v3[1,1]
		
		assert `b3'==`b2'
		assert `v3'==`v2'
	
		//tvc
		merlin (stime trt trt#fp(stime, power(0)) `baseline', timevar(stime) family(`dist', `df' failure(died)))  , evaltype(gf0)
		mat b1 = e(b)
		mat v1 = e(V)
		local b1 : display %9.5f b1[1,1]
		local v1 : display %9.5f v1[1,1]
		
		merlin (stime trt trt#fp(stime, power(0)) `baseline', timevar(stime) family(`dist', `df' failure(died)))  , evaltype(gf1)
		mat b2 = e(b)
		mat v2 = e(V)
		local b2 : display %9.5f b2[1,1]
		local v2 : display %9.5f v2[1,1]
		
		merlin (stime trt trt#fp(stime, power(0)) `baseline', timevar(stime) family(`dist', `df' failure(died)))  , evaltype(gf2)
		mat b3 = e(b)
		mat v3 = e(V)
		local b3 : display %9.5f b3[1,1]
		local v3 : display %9.5f v3[1,1]
		
		assert `b3'==`b2'
		assert `v3'==`v2'
		
	}
	
}
