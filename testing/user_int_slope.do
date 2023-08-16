//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 500
gen id1 = _n
expand 5
bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor1 = (1,0.5\0.5,1)
drawnorm u1 u2, means(0 0) sds(0.5 0.5) corr(cor1)
//bys id1 (id2) : gen u1 = rnormal() if _n==1
bys id1 (id2) : replace u1 = u1[1]
//bys id1 (id2) : gen u2 = rnormal() if _n==1
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (-0.5+u2) * trt

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(trtui 1 u1 1) maxt(5)
stset stime, f(dead)

mata:
real matrix logl_function(transmorphic gml)	
{
	y1 		= merlin_util_depvar(gml)
	linpred = merlin_util_xzb(gml)
	gam 	= exp(merlin_util_xb(gml,1))

	//calculate and return the log-likelihood
	return(y1[,2] :* (linpred :+ log(gam) :+ (gam:-1):*log(y1[,1]) :+ log(y1[,1])) :- exp(linpred):*y1[,1]:^gam)
}

real matrix weib_haz(	transmorphic gml, 		///
						real matrix t)
{
	linpred = merlin_util_xzb(gml)
	gam 	= exp(merlin_util_ap(gml,1))

	//calculate and return the hazard function
	return(log(exp(linpred) :* gam :* t :^ (gam:-1)))
}

real matrix weib_chaz(	transmorphic gml, 		///
						real matrix t)	
{
	linpred = merlin_util_xzb(gml)
	gam 	= exp(merlin_util_xb(gml,1))

	//calculate and return the cumulative hazard function
	return(exp(linpred) :* t :^ gam)
}

real matrix poly_haz(	transmorphic gml, 		///
						real matrix t)	
{
	linpred = merlin_util_xzb(gml)
	gam1 = merlin_util_xb(gml,1)
	gam2 = merlin_util_xb(gml,2)

	//calculate and return the hazard function
	return(exp(linpred :+ gam1 :* t :+ gam2 :* t:^2))
}

real matrix rcs_haz(transmorphic gml, 		///
					real matrix t)	
{
	linpred 	= merlin_util_xzb(gml)
	xb 			= merlin_util_xb(gml,1)\merlin_util_xb(gml,2)
	knots 		= (-5,-1,1.609)

	//calculate and return the hazard function
	return(exp(linpred :+ merlin_rcs(log(t),knots) * xb))
}

real matrix rcs_haz2(transmorphic gml, 		///
					real matrix t)	
{
	linpred 	= merlin_util_xzb(gml)
	xb 			= merlin_util_xb(gml,1)\merlin_util_xb(gml,2)
	knots 		= (-5.08854703,.8450899120,1.6091114621)
	age = st_data(.,"age")
	
	return(exp(linpred :+ merlin_rcs(log(t:+age),knots) * xb))
}

end

// merlin (stime trt trt#M2[id1] M1[id1] , family(weib, failure(dead))), adaptopts(log)
// merlin (stime trt trt#M2[id1] M1[id1] , family(user, llfunction(logl_function) nap(1) failure(dead))),

merlin (stime trt M1[id1]@1 , family(user, loghfunction(weib_haz) nap(1) failure(dead)))

predict s1, eta marginal
range tvar 0 10 10
predict h2, timelost marginal timevar(tvar)

// merlin (stime trt trt#M2[id1] M1[id1] , family(user, hfunction(weib_haz) chfunction(weib_chaz) nap(1) failure(dead))), 
// merlin (stime trt trt#M2[id1] M1[id1] , family(user, hfunction(poly_haz) nap(2) failure(dead))),
// merlin (stime trt trt#M2[id1] M1[id1] , family(user, hfunction(rcs_haz) nap(2)  failure(dead)))
// gen age = 0
// merlin (stime trt trt#M2[id1] M1[id1] , family(user, hfunction(rcs_haz2) nap(2)  failure(dead))), 




