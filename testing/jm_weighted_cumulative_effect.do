//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/megenreg/Stata"
adopath ++ "./megenreg"
clear all

do ./build/buildmlib.do
mata mata clear

use https://www.mjcrowther.co.uk/data/jm_example.dta, clear
gen trtage = trt*age
timer clear
timer on 1
/*
megenreg 	(logb fp(1)@l1 fp(1)#M2[id] M1[id], family(gaussian) timevar(time)) 			///
			(stime trt age trtage iEV[logb]@alpha, family(weib,  failure(died)))	 	///
			, intpoints(3) intmethod(gh) cov(unstr)

timer off 1
mat b = e(b)
mat b = b,0
*/
//weights cumulative - needs better starting values
mata:			
real matrix userhaz(gml, t)
{
	real matrix		linpred, logh, gq, nodes, weights, res
	real scalar 	gam1, alpha, sd, maxt, nq

	linpred = megenreg_util_xzb(gml,t)
	gam1 	= exp(megenreg_util_xb(gml,1))
	alpha 	= megenreg_util_xb(gml,2)
	sd		= exp(megenreg_util_xb(gml,3))
	sd2		= sd^2
	
	logh 	= linpred :+ log(gam1) :+ (gam1:-1):*log(t) 

	nq 		= 15
	gq 		= megenreg_gq(nq,"legendre")
	nt		= rows(t)
	
	nodes 	= t :* J(nt,1,gq[,1]'):/2 :+ t:/2
	weights = t :* J(nt,1,gq[,2]'):/2
	
	maxt 	= max(megenreg_util_timevar(gml))
	mvnmaxt = mvnormalcv(0,maxt,0,sd2)
	
	res 	= J(nt,1,0)
	
	for (q=1;q<=nq;q++) {
		res = res :+ weights[,q] :* normalden(t:-nodes[,q],sd) :* megenreg_util_xzb_mod(gml,1,nodes[,q]) :/ mvnmaxt
	}

	logh = logh :+ alpha :* res
	return(logh)
}
end
timer on 2

megenreg 	(logb fp(1)@l1 fp(1)#M2[id] M1[id], family(gaussian) timevar(time)) 								///
			(stime trt age trtage, family(user, loghfunc(userhaz) nap(3)  failure(died)) timevar(stime)) 	///
			, //intpoints(3) intmethod(gh) //cov(unstr) 
			
mat sv = e(b)
megenreg 	(logb fp(1)@l1 fp(1)#M2[id] M1[id], family(gaussian) timevar(time)) 								///
			(stime trt age trtage, family(user, hfunc(userhaz) nap(3)  failure(died)) timevar(stime)) 	///
			, from(sv) cov(unstr) 


timer off 2
timer list

