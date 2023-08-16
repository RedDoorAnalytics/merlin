//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

do ./build/buildmlib.do
mata mata clear

use https://www.mjcrowther.co.uk/data/jm_example.dta, clear

timer clear
timer on 1
merlin 	(stime trt 												///
				EV[logb]										///
				, timevar(stime) family(weib, failure(died))) 	///
		(logb 	fp(time,pow(1))									///
				fp(time,pow(1))#M2[id]@1						///
						M1[id]@1								///
				, timevar(time) family(gaussian)) 				///
		, trace
timer off 1
timer list


gen tp = round(prothrombin)
merlin 	(tp time , family(nbinomial) timevar(time))  	///
		, random

merlin 	(logb rcs(time,df(3) orthog) , family(gaussian) timevar(time))  	///
		, zeros 
merlin 	(stime trt, family(weibull, failure(died)))  	///
		, zeros
merlin 	(logb time M1[id]@1, family(gaussian))  	///
		, inpoints(7) intmethod(gh) zeros hessian
merlin 	(logb time time#M2[id]@1 M1[id]@1, family(gaussian))  	///
		, inpoints(35) intmethod(gh) zeros 

merlin 	(logb time M1[id]@1 time#M2[id]@1, family(gaussian) timevar(time))  	///
		(stime trt XB[logb], family(weibull, failure(died)) timevar(stime))							///
		, inpoints(7) intmethod(gh) cov(diag)		zeros
		
merlin 	(stime trt, family(rp, df(3) failure(died)))  	///
		, 

merlin 	(stime trt trt#rcs(stime,df(1) log event orthog), family(rp, df(3) failure(died)) timevar(stime))  	///
		, 
		
mata:
real matrix lev1_logl(gml)
{
     y        	= merlin_util_depvar(gml)               //response
     xb 		= merlin_util_xzb(gml)                  //main lin. pred.
     residxb 	= exp(merlin_util_xzb_mod(gml,2))       //lev1 lin. pred
	 logl		= lnnormalden(y,xb,sqrt(residxb))		//logl
     return(logl)
}
end

timer clear
timer on 1
merlin 	(logb time#M2[id]@1 time M1[id]@1, family(user, llfunction(lev1_logl)))  	///
		(M3[id]@1, family(null))													///
		(stime trt M1[id] M2[id] M3[id], family(weibull, failure(died)))							///
		, inpoints(15) intmethod(gh) cov(unstr)
timer off 1
timer list
//5mins
timer clear
timer on 1
merlin 	(stime trt 												///
				EV[logb]												///
				, timevar(stime) family(weib, failure(died))) 	///
		(logb 	fp(time,pow(1))									///
				fp(time,pow(1))#M2[id]@1							///
													///
				, timevar(time) family(gaussian)) 					///
		, intmethod(gh) intpoints(55)
timer off 1
timer list
cap drop tvar
range tvar 0 10 100

//touse needs work in posting

predict m1, mu outcome(2)

// predict s1, survival fixedonly  at(trt 1)
// predict s2, survival marginal at(trt 1)
predict rmst1, rmst fixedonly timevar(tvar) //at(trt 1)
predict rmst2, rmst marginal timevar(tvar)

//predictnl rmst3 = predict(rmst marginal), ci(rmst2_lci rmst2_uci)
// scatter s1 s2 stime2
