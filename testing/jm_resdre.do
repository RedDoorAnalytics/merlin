
set seed 7254
pr drop _all
clear
set obs 1000
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal(0,0.5)
gen trt = runiform()>0.5

//simulate form weibull linking both random intercepts
survsim stime died, hazard(0.1:*1.2:*#t:^0.2 :* 	///
	exp(0.2 :* u1  :- 0.2 * u2 )) cov(trt -0.5) maxt(5)
expand 5
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u1 + 0.1*time + 0.5*trt		//lin pred for continuous response
gen sd = exp(log(0.5)+u2)				//exp of lin pred for log residual standard deviation
gen y = rnormal(xb,sd)

bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1

//logL funciton of gaussian response, with complex residual variance
mata:
real matrix logl_function(transmorphic gml)
{
	y1 		= megenreg_util_depvar(gml)
	expval 	= megenreg_util_xzb(gml)
	sd 		= exp(megenreg_util_xzb_mod(gml,3))
	return(lnnormalden(y1, expval, sd))
}
end

//fit megenreg model, combining Weibull, user function for gaussian, and extra null model for log of residual standard dev.
megenreg 	(stime trt M1[id1]@alp1 M2[id1]@alp2 , 	family(weib, failure(died))) 							///
			(y trt fp(1)@lam1 M1[id1] , 			family(user, llfunction(logl_function)) timevar(time)) 	///
			(M2[id1] , 								family(null))  											///
			, intmethod(gh) 


