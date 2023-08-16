*! version 1.0.0 ?????2016

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
local TR 		transmorphic
local RS 		real scalar
local RC 		real colvector
local SS 		string scalar
local PS 		pointer scalar
local RR 		real rowvector
local RM 		real matrix
local PC 		pointer colvector
local PM 		pointer matrix
local SC 		string colvector

version 14.1

mata:

`RM' merlin_logl_aft(	`gml' gml, 	///
						`RM' G, 	///
						`RM' H)
{
	model 		= gml.model
	y 			= merlin_util_depvar(gml)
	Nobs		= merlin_get_nobs(gml,model)
	hast		= gml.istimedep[model,1]
	haslt		= gml.hasltrunc[model]
	logl 		= J(Nobs,1,0)
	
	//core
	xb			= merlin_util_xzb(gml)

	//===========================================================================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 			= merlin_get_surv_index(gml)
			logl[index1,] 	= xb[index1,] :+ log(y[index1,1])
		}

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		Nobs2 		= merlin_get_nobs(gml)
		if (Nobs2) {
			index2 	= merlin_get_surv_index(gml)
			if (hast) {
				Ngq 	= gml.chip
				chq2 	= J(Nobs2,Ngq,0)
				gq 		= merlin_gq(Ngq,"legendre")
				if (!haslt) qp2	= y[index2,1] :/ 2 :* J(Nobs2,1,gq[,1]') :+ y[index2,1]:/2
				else 		qp2	= (y[index2,1]:-y[index2,3]) :/ 2 :* J(Nobs2,1,gq[,1]') :+ (y[index2,1]:+y[index2,3]):/2
				for (q=1;q<=Ngq;q++) {
					chq2[,q] = exp(merlin_util_xzb(gml,qp2[,q]))
				}
				if (!haslt) logl[index2] = logl[index2] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
				else 		logl[index2] = logl[index2] :- (y[index2,1]:-y[index2,3]):/2 :* (chq2 * gq[,2]) 
			}
			else {
				logl[index2] = logl[index2] :- exp(xb[index2,]) :* y[index2,1]
				if (haslt) logl[index2] = logl[index2] :+ exp(xb[index2,]) :* y[index2,3]
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		Nobs3 	= merlin_get_nobs(gml)
		if (Nobs3) {
			//exit times
			index3 = merlin_get_surv_index(gml)
			if (hast) {
				Ngq 	= gml.chip
				chq3 	= J(Nobs3,Ngq,.)
				gq 		= merlin_gq(Ngq,"legendre")
				qp3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,1]') :+ y[index3,1]:/2
				for (q=1;q<=Ngq;q++) {
					chq3[,q] = exp(merlin_util_xzb(gml,qp3[,q]))
				}
				ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
				logl[index3,] = 1:-exp(-ch3)
			}
			else logl[index3,] = -exp(-exp(xb[index3,]) :* y[index3,1])
			
			//entry times
			gml.survind = 5
			if (gml.hasltrunc[model]) 	tind = 4
			else 						tind = 3
			index5 	= merlin_get_surv_index(gml)
			if (hast) {
				Nobs5 	= merlin_get_nobs(gml)
				chq5 	= J(Nobs5,Ngq,.)
				qp5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,1]') :+ y[index5,tind]:/2
				for (q=1;q<=Ngq;q++) {
					chq5[,q] = exp(merlin_util_xzb(gml,qp5[,q]))
				}
				ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
				logl[index5,] = logl[index5,] :- (1:-exp(-ch5))
			}
			else logl[index5,] = logl[index5,] :+ exp(-exp(xb[index5,]) :* y[index5,tind])
			//logL
			logl[index3,] = log(logl[index3,])
		}

	if (gml.todo==0) return(logl)
	
	
}

`RM' merlin_aft_h(`gml' gml, `RC' t)
{	
	smp  			= c("epsdouble")
	hstep 			= J(rows(t),1,1)
	index 			= selectindex(abs(t):<=1)
	hstep[index] 	= abs(t)[index]
	hstep 			= hstep :* smp :^(1/3)
	lh 				= merlin_aft_ch(gml,t :+ hstep)
	rh 				= merlin_aft_ch(gml,t :- hstep)
	return((lh :- rh):/(2:*hstep))
}

`RM' merlin_aft_logh(`gml' gml, `RC' t)
{	
	haz = merlin_aft_h(gml,t)
	if (gml.hasbh[gml.model,1]) {
		haz = haz :+ asarray(gml.bhazards,gml.model)[asarray(gml.surv_index,(gml.model,1))]
	}
	return(log(haz))
}

`RM' merlin_aft_ch(`gml' gml,`RC' t)
{
	mod		= gml.model
	xb		= -merlin_util_xzb(gml,t)
	basis 	= log(t:*exp(xb))
	nc		= cols(basis)

	brcs		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	if (hasorthog) rmat = asarray(gml.distancb,(mod,5))

	xbs = dxbs = J(merlin_get_nobs(gml),nc,.)
	if (hasorthog) { 
		for (i=1;i<=nc;i++) {
			xbs[,i]  = merlin_rcs(basis[,i],knots,0,rmat) * brcs
// 			dxbs[,i] = merlin_rcs(basis[,i],knots,1,rmat) * brcs
		}
	}
	else {
		for (i=1;i<=nc;i++) {
			xbs[,i]  = merlin_rcs(basis[,i],knots,0) * brcs
// 			dxbs[,i] = merlin_rcs(basis[,i],knots,1) * brcs
		}	
	}
	xbs = xbs :+ asarray(gml.distancb,(mod,2))
	return(exp(xbs))
}

`RM' merlin_aft_cdf(`gml' gml,`RC' t)
{
	return(1:-merlin_aft_s(gml,t))
}

`RM' merlin_aft_s(`gml' gml, `RC' t)
{
	return(exp(-merlin_aft_ch(gml,t)))
}

end


