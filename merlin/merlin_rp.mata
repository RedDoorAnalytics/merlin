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

/*
-> functions for family(rp, ...)
*/

mata:

`RM' merlin_logl_rp(`gml' gml, 	///
                `RM' G, 	///
                `RM' H)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml,model)
	hast	= gml.istimedep[model,1]
	haslt	= gml.hasltrunc[model]
	hasbh	= gml.hasbh[model,]
	logl 	= J(Nobs,1,0)

	//core
	rcs	= asarray(gml.distancb,(model,2))
	drcs	= asarray(gml.distancb,(model,6))
	brcs 	= asarray(gml.distancb,(model,1))

	xb	= merlin_util_xzb(gml) :+ rcs * brcs
	expxb	= exp(xb)

	//===================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		Nobs1	= merlin_get_nobs(gml,model)
		if (Nobs1) {
			index1 		= merlin_get_surv_index(gml)
			if (hast)	dxb = merlin_util_xzb_deriv(gml)  :+ (drcs[index1,] * brcs) :/ y[index1,1]
			else 		dxb = (drcs[index1,] * brcs) :/ y[index1,1]
			logl[index1,] 	= xb[index1] :+ log(dxb) 
			if (hasbh[1]) {
				totalh1		= exp(logl[index1,]) :+ merlin_util_bhazard(gml)
				logl[index1,] 	= log(totalh1)
			}
			if (hasbh[2]) {
				logl[index1,] = log(exp(xb[index1]) :* merlin_util_bhazard(gml) :+ exp(logl[index1]) :* merlin_util_bHazard(gml))
			}
		}

		//exactly observed events and/or right censoring -> CH function
		gml.survind = 2
		Nobs2		= merlin_get_nobs(gml,model)
		if (Nobs2) {
			index2 	= merlin_get_surv_index(gml)
			if (hasbh[2])	logl[index2] 	= logl[index2] :- expxb[index2] :* merlin_util_bHazard(gml)
			else		logl[index2] 	= logl[index2] :- expxb[index2]
			if (haslt) {
				gml.survind 	= 4
				index4 		= merlin_get_surv_index(gml)
				rcs0 		= asarray(gml.distancb,(model,8))
				if (hast) 	xb0 = merlin_util_xzb(gml,y[index4,3])
				else		xb0 = merlin_util_xzb(gml)
				expxb0 		= exp(xb0 :+ rcs0 * brcs)
				/*if (hasbh[2])	logl[index4] 	= logl[index4] :+ expxb0 :* merlin_util_bHazard(gml)
				else*/		logl[index4] 	= logl[index4] :+ expxb0
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		Nobs3 		= merlin_get_nobs(gml,model)
		if (Nobs3) {
			index3 	= merlin_get_surv_index(gml)
			if (haslt) 	tind = 4
			else 		tind = 3
			if (gml.hastmat) {
				result 		= J(Nobs3,1,0)
				//find where to go from starting state of gml.model transition
				start		= sum((1::rows(gml.tmat)) :* rowsum(model:==gml.tmat))
				Ns 		= sum(gml.tmat[start,]:!=.)
				transind 	= select(1::rows(gml.tmat),gml.tmat[start,]':!=.)'
				posstates 	= gml.tmat[start,transind]
				Ngq 		= gml.chip
				gq 		= merlin_gq(Ngq,"legendre")
				//left 0s handled fine with numerical integration
				scale 		= (y[index3,1]:-y[index3,tind]) :/ 2
				qp		= scale :* J(Nobs3,1,gq[,1]') :+ (y[index3,1]:+y[index3,tind]):/2
				cr 		= J(Nobs3,Ngq,0)
				for (q=1;q<=Ngq;q++) {
					for (i=1;i<=Ns;i++) {
						mod = posstates[i]
						if (gml.issurv[mod]) {
							if (mod==model) {
								cr[,q] = cr[,q] :+ (*gml.Plogh[mod])(gml,qp[,q]) :- (*gml.Pch[mod])(gml,qp[,q])
							}
							else {
								gml.model 	= mod
								cr[,q] 		= cr[,q] :- (*gml.Pch[mod])(gml,qp[,q])
								gml.model 	= model
							}	
						}
					}
				}
				logqw		 = log(scale :* J(Nobs3,1,gq[,2]'))
				logl[index3] = log(quadrowsum(exp(cr :+ logqw),1))
			}
			else {
				//exit times
				logl[index3,] 	= 1:-exp(-expxb[index3])
				//entry times
				gml.survind = 5
				index5 	= merlin_get_surv_index(gml)
				rcsl0	= asarray(gml.distancb,(model,7))
				xbl0	= merlin_util_xzb(gml,y[index5,tind]) :+ rcsl0 * brcs
				logl[index5,] = logl[index5,] :- (1:-exp(-exp(xbl0)))
				//logL
				logl[index3,] = log(logl[index3,])
			}
		}
	
	if (gml.haspenalty) logl = logl :- merlin_get_penalty(gml):/gml.N

	if (gml.todo==0) return(logl)

	//===================================================================//
	// score
	
		gml.survind = 0
	
		//indexes for covariates and baseline spline equations
		NHbs 	= asarray(gml.NHbs,model)
		sindex1 = (1..NHbs[1]) :+ gml.skip[model]
		sindex2	= ((NHbs[1]+1)..(NHbs[1]+NHbs[2])) :+ gml.skip[model]

		//core
		x  = merlin_util_xz(gml)

		//exactly observed events -> hazard function
		if (Nobs1) {		
			gml.survind = 1
			if (hasbh[1]) {
				if (hast) {
					dx1 		  = merlin_util_xz_deriv(gml)
					G[index1,sindex1] = expxb[index1,] :* (x[index1,] :* dxb :+ dx1) :/ totalh1
					G[index1,sindex2] = expxb[index1,] :* (drcs[index1,]:/ y[index1,1] :+ rcs[index1,] :* dxb) :/ totalh1
				}
				else {
					G[index1,sindex1] = x[index1,] :* expxb[index1,] :* (drcs[index1,] * brcs) :/ y[index1,1] :/ totalh1
					G[index1,sindex2] = expxb[index1,] :* (drcs[index1,] :+ rcs[index1,] :* (drcs[index1,] * brcs)) :/ totalh1 :/ y[index1,1]
				}
			}
			else {
				if (hast) {
                                        dx1 = merlin_util_xz_deriv(gml)
                                        G[index1,sindex1] = x[index1,] :+ dx1 :/ dxb
				}
				else 	G[index1,sindex1] = x[index1,] 
				G[index1,sindex2] = rcs[index1,] :+ drcs[index1,] :/ dxb :/ y[index1,1]
			}
		}  

		//exactly observed events and/or right censoring -> survival function
		if (Nobs2) {
			gml.survind = 2
			G[index2,sindex1] 		= G[index2,sindex1] :- x[index2,] :* expxb[index2]
			G[index2,sindex2] 		= G[index2,sindex2] :- rcs[index2,] :* expxb[index2]
			if (haslt) {
				gml.survind 		= 4
				x0 			= merlin_util_xz(gml,y[index4,3])
				G[index4,sindex1] 	= G[index4,sindex1] :+ x0 :* expxb0
				G[index4,sindex2] 	= G[index4,sindex2] :+ rcs0 :* expxb0
			} 
		}

		//interval censoring -> cdf function
		if (Nobs3) {
			//exit times
			expexpxb3 		= exp(-expxb[index3]) :* expxb[index3]
			G[index3,sindex1] 	= expexpxb3 :* x[index3,]
			G[index3,sindex2] 	= expexpxb3 :* rcs[index3,]
			//entry times
			gml.survind 		= 5
			xl0 			= merlin_util_xz(gml,y[index5,tind])
			expexpxbl05 		= exp(-exp(xbl0)) :* exp(xbl0)
			G[index5,sindex1] 	= G[index5,sindex1] :- expexpxbl05 :* xl0
			G[index5,sindex2] 	= G[index5,sindex2] :- expexpxbl05 :* rcsl0
			explogl3 			= exp(logl[index3,])
			sindex 				= sindex1,sindex2
			G[index3,sindex] 	= G[index3,sindex] :/ explogl3
		}

	if (gml.haspenalty) {
		cmpbix	= asarray(gml.CmpXBIndex,(gml.model,2))
		bindex 	= 1..(max(cmpbix)-1)
		G[,bindex] = G[,bindex] :- merlin_get_deriv_penalty(gml):/gml.N
	}
		
	if (gml.todo==1) return(logl)

	//===================================================================//
	// hessian

		//core
		dxb2 = dxb:^2

		//x indices
		xindex1 = xindex2 = J(1,0,.)
		for (i=1; i<=NHbs[1]; i++) {
			refind = 1
			while (refind<=i) {
				xindex1 = xindex1,i
				xindex2 = xindex2,refind
				refind++
			}
		}
		Hx = J(Nobs,cols(xindex1),0)
		
		//x,rcs indices
		xindex3 = xindex4 = J(1,0,.)
		for (i=1; i<=NHbs[2]; i++) {
			for (j=1;j<=NHbs[1];j++) {
				xindex3 = xindex3,i
				xindex4 = xindex4,j
			}
		}
		Hxrcs = J(Nobs,cols(xindex3),0)
		
		//rcs indices
		xindex5 = xindex6 = J(1,0,.)
		for (i=1; i<=NHbs[2]; i++) {
			refind = 1
			while (refind<=i) {
				xindex5 = xindex5,i
				xindex6 = xindex6,refind
				refind++
			}
		}
		Hrcs = J(Nobs,cols(xindex5),0)
	
		//exactly observed events -> hazard function
		if (Nobs1) {
                        if (hast) {
				//xb
				Hx[index1,] 	= -dx1[,xindex1] :* dx1[,xindex2] :/ dxb2
				//xbrcs
				Hxrcs[index1,] 	= -drcs[index1,xindex3] :* dx1[,xindex4] :/ y[index1,1] :/ dxb2
			}
			//rcs
			Hrcs[index1,] 		= -(drcs[index1,xindex5] :* drcs[index1,xindex6] :/ y[index1,1]:^2) :/ dxb2
		}  

		//exactly observed events and/or right censoring -> survival function
		if (Nobs2) {
			//xb
			Hx[index2,] 	= Hx[index2,] :- x[index2,xindex1] :* x[index2,xindex2] :* expxb[index2]
			//xbrcs
			Hxrcs[index2,] 	= Hxrcs[index2,] :- rcs[index2,xindex3] :* x[index2,xindex4]  :* expxb[index2]
			//rcs
			Hrcs[index2,] 	= Hrcs[index2,] :- rcs[index2,xindex5] :* rcs[index2,xindex6] :* expxb[index2]
			if (haslt) {
				gml.survind = 4
				//xb
				Hx[index4,] 	= Hx[index4,] :+ x0[,xindex1] :* x0[,xindex2] :* expxb0
				//xbrcs
				Hxrcs[index4,] 	= Hxrcs[index4,] :+ rcs0[,xindex3] :* x0[,xindex4]  :* expxb0
				//rcs
				Hrcs[index4,] 	= Hrcs[index4,] :+ rcs0[,xindex5] :* rcs0[,xindex6] :* expxb0
			}
		}
		
		//interval censoring -> cdf function
		if (Nobs3) {
			explogl5 = exp(logl[index5])
			//xb - exit times
			Hx[index3,]		=  -expexpxb3 :* x[index3,xindex1] :* G[index3,xindex2] :/ explogl3 :+ x[index3,xindex1] :* expexpxb3 :* (x[index3,xindex2] :- expxb[index3] :* x[index3,xindex2]) :/ explogl3
			//xb - entry times
			Hx[index5,]		=  Hx[index5,] :- (-expexpxbl05 :* xl0[,xindex1] :* G[index5,xindex2] :/ explogl5 :+ 1:/ explogl5 :* xl0[,xindex1] :* expexpxbl05 :* (xl0[,xindex2] :- exp(xbl0) :* xl0[,xindex2]))			
			//x,rcs - exit times
			Hxrcs[index3,]	=  -expexpxb3 :* rcs[index3,xindex3] :* G[index3,xindex4] :/ explogl3 :+ rcs[index3,xindex3] :* expexpxb3 :* (x[index3,xindex4] :- expxb[index3] :* x[index3,xindex4]) :/ explogl3
			//x,rcs - entry times
			Hxrcs[index5,]	=  Hxrcs[index5,] :- (-expexpxbl05 :* rcsl0[,xindex3] :* G[index5,xindex4] :/ explogl5 :+ 1 :/ explogl5 :* rcsl0[,xindex3] :* expexpxbl05 :* (xl0[,xindex4] :- exp(xbl0) :* xl0[,xindex4]))
			//rcs - exit times
			Hrcs[index3,]	=  -expexpxb3 :* rcs[index3,xindex5] :* G[index3,NHbs[1]:+xindex6] :/ explogl3 :+ rcs[index3,xindex5] :* expexpxb3 :* (rcs[index3,xindex6] :- expxb[index3] :* rcs[index3,xindex6]) :/ explogl3
			//rcs - entry times
			Hrcs[index5,]	=  Hrcs[index5,] :- (-expexpxbl05 :* rcsl0[,xindex5] :* G[index5,NHbs[1]:+xindex6] :/ explogl5 :+ 1:/ explogl5 :* rcsl0[,xindex5] :* expexpxbl05 :* (rcsl0[,xindex6] :- exp(xbl0) :* rcsl0[,xindex6]))
		}

		//build H
		Hxsum 		= quadcolsum(Hx,1)
		Hxrcssum	= quadcolsum(Hxrcs,1)
		Hrcssum 	= quadcolsum(Hrcs,1)
		
		Hm	= J(NHbs[1],NHbs[1],.)	
		el 	= 1
		for (e1=1;e1<=NHbs[1];e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) Hm[e1,e1] = Hxsum[el++]
				else 		Hm[e2,e1] = Hm[e1,e2] = Hxsum[el++]
				e2++
			}
		}
		
		Hoff = J(NHbs[2],NHbs[1],.)
		el = 1				
		for (e1=1;e1<=NHbs[2];e1++) {
			for (e2=1;e2<=NHbs[1];e2++) {
				Hoff[e1,e2] = Hxrcssum[el++]
			}					
		} 
		
		H2	= J(NHbs[2],NHbs[2],.)	
		el 	= 1
		for (e1=1;e1<=NHbs[2];e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) H2[e1,e1] = Hrcssum[el++]
				else 		H2[e2,e1] = H2[e1,e2] = Hrcssum[el++]
				e2++
			}
		}
		
		Hm = Hm,Hoff'\Hoff,H2
		
		Hindex1 = 1 + gml.skip[model] 
		Hindex2 = NHbs[1] + NHbs[2] + gml.skip[model] 
		H[|Hindex1,Hindex1\Hindex2,Hindex2|] = Hm
		
		return(logl)
}

`RM' merlin_rp_logch(`gml' gml, `RC' t, | `RC' t0)
{
	brcs 	= asarray(gml.distancb,(gml.model,1))
	knots 	= asarray(gml.distancb,(gml.model,3))
	hasorth = asarray(gml.distancb,(gml.model,4))
	if (hasorth) {
		rmat = asarray(gml.distancb,(gml.model,5))
		if (args()==2) 	return(merlin_util_xzb(gml,t) :+ merlin_rcs(log(t),knots,0,rmat) * brcs )
		else 			return(merlin_util_xzb(gml,t,t0) :+ merlin_rcs(log(t),knots,0,rmat) * brcs )
	}
	else {
		if (args()==2) 	return(merlin_util_xzb(gml,t) :+ merlin_rcs(log(t),knots,0) * brcs )
		else			return(merlin_util_xzb(gml,t,t0) :+ merlin_rcs(log(t),knots,0) * brcs )
	}
}

`RM' merlin_rp_ch(`gml' gml, `RC' t, | `RC' t0)
{
	if (args()==2) 	return(exp(merlin_rp_logch(gml,t)))
	else 			return(exp(merlin_rp_logch(gml,t,t0)))
}

`RM' merlin_rp_s(`gml' gml, `RC' t)
{
	return(exp(-merlin_rp_ch(gml,t)))
}

`RM' merlin_rp_logh(`gml' gml, `RC' t)
{
	brcs 	= asarray(gml.distancb,(gml.model,1))
	knots 	= asarray(gml.distancb,(gml.model,3))
	hasorth = asarray(gml.distancb,(gml.model,4))
	logt    = log(t)
	if (hasorth) {
		rmat = asarray(gml.distancb,(gml.model,5))
		logh = merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0,rmat) * brcs 
		logh = logh :+ log((merlin_rcs(logt,knots,1,rmat) * brcs) :/ t :+ merlin_util_xzb_deriv(gml,t))
	}
	else {
		logh = merlin_util_xzb(gml,t) :+ merlin_rcs(logt,knots,0) * brcs 
		logh = logh :+ log((merlin_rcs(logt,knots,1) * brcs) :/ t:+ merlin_util_xzb_deriv(gml,t) )
	}
	if (gml.hasbh[gml.model,1]) {
		logh = log(exp(logh) :+ merlin_util_bhazard(gml))
	}
	return(logh)
}

`RM' merlin_rp_h(`gml' gml, `RC' t)
{	
	return(exp(merlin_rp_logh(gml,t)))
}


`RM' merlin_rp_cdf(`gml' gml, `RC' t)
{	
	return(1 :- merlin_rp_s(gml,t))
}

end
