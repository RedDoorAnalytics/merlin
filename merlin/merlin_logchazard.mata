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

/*
-> functions for family(logchazard, ...)
*/

`RM' merlin_logl_logchazard(`gml' gml, 	///
							`RM' G, 	///
							`RM' H)
{	
	model 		= gml.model
	y 			= merlin_util_depvar(gml)
	Nobs		= merlin_get_nobs(gml,model)
	hast		= gml.istimedep[model,1]
	haslt		= gml.hasltrunc[model]
	hasbh		= gml.hasbh[model,1]
	logl 		= J(Nobs,1,0)

	//core
	xb			= merlin_util_xzb(gml)
	expxb		= exp(xb)

	//===========================================================================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 			= merlin_get_surv_index(gml)
			dxb				= merlin_util_xzb_deriv(gml)
			logl[index1,] 	= xb[index1] :+ log(dxb)
			if (hasbh) logl[index1,] = log(exp(logl[index1,]) :+ merlin_util_bhazard(gml))
		}

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (merlin_get_nobs(gml,model)) {
			index2 			= merlin_get_surv_index(gml)
			logl[index2] 	= logl[index2] :- expxb[index2]
			if (haslt) {
				gml.survind	= 4
				index4 		= merlin_get_surv_index(gml)
				expxb0 		 = exp(merlin_util_xzb(gml,y[index4,3]))
				logl[index4] = logl[index4] :+ expxb0
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		if (merlin_get_nobs(gml,model)) {
			//exit times
			index3 = merlin_get_surv_index(gml)
			logl[index3,] = 1:-exp(-expxb[index3])
			//entry times
			gml.survind = 5
			if (gml.hasltrunc[model]) 	tind = 4
			else 						tind = 3
			index5 	= merlin_get_surv_index(gml)
			xbl0	= merlin_util_xzb(gml,y[index5,tind])
			logl[index5,] = logl[index5,] :- (1:-exp(-exp(xbl0)))
			//logL
			logl[index3,] = log(logl[index3,])
		}

		//left truncation handled externally

	if (gml.todo==0) return(logl)

	//===========================================================================================================================//
	// score
	
		gml.survind = 0
		
		//indexes for covariates and baseline spline equations
		NHbs 		= asarray(gml.NHbs,model)
		sindex 		= (1..NHbs[1]) :+ gml.skip[model]
		
		//core
		x  			= merlin_util_xz(gml)

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			if (hast) {
					dx = merlin_util_xz_deriv(gml)
					G[index1,sindex] = x[index1,] :+ dx :/ dxb
			}
			else 	G[index1,sindex] = x[index1,] 
		}  

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (merlin_get_nobs(gml,model)) {
			G[index2,sindex] = G[index2,sindex] :- x[index2,] :* expxb[index2]
			if (haslt) {
				x0 = merlin_util_xz(gml,y[index4,3])
				G[index4,sindex] = G[index4,sindex] :+ x0 :* expxb0
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		if (merlin_get_nobs(gml,model)) {
			//exit times
			expexpxb3 = exp(-expxb[index3]) :* expxb[index3]
			G[index3,sindex] = expexpxb3 :* x[index3,]
			//entry times
			gml.survind = 5
			xl0 = merlin_util_xz(gml,y[index5,tind])
			expexpxbl05 = exp(-exp(xbl0)) :* exp(xbl0)
			G[index5,sindex] = G[index5,sindex] :- expexpxbl05 :* xl0
			explogl3 = exp(logl[index3,])
			G[index3,sindex] = G[index3,sindex] :/ explogl3
		}

	if (gml.todo==1) return(logl)
	
	//===========================================================================================================================//
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
		
		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			if (hast) {
				//xb
				Hx[index1,] = -dx[,xindex1] :* dx[,xindex2] :/ dxb2
			}
			//else 0
		}  

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (merlin_get_nobs(gml,model)) {
			//xb
			Hx[index2,] 	= Hx[index2,] :- x[index2,xindex1] :* x[index2,xindex2] :* expxb[index2]
			if (haslt) Hx[index4,] 	= Hx[index4,] :+ x0[,xindex1] :* x0[,xindex2] :* expxb0
		}
		
		//interval censoring -> cdf function
		gml.survind = 3
		if (merlin_get_nobs(gml,model)) {
			explogl5 = exp(logl[index5])
			//xb - exit times
			Hx[index3,]		=  -expexpxb3 :* x[index3,xindex1] :* G[index3,xindex2] :/ explogl3 :+ x[index3,xindex1] :* expexpxb3 :* (x[index3,xindex2] :- expxb[index3] :* x[index3,xindex2]) :/ explogl3
			//xb - entry times
			Hx[index5,]		=  Hx[index5,] :- (-expexpxbl05 :* xl0[,xindex1] :* G[index5,xindex2] :/ explogl5 :+ 1:/ explogl5 :* xl0[,xindex1] :* expexpxbl05 :* (xl0[,xindex2] :- exp(xbl0) :* xl0[,xindex2]))
		}

		//build H
		Hxsum = quadcolsum(Hx,1)
		
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
		
		Hindex1 = 1 + gml.skip[model] 
		Hindex2 = NHbs[1] + gml.skip[model] 
		H[|Hindex1,Hindex1\Hindex2,Hindex2|] = Hm
		
		return(logl)
}

`RM' merlin_logchazard_h(`gml' gml, `RC' t)
{	
	return(exp(merlin_logchazard_logh(gml,t)))
}

`RM' merlin_logchazard_logh(`gml' gml, `RC' t)
{	
	logh = merlin_util_xzb(gml,t) :+ log(merlin_util_xzb_deriv(gml,t))
	if (gml.hasbh[gml.model,1]) {
		logh = log(exp(logh) :+ merlin_util_bhazard(gml))
	}
	return(logh)
}

`RM' merlin_logchazard_ch(`gml' gml, `RC' t, | `RC' t0)
{	
	if (args()==2)	return(exp(merlin_util_xzb(gml,t)))
	else 			return(exp(merlin_util_xzb(gml,t,t0)))
}

`RM' merlin_logchazard_s(`gml' gml, `RC' t)
{
	return(exp(-merlin_logchazard_ch(gml,t)))
}

`RM' merlin_logchazard_cdf(`gml' gml, `RC' t)
{	
	return(1 :- merlin_logchazard_s(gml,t))
}

end
