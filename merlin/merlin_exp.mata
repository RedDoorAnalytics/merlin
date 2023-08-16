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

`RM' merlin_logl_exp(	`gml' gml, 	///
						`RM' G, 	///
						`RM' H)
{	
	model 	= gml.model
	y       = merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml,model)
	hast	= gml.istimedep[model,1]
	haslt	= gml.hasltrunc[model]
	hasbh	= gml.hasbh[model,]
	logl 	= J(Nobs,1,0)

	//core
	xb	= merlin_util_xzb(gml)

	//===========================================================================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 	= merlin_get_surv_index(gml)
			if 	(hasbh[1]) 	logl[index1,] = log(exp(xb[index1,]) :+ merlin_util_bhazard(gml))
			else if (hasbh[2]) 	logl[index1,] = log(exp(xb[index1,]) :* merlin_util_bhazard(gml))
			else			logl[index1,] = xb[index1,]
		}

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		Nobs2 		= merlin_get_nobs(gml)
		if (Nobs2) {
			index2 	= merlin_get_surv_index(gml)
			if (hast) {
				Ngq 	= gml.chip
				chq2 	= J(Nobs2,Ngq,0)
				gq 	= merlin_gq(Ngq,"legendre")
				if (!haslt) qp2	= y[index2,1] :/ 2 :* J(Nobs2,1,gq[,1]') :+ y[index2,1]:/2
				else 	    qp2	= (y[index2,1]:-y[index2,3]) :/ 2 :* J(Nobs2,1,gq[,1]') :+ (y[index2,1]:+y[index2,3]):/2
				for (q=1;q<=Ngq;q++) {
					chq2[,q] = exp(merlin_util_xzb(gml,qp2[,q]))
				}
				if (!haslt) logl[index2] = logl[index2] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
				else 		logl[index2] = logl[index2] :- (y[index2,1]:-y[index2,3]):/2 :* (chq2 * gq[,2]) 
			}
			else {
				logl[index2] = logl[index2] :- exp(xb[index2,]) :* y[index2,1]
				if (haslt) {
					gml.survind = 4
					index4 		= merlin_get_surv_index(gml)
					logl[index4] = logl[index4] :+ exp(xb[index4,]) :* y[index4,3]
				}
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

		//left truncation handled externally

	if (gml.todo==0) return(logl)

	//===========================================================================================================================//
	// score
	
		gml.survind = 0
		
		//indexes for covariates and baseline spline equations
		NHbs 	= asarray(gml.NHbs,model)
		
		//core
		sindex	= (1..NHbs[1]) :+ gml.skip[model]
		x  	= merlin_util_xz(gml)

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			G[index1,sindex] = x[index1,] 
		}  

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (merlin_get_nobs(gml,model)) {
			if (hast) {
				dch2 	= J(Nobs2,1,0)
				if (!haslt) qw2 = y[index2,1]:/2 :* J(Nobs2,1,gq[,2]')
				else		qw2 = (y[index2,1]:-y[index2,3]):/2 :* J(Nobs2,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dch2 = dch2 :+ chq2[,q] :* merlin_util_xz(gml,qp2[,q]) :* qw2[,q]
				}
				G[index2,sindex] = G[index2,sindex] :- dch2
			}
			else {
				G[index2,sindex] = G[index2,sindex] :- exp(xb[index2,]) :* y[index2,1] :* x[index2,]
				if (haslt) G[index4,sindex] = G[index4,sindex] :+ exp(xb[index4,]) :* y[index4,3] :* x[index4,]
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		if (Nobs3) {	
			//exit times
			if (hast) {
				dchq3 	= J(Nobs3,1,0)
				qw3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dchq3 = dchq3 :+ chq3[,q] :* qw3[,q] :* merlin_util_xz(gml,qp3[,q])
				}
				G[index3,sindex] = exp(-ch3) :* dchq3
			}
			else G[index3,sindex] = exp(-exp(xb[index3,]) :* y[index3,1]) :* exp(xb[index3,]) :* y[index3,1] :* merlin_util_xz(gml)
			
			//entry times
			gml.survind = 5
			if (hast) {
				dchq5 	= J(Nobs5,1,0)
				qw5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dchq5 = dchq5 :+ chq5[,q] :* qw5[,q] :* merlin_util_xz(gml,qp5[,q])
				}
				G[index5,sindex] = G[index5,sindex] :- exp(-ch5) :* dchq5
			}
			else G[index5,sindex] = G[index5,sindex] :- exp(-exp(xb[index5,]) :* y[index5,tind]) :* exp(xb[index5,]) :* y[index5,tind] :* merlin_util_xz(gml)
			G[index3,sindex] = G[index3,sindex] :/ exp(logl[index3,])
		}

	if (gml.todo==1) return(logl)
	
	//===========================================================================================================================//
	// hessian

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
		// 0
		
		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (Nobs2) {
			//xb
			if (hast) {
				dch 	= J(Nobs2,1,0)
				for (q=1;q<=Ngq;q++) {
					x = merlin_util_xz(gml,qp2[,q])
					dch = dch :+ chq2[,q] :* qw2[,q] :* x[,xindex1] :* x[,xindex2]
				}
				Hx[index2,] 	= Hx[index2,] :- dch
			}
			else {
				Hx[index2,] = Hx[index2,] :- exp(xb[index2,]) :* y[index2,1] :* x[index2,xindex1] :* x[index2,xindex2]
				if (haslt) Hx[index4,] = Hx[index4,] :+ exp(xb[index4,]) :* y[index4,3] :* x[index4,xindex1] :* x[index4,xindex2]
			}
		}
		
		//interval censoring -> cdf function
		gml.survind = 3
		if (merlin_get_nobs(gml,model)) {
			//xb - exit times
			dchq32 	= J(Nobs3,1,0)
			for (q=1;q<=Ngq;q++) {
				x3 = merlin_util_xz(gml,qp3[,q])
				dchq32 = dchq32 :+ chq3[,q] :* qw3[,q] :* x3[,xindex1] :* x3[,xindex2]
			}
			
			G13 = exp(-ch3) :* dchq3
			Hx[index3,]	= (exp(-ch3) :* dchq32 :- G13[,xindex1] :* dchq3[,xindex2]) :/ exp(logl[index3]) :- G13[,xindex1] :* G13[,xindex2] :/ exp(logl[index3]) :/ exp(logl[index3])
			
			//xb - entry times
			gml.survind = 5
			dchq52 	= J(Nobs5,1,0)
			for (q=1;q<=Ngq;q++) {
				x5 = merlin_util_xz(gml,qp5[,q])
				dchq52 = dchq52 :+ chq5[,q] :* qw5[,q] :* x5[,xindex1] :* x5[,xindex2]
			}
			G15 = exp(-ch5) :* dchq5
			Hx[index5,]	= Hx[index5,] :- (exp(-ch5) :* dchq52 :- G15[,xindex1] :* dchq5[,xindex2]) :/ exp(logl[index5]) :- G15[,xindex1] :* G15[,xindex2] :/ exp(logl[index5]) :/ exp(logl[index5])
							
		}

		//build H
		Hxsum 	= quadcolsum(Hx,1)
		Hm		= J(NHbs[1],NHbs[1],.)	
		el 		= 1
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

`RM' merlin_exp_logh(`gml' gml,`RC' t)
{
	logh 	= merlin_util_xzb(gml,t)
	if (gml.hasbh[gml.model,1]) {
		logh = log(exp(logh) :+ asarray(gml.bhazards,gml.model)[asarray(gml.surv_index,(gml.model,1))])
	}
	return(logh)
}

`RM' merlin_exp_ch(`gml' gml, `RC' t, | `RC' t0)
{
	hast0 = args()==3
	if (gml.NI[gml.model]) {
		nobs 	= rows(t)
		ch 		= J(nobs,1,0)
		Ngq 	= gml.chip
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
		loghazq = log(t:/2 :* J(nobs,1,gq[,2]'))
		
		if (hast0) {
			for (q=1;q<=Ngq;q++) {
				ch = ch :+ exp(merlin_util_xzb(gml,qp[,q],t0) :+ loghazq[,q])
			}
		}
		else {
			for (q=1;q<=Ngq;q++) {
				ch = ch :+ exp(merlin_util_xzb(gml,qp[,q]) :+ loghazq[,q])
			}
		}
		
		return(ch)
		
	}
	else return(exp(merlin_util_xzb(gml,t)) :* t)
}

`RM' merlin_exp_cdf(`gml' gml,`RC' t)
{
	return(1:-exp(-merlin_exp_ch(gml,t)))
}

`RM' merlin_exp_s(`gml' gml, | `RC' t)
{
	return(exp(-merlin_exp_ch(gml,t))) 
}


`RM' merlin_exp_expval(`gml' gml, | `RC' t)
{
	if (args()==1) 	return(1:/exp(merlin_util_xzb(gml)))
	else	 		return(1:/exp(merlin_util_xzb(gml,t)))
}


end
