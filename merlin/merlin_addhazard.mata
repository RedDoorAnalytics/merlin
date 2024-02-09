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
-> functions for family(addhazard, ...)
*/


`RM' merlin_logl_addhazard(	`gml' gml, 	///
                                `RM' G, 	///
                                `RM' H)
{	
	model 		= gml.model
	y 		= merlin_util_depvar(gml)
	Nobs		= merlin_get_nobs(gml,model)
	hast		= gml.istimedep[model,1]
	haslt		= gml.hasltrunc[model]
	hasbh		= gml.hasbh[model,]
	logl 		= J(Nobs,1,0)

	//====================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		Nobs1		= merlin_get_nobs(gml,model)
		if (Nobs1) {
			index1 			= merlin_get_surv_index(gml)
			xb1				= merlin_util_xzb(gml)
			if 		(hasbh[1]) 	xb1 = xb1 :+ merlin_util_bhazard(gml)
			else if (hasbh[2]) 	xb1 = xb1 :* merlin_util_bhazard(gml)
			logl[index1,] = log(xb1)
		}

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		Nobs2 		= merlin_get_nobs(gml)
		if (Nobs2) {
			index2 	= merlin_get_surv_index(gml)
			Ngq 	= gml.chip
			ch 		= J(Nobs2,Ngq,.)
			gq 		= merlin_gq(Ngq,"legendre")
			if (!haslt) qp = y[index2,1] :/ 2 :* J(Nobs2,1,gq[,1]') :+ y[index2,1]:/2
			else		qp = (y[index2,1] :- y[index2,3]) :/ 2 :* J(Nobs2,1,gq[,1]') :+ (y[index2,1]:+y[index2,3]):/2
			for (q=1;q<=Ngq;q++) {
				ch[,q] = merlin_util_xzb(gml,qp[,q])
			}
			if (!haslt) logl[index2] = logl[index2] :- y[index2,1] :/ 2 :* (ch * gq[,2])
			else 		logl[index2] = logl[index2] :- (y[index2,1]:-y[index2,3]) :/ 2 :* (ch * gq[,2])
		}

		//interval censoring -> cdf function
		gml.survind = 3
		Nobs3 		= merlin_get_nobs(gml)
		if (Nobs3) {
			//exit times
			index3 	= merlin_get_surv_index(gml)
			Ngq 	= gml.chip
			chq3 	= J(Nobs3,Ngq,.)
			gq 		= merlin_gq(Ngq,"legendre")
			qp3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,1]') :+ y[index3,1]:/2
			for (q=1;q<=Ngq;q++) {
				chq3[,q] = merlin_util_xzb(gml,qp3[,q])
			}
			ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
			logl[index3,] = 1:-exp(-ch3)
			
			//entry times
			gml.survind = 5
			if (gml.hasltrunc[model]) 	tind = 4
			else 						tind = 3
			index5 	= merlin_get_surv_index(gml)
			Nobs5 	= merlin_get_nobs(gml)
			chq5 	= J(Nobs5,Ngq,.)
			qp5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,1]') :+ y[index5,tind]:/2
			for (q=1;q<=Ngq;q++) {
				chq5[,q] = merlin_util_xzb(gml,qp5[,q])
			}
			ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
			logl[index5,] = logl[index5,] :- (1:-exp(-ch5))
			
			//logL
			logl[index3,] = log(logl[index3,])
		}

		//left truncation handled externally

	if (gml.todo==0) return(logl)

	//====================================================================//
	// score

		gml.survind = 0
		
		//indexes for covariates and baseline spline equations
		NHbs 		= asarray(gml.NHbs,model)
		sindex 		= (1..NHbs[1]) :+ gml.skip[model]
		
		//core
		x			= merlin_util_xz(gml)

		//exactly observed events -> log hazard function
		if (Nobs1) G[index1,sindex] = x[index1,] :/ xb1

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (Nobs2) {
			dch 	= J(Nobs2,1,0)
			if (!haslt) qw = y[index2,1] :/ 2 :* J(Nobs2,1,gq[,2]')
			else		qw = (y[index2,1] :- y[index2,3]) :/ 2 :* J(Nobs2,1,gq[,2]')
			for (q=1;q<=Ngq;q++) {
				dch = dch :+ qw[,q] :* merlin_util_xz(gml,qp[,q])
			}
			G[index2,sindex] = G[index2,sindex] :- dch
		}

		//interval censoring -> cdf function
		gml.survind = 3
		if (Nobs3) {	
			//exit times
			dchq3 	= J(Nobs3,1,0)
			qw3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,2]')
			for (q=1;q<=Ngq;q++) {
				dchq3 = dchq3 :+ qw3[,q] :* merlin_util_xz(gml,qp3[,q])
			}
			G[index3,sindex] = exp(-ch3) :* dchq3

			//entry times
			gml.survind = 5
			dchq5 	= J(Nobs5,1,0)
			qw5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,2]')
			for (q=1;q<=Ngq;q++) {
				dchq5 = dchq5 :+ qw5[,q] :* merlin_util_xz(gml,qp5[,q])
			}
			G[index5,sindex] = G[index5,sindex] :- exp(-ch5) :* dchq5
			G[index3,sindex] = G[index3,sindex] :/ exp(logl[index3,])
		}

	if (gml.todo==1) return(logl)
	
	//====================================================================//
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
		if (Nobs1) Hx[index1,] = -x[index1,xindex1] :* x[index1,xindex2] :/ (xb1:^2)
		
		//exactly observed events and/or right censoring -> survival function
		//0
		
		//interval censoring -> cdf function
		gml.survind = 3
		if (merlin_get_nobs(gml,model)) {
			//xb - exit times
			dchq32 	= J(Nobs3,1,0)
			for (q=1;q<=Ngq;q++) {
				x3 = merlin_util_xz(gml,qp3[,q])
				dchq32 = dchq32 :+ qw3[,q] :* x3[,xindex1] :* x3[,xindex2]
			}
			
			G13 = exp(-ch3) :* dchq3
			Hx[index3,]	= (exp(-ch3) :* dchq32 :- G13[,xindex1] :* dchq3[,xindex2]) :/ exp(logl[index3]) :- G13[,xindex1] :* G13[,xindex2] :/ exp(logl[index3]) :/ exp(logl[index3])
			
			//xb - entry times
			gml.survind = 5
			dchq52 	= J(Nobs5,1,0)
			for (q=1;q<=Ngq;q++) {
				x5 = merlin_util_xz(gml,qp5[,q])
				dchq52 = dchq52 :+ qw5[,q] :* x5[,xindex1] :* x5[,xindex2]
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

`RM' merlin_addhazard_logh(`gml' gml, `RC' t)
{
	logh = log(merlin_util_xzb(gml,t))
	if (gml.hasbh[gml.model,1]) {
		logh = log(exp(logh) :+ merlin_util_bhazard(gml))
	}
	return(logh)
}

`RM' merlin_addhazard_ch(`gml' gml, `RC' t, | `RC' t0)
{
	nobs 	= rows(t)
	ch 	= J(nobs,1,0)
	Ngq 	= gml.chip

	gq 	= merlin_gq(Ngq,"legendre")
	qp	= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
	qw	= t :/ 2 :* J(nobs,1,gq[,2]')
	if (args()==2) {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ merlin_util_xzb(gml,qp[,q]) :* qw[,q]
		}
	}
	else {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ merlin_util_xzb(gml,qp[,q],t0) :* qw[,q]
		}
	}
	return(ch)
}

`RM' merlin_addhazard_s(`gml' gml, `RC' t)
{
	return(exp(-merlin_addhazard_ch(gml,t)))
}

`RM' merlin_addhazard_cdf(`gml' gml, `RC' t)
{
	return(1:-merlin_addhazard_s(gml,t))
}

end
