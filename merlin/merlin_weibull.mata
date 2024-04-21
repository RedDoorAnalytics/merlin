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

`RM' merlin_logl_weibull(`gml' gml, 	///
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
	xb	= merlin_util_xzb(gml)
	gam	= merlin_util_dap(gml,1)

	//===================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 			= merlin_get_surv_index(gml)
			logl[index1] 	= xb[index1] :+ log(gam) :+ (gam - 1) :* log(y[index1,1])
			if (hasbh[1]) logl[index1,] = log(exp(logl[index1,]) :+ merlin_util_bhazard(gml))
			if (hasbh[2]) logl[index1,] = log(exp(logl[index1,]) :* merlin_util_bhazard(gml))
		}

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		Nobs2 		= merlin_get_nobs(gml)
		if (Nobs2) {
			index2 	= merlin_get_surv_index(gml)
			if (hast | hasbh[2]) {
				Ngq 	= gml.chip
				chq2 	= J(Nobs2,Ngq,0)
				gq 		= merlin_gq(Ngq,"legendre")
				if (!haslt) qp2 = y[index2,1] :/ 2 :* J(Nobs2,1,gq[,1]') :+ y[index2,1]:/2
				else 		qp2 = (y[index2,1] :- y[index2,3]) :/ 2 :* J(Nobs2,1,gq[,1]') :+ (y[index2,1] :+ y[index2,3]):/2
				for (q=1;q<=Ngq;q++) {
					chq2[,q] = merlin_util_xzb(gml,qp2[,q]) 
				}
				if (hasbh[2]) chq2 = chq2 :+ log(merlin_util_bhazard(gml))
				chq2 = exp(chq2 :+ log(gam) :+ (gam-1) :* log(qp2))
				if (!haslt)	logl[index2] = logl[index2] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
				else 		logl[index2] = logl[index2] :- (y[index2,1]:-y[index2,3]):/2 :* (chq2 * gq[,2]) 
			}
			else {
				logl[index2] = logl[index2] :- exp(xb[index2,]) :* y[index2,1] :^ gam
				if (haslt) {
					gml.survind 	= 4
					index4 		= merlin_get_surv_index(gml)
					logl[index4] 	= logl[index4] :+ exp(xb[index4,]) :* y[index4,3] :^ gam
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
					chq3[,q] = exp(merlin_util_xzb(gml,qp3[,q])) :* gam :* qp3[,q] :^(gam-1)
				}
				ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
				logl[index3,] = -exp(-ch3)
			}
			else logl[index3,] = -exp(-exp(xb[index3,]) :* y[index3,1] :^gam)
			
			//entry times
			gml.survind = 5
			if (gml.hasltrunc[model]) 	tind = 4
			else 						tind = 3
			index5 = merlin_get_surv_index(gml)
			if (hast) {
				Nobs5 	= merlin_get_nobs(gml)
				chq5 	= J(Nobs5,Ngq,.)
				qp5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,1]') :+ y[index5,tind]:/2
				for (q=1;q<=Ngq;q++) {
					chq5[,q] = exp(merlin_util_xzb(gml,qp5[,q])) :* gam :* qp5[,q] :^(gam-1)
				}
				ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
				logl[index5,] = logl[index5,] :+ exp(-ch5)
			}
			else logl[index5,] = logl[index5,] :+ exp(-exp(xb[index5,]) :* y[index5,tind] :^gam)
			//logL
			logl[index3,] = log(logl[index3,])
		}

	if (gml.todo==0) return(logl)

	//===================================================================//
	// score
	
		gml.survind = 0
		
		//indexes for covariates and baseline equations
		NHbs 	= asarray(gml.NHbs,model)
		sindex1 = (1..NHbs[1]) :+ gml.skip[model]
		sindex2	= ((NHbs[1]+1)..(NHbs[1]+NHbs[2])) :+ gml.skip[model]
		
		//core
		x  	= merlin_util_xz(gml)

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			G[index1,sindex1] = x[index1,] 
			G[index1,sindex2] = 1 :+ gam :* log(y[index1,1])
		}  

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (merlin_get_nobs(gml,model)) {
			if (hast) {
				dch2 = dch22 = J(Nobs2,1,0)
				if (!haslt) qw2  = y[index2,1]:/2 :* J(Nobs2,1,gq[,2]')
				else 		qw2  = (y[index2,1]:-y[index2,3]):/2 :* J(Nobs2,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dchq2temp 	= chq2[,q] :* qw2[,q]
					dch2 		= dch2 :+ dchq2temp :* merlin_util_xz(gml,qp2[,q])
					dch22 		= dch22 :+ dchq2temp :* (gam:*log(qp2[,q]) :+ 1)
				}
				G[index2,sindex1] = G[index2,sindex1] :- dch2
				G[index2,sindex2] = G[index2,sindex2] :- dch22
			}
			else {
				Gtemp = exp(xb[index2,]) :* y[index2,1] :^ gam
				G[index2,sindex1] = G[index2,sindex1] :- Gtemp :* x[index2,]
				G[index2,sindex2] = G[index2,sindex2] :- Gtemp :* gam :* log(y[index2,1])
				if (haslt) {
					Gtemp = exp(xb[index4,]) :* y[index4,3] :^ gam
					G[index4,sindex1] = G[index4,sindex1] :+ Gtemp :* x[index4,]
					G[index4,sindex2] = G[index4,sindex2] :+ Gtemp :* gam :* log(y[index4,3])
				}
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		if (Nobs3) {	

			//exit times
			if (hast) {
				dchq3 	= dchq32 	= J(Nobs3,1,0)
				qw3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dchq3temp 	= chq3[,q] :* qw3[,q]
					dchq3 		= dchq3 :+ dchq3temp :* merlin_util_xz(gml,qp3[,q])
					dchq32		= dchq32 :+ dchq3temp :* (gam:*log(qp3[,q]) :+ 1)
				}
				G[index3,sindex1] = exp(-ch3) :* dchq3
				G[index3,sindex2] = exp(-ch3) :* dchq32
			}
			else {
				Gtemp = exp(-exp(xb[index3,]) :* y[index3,1] :^ gam) :* exp(xb[index3,]) :* y[index3,1] :^ gam
				G[index3,sindex1] = Gtemp :* merlin_util_xz(gml)
				G[index3,sindex2] = Gtemp :* gam :* log(y[index3,1])
			}

			//entry times
			gml.survind = 5
			if (hast) {
				dchq5 	= dchq52 	= J(Nobs5,1,0)
				qw5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dchq5temp 	= chq5[,q] :* qw5[,q]
					dchq5 		= dchq5 :+ dchq5temp :* merlin_util_xz(gml,qp5[,q])
					dchq52 		= dchq52 :+ dchq5temp :* (gam:*log(qp5[,q]) :+ 1)
				}
				G[index5,sindex1] = G[index5,sindex1] :- exp(-ch5) :* dchq5
				G[index5,sindex2] = G[index5,sindex2] :- exp(-ch5) :* dchq52
			}
			else {
				Gtemp = exp(-exp(xb[index5,]) :* y[index5,tind] :^ gam) :* exp(xb[index5,]) :* y[index5,tind] :^ gam
				G[index5,sindex1] = G[index5,sindex1] :- Gtemp :* merlin_util_xz(gml)
				G[index5,sindex2] = G[index5,sindex2] :- Gtemp :* gam :* log(y[index5,tind])
			}
			sindex 				= sindex1,sindex2
			G[index3,sindex] 	= G[index3,sindex] :/ exp(logl[index3,])
		}

	if (gml.todo==1) return(logl)
	
	//===================================================================//
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
		
		//x,gamma indices
		xindex3 = xindex4 = J(1,0,.)
		for (i=1; i<=NHbs[2]; i++) {
			for (j=1;j<=NHbs[1];j++) {
				xindex3 = xindex3,i
				xindex4 = xindex4,j
			}
		}
		Hxg = J(Nobs,cols(xindex3),0)
		
		//gamma indices
		xindex5 = xindex6 = J(1,0,.)
		for (i=1; i<=NHbs[2]; i++) {
			refind = 1
			while (refind<=i) {
				xindex5 = xindex5,i
				xindex6 = xindex6,refind
				refind++
			}
		}
		Hg = J(Nobs,1,0)

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {	
			// x - 0
			//xg - 0
			//g
			Hg[index1] = gam :* log(y[index1,1])
		}

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (Nobs2) {
			//xb
			if (hast) {
				dch = dchxg = J(Nobs2,1,0)
				logqp2	= log(qp2)
				if (!haslt) hazq = y[index2,1] :/ 2 :* J(Nobs2,1,gq[,2]') :* qp2 :^ (gam-1) :* gam :* (gam:*logqp2 :+ 1)
				else 		hazq = (y[index2,1] :- y[index2,3]) :/ 2 :* J(Nobs2,1,gq[,2]') :* qp2 :^ (gam-1) :* gam :* (gam:*logqp2 :+ 1)
				hazq2	= logqp2 :* gam :* qp2:^(gam:-1) :* ((gam:^2 :* logqp2 :+ 2:*gam) :+ gam  :+ 1:/logqp2)
				for (q=1;q<=Ngq;q++) {
					x = merlin_util_xz(gml,qp2[,q])
					dch = dch :+ chq2[,q] :* qw2[,q] :* x[,xindex1] :* x[,xindex2]
					expxb = exp(merlin_util_xzb(gml,qp2[,q]))
					dchxg = dchxg :+ expxb :* hazq[,q] :* merlin_util_xz(gml,qp2[,q])
					hazq2[,q] = hazq2[,q] :* expxb 
				}
				Hx[index2,] 	= -dch
				Hxg[index2,] 	= -dchxg
				if (!haslt) Hg[index2] = Hg[index2] :- y[index2,1]:/2 :* (hazq2 * gq[,2])
				else 		Hg[index2] = Hg[index2] :- (y[index2,1]:-y[index2,3]):/2 :* (hazq2 * gq[,2])
			}
			else {
				Htemp 		 = exp(xb[index2,]) :* y[index2,1] :^ gam
				Hx[index2,]  = -Htemp :* x[index2,xindex1] :* x[index2,xindex2]
				Hxg[index2,] = -Htemp :* gam :* log(y[index2,1]) :* x[index2,xindex4]
				Hg[index2] 	 = Hg[index2] :- log(y[index2,1]) :* gam :* Htemp :* (gam :* log(y[index2,1]) :+ 1)
				if (haslt) {
					Htemp 		 = exp(xb[index4,]) :* y[index4,3] :^ gam
					Hx[index4,]  = Hx[index4,] :+ Htemp :* x[index4,xindex1] :* x[index4,xindex2]
					Hxg[index4,] = Hxg[index4,] :+ Htemp :* gam :* log(y[index4,3]) :* x[index4,xindex4]
					Hg[index4] 	 = Hg[index4] :+ log(y[index4,3]) :* gam :* Htemp :* (gam :* log(y[index4,3]) :+ 1)
				}
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
		Hxsum 		= quadcolsum(Hx,1)
		Hxgsum	= quadcolsum(Hxg,1)
		Hgsum 	= quadcolsum(Hg,1)
		
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
				Hoff[e1,e2] = Hxgsum[el++]
			}					
		} 
		
		H2	= J(NHbs[2],NHbs[2],.)	
		el 	= 1
		for (e1=1;e1<=NHbs[2];e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) H2[e1,e1] = Hgsum[el++]
				else 		H2[e2,e1] = H2[e1,e2] = Hgsum[el++]
				e2++
			}
		}
		
		Hm = Hm,Hoff'\Hoff,H2
		
		Hindex1 = 1 + gml.skip[model] 
		Hindex2 = NHbs[1] + NHbs[2] :+ gml.skip[model] 
		H[|Hindex1,Hindex1\Hindex2,Hindex2|] = Hm
		
		return(logl)
}

`RM' merlin_logl_weibull_ml(`gml' gml)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml,model)
	hast	= gml.istimedep[model,1]
        haslt	= gml.hasltrunc[model]
	hasbh	= gml.hasbh[model,1]
	logl 	= J(Nobs,gml.ndim[gml.Nrelevels],0)

	//core
	xb	= merlin_util_xzb(gml)
	gam	= merlin_util_dap(gml,1)

	//====================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 		= merlin_get_surv_index(gml)
			logl[index1,] 	= xb[index1,] :+ log(gam) :+ (gam:-1) :* log(y[index1,1])
			if (hasbh) logl[index1,] = log(exp(logl[index1,]) :+ merlin_util_bhazard(gml))
		}

		//exactly observed events and/or right censoring 
                // -> survival function
		gml.survind = 2
		Nobs2 		= merlin_get_nobs(gml)
		if (Nobs2) {
			index2 	= merlin_get_surv_index(gml)
			if (hast) {
				Ngq 	= gml.chip
				chq2 	= J(Nobs2,1,0)
				gq 	= merlin_gq(Ngq,"legendre")
                                
                                if (!haslt) {
                                        qp2 = y[index2,1] :/ 2 :* 
                                                J(Nobs2,1,gq[,1]') :+ 
                                                y[index2,1] :/ 2
                                        loghazq2 = log(y[index2,1])
                                }
				else {
                                        qp2 = (y[index2,1] :- y[index2,3]) :/ 
                                                2 :* J(Nobs2,1,gq[,1]') :+ 
                                                (y[index2,1] :+ 
                                                        y[index2,3]) :/ 2
                                        loghazq2 = log(y[index2,1] :+ 
                                                y[index2,3])
                                }

				loghazq2 = loghazq2 :+ log(gam) :+ (gam-1) :* 
                                        log(qp2) :+ J(Nobs2,1,log(gq[,2]')) :- 
                                        log(2)

				for (q=1;q<=Ngq;q++) {
					chq2 = chq2 :+ 
                                          exp(merlin_util_xzb(gml,qp2[,q]) :+ 
                                          loghazq2[,q])
				}

				logl[index2,] = logl[index2,] :- chq2
			}
			else {
                                logl[index2,] = logl[index2,] :- 
                                        exp(xb[index2,]) :* y[index2,1] :^ gam
                                if (haslt) {
					gml.survind 	= 4
					index4 	= merlin_get_surv_index(gml)
					logl[index4,] 	= logl[index4,] :+ 
                                                exp(xb[index4,]) :* 
                                                y[index4,3] :^ gam
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
				gq 	= merlin_gq(Ngq,"legendre")
				qp3	= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,1]') :+ y[index3,1]:/2
				for (q=1;q<=Ngq;q++) {
					chq3[,q] = exp(merlin_util_xzb(gml,qp3[,q])) :* gam :* qp3[,q] :^(gam-1)
				}
				ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
				logl[index3,] = -exp(-ch3)
			}
			else logl[index3,] = -exp(-exp(xb[index3,]) :* y[index3,1] :^gam)
			
			//entry times
			gml.survind = 5
			if (gml.hasltrunc[model]) 	tind = 4
			else 				tind = 3
			index5 = merlin_get_surv_index(gml)
			if (hast) {
				Nobs5 	= merlin_get_nobs(gml)
				chq5 	= J(Nobs5,Ngq,.)
				qp5	= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,1]') :+ y[index5,tind]:/2
				for (q=1;q<=Ngq;q++) {
					chq5[,q] = exp(merlin_util_xzb(gml,qp5[,q])) :* gam :* qp5[,q] :^(gam-1)
				}
				ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
				logl[index5,] = logl[index5,] :+ exp(-ch5)
			}
			else logl[index5,] = logl[index5,] :+ exp(-exp(xb[index5,]) :* y[index5,tind] :^gam)
			//logL
			logl[index3,] = log(logl[index3,])
		}

		//left truncation handled externally
	
	return(logl)
}


`RM' merlin_weibull_logh(`gml' gml,`RC' t)
{
	gam 	= asarray(gml.distancb,(gml.model,1))
	logh 	= merlin_util_xzb(gml,t) :+ log(gam) :+ (gam - 1) :* log(t)
	if (gml.hasbh[gml.model,1]) {
		logh = log(exp(logh) :+ asarray(gml.bhazards,gml.model)[asarray(gml.surv_index,(gml.model,1))])
	}
	return(logh)
}

`RM' merlin_weibull_ch(`gml' gml,`RC' t, | `RC' t0)
{
	gam 	= asarray(gml.distancb,(gml.model,1))
	hast0 	= args()==3
	
	if (gml.NI[gml.model]) {
		nobs 	= rows(t)
		ch 	= J(nobs,1,0)
		Ngq 	= gml.chip
		gq 	= merlin_gq(Ngq,"legendre")
		qp	= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
		loghazq = log(gam) :+ (gam-1) :* log(qp) :+ log(t:/2 :* J(nobs,1,gq[,2]'))
		
		if (!hast0) {
			for (q=1;q<=Ngq;q++) {
				ch = ch :+ exp(merlin_util_xzb(gml,qp[,q]) :+ loghazq[,q])
			}
		}
		else {
			for (q=1;q<=Ngq;q++) {
				ch = ch :+ exp(merlin_util_xzb(gml,qp[,q],t0) :+ loghazq[,q])
			}
		}
		return(ch)
	}
	else {
		if (!hast0)	return(exp(merlin_util_xzb(gml,t)) :* t :^ gam)
		else 		return(exp(merlin_util_xzb(gml,t,t0)) :* t :^ gam)
	}
}

`RM' merlin_weibull_cdf(`gml' gml,`RC' t)
{
	return(1:-merlin_weibull_s(gml,t))
}

`RM' merlin_weibull_s(`gml' gml, `RC' t)
{
	return(exp(-merlin_weibull_ch(gml,t))) 	
}

`RM' merlin_weibull_expval(`gml' gml, `RC' t)
{
	
}

end
