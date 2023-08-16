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

mata:

`RM' merlin_logl_gompertz(	`gml' gml, 	///
                                `RM' G, 	///
                                `RM' H)
{	
	model 		= gml.model
	y 			= merlin_util_depvar(gml)
	Nobs		= merlin_get_nobs(gml,model)
	hast		= gml.istimedep[model,1]
	haslt		= gml.hasltrunc[model]
	hasbh		= gml.hasbh[model,]
	logl 		= J(Nobs,1,0)

	//core
	xb			= merlin_util_xzb(gml)
	gam			= merlin_util_dap(gml,1)

	//===========================================================================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 			= merlin_get_surv_index(gml)
			logl[index1,] 	= xb[index1,] :+ gam :* y[index1,1]
			if (hasbh[1])	logl[index1,] = log(exp(logl[index1,]) :+ merlin_util_bhazard(gml))
			if (hasbh[2]) 	logl[index1,] = log(exp(logl[index1,]) :* merlin_util_bhazard(gml))
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
					chq2[,q] = exp(merlin_util_xzb(gml,qp2[,q]) :+ gam :* qp2[,q]) 
				}
				if (!haslt) logl[index2] = logl[index2] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
				else 		logl[index2] = logl[index2] :- (y[index2,1]:-y[index2,3]):/2 :* (chq2 * gq[,2]) 
			}
			else {
				logl[index2] = logl[index2] :- exp(xb[index2,]) :* (1/gam) :*(exp(gam:*y[index2,1]):-1)
				if (haslt) logl[index2] = logl[index2] :+ exp(xb[index2,]) :* (1/gam) :*(exp(gam:*y[index2,3]):-1)
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
					chq3[,q] = exp(merlin_util_xzb(gml,qp3[,q]) :+ gam :* qp3[,q])
				}
				ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
				logl[index3,] = 1:-exp(-ch3)
			}
			else logl[index3,] = -exp(-exp(xb[index3,]):/gam :* (exp(gam:*y[index3,1]):-1))
			
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
					chq5[,q] = exp(merlin_util_xzb(gml,qp5[,q]) :+ gam :* qp5[,q])
				}
				ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
				logl[index5,] = logl[index5,] :- (1:-exp(-ch5))
			}
			else logl[index5,] = logl[index5,] :+ exp(-exp(xb[index5,]):/gam :* (exp(gam:*y[index5,tind]):-1))
			//logL
			logl[index3,] = log(logl[index3,])
		}

	if (gml.todo==0) return(logl)

	//===========================================================================================================================//
	// score

		gml.survind = 0
		
		//indexes for covariates and baseline equations
		NHbs 		= asarray(gml.NHbs,model)
		sindex1 	= (1..NHbs[1]) :+ gml.skip[model]
		sindex2		= ((NHbs[1]+1)..(NHbs[1]+NHbs[2])) :+ gml.skip[model]
		
		//core
		x  			= merlin_util_xz(gml)

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			G[index1,sindex1] = x[index1,] 
			G[index1,sindex2] = y[index1,1]
		}  

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (merlin_get_nobs(gml,model)) {
			if (hast) {
				dch2 = dch22 = J(Nobs2,1,0)
				if (!haslt) qw2  = y[index2,1]:/2 :* J(Nobs2,1,gq[,2]')
				else 		qw2  = (y[index2,1]:-y[index2,3]):/2 :* J(Nobs2,1,gq[,2]')
				dchq2temp 	= chq2 :* qw2
				for (q=1;q<=Ngq;q++) {	
					dch2 		= dch2 :+ dchq2temp[,q] :* merlin_util_xz(gml,qp2[,q])
					dch22 		= dch22 :+ dchq2temp[,q] :* qp2[,q]
				}
				G[index2,sindex1] = G[index2,sindex1] :- dch2
				G[index2,sindex2] = G[index2,sindex2] :- dch22
			}
			else {
				Gtemp = exp(xb[index2,]) :/ gam 
				G[index2,sindex1] = G[index2,sindex1] :- Gtemp :* (exp(gam :* y[index2,1]) :-1) :* x[index2,]
				G[index2,sindex2] = G[index2,sindex2] :- Gtemp :* (exp(gam :* y[index2,1]) :* (y[index2,1] :- 1:/gam) :+ 1:/gam)
				if (haslt) {
					G[index2,sindex1] = G[index2,sindex1] :+ Gtemp :* (exp(gam :* y[index2,3]) :-1) :* x[index2,]
					G[index2,sindex2] = G[index2,sindex2] :+ Gtemp :* (exp(gam :* y[index2,3]) :* (y[index2,3] :- 1:/gam) :+ 1:/gam)
				}
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		if (Nobs3) {

			//exit times
			if (hast) {
				dchq31 	= dchq32 	= J(Nobs3,1,0)
				qw3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,2]')
				dchq3 	= chq3 :* qw3
				for (q=1;q<=Ngq;q++) {
					dchq31 		= dchq31 :+ dchq3[,q] :* merlin_util_xz(gml,qp3[,q])
					dchq32		= dchq32 :+ dchq3[,q] :* qp3[,q]
				}
				G[index3,sindex1] = exp(-ch3) :* dchq31
				G[index3,sindex2] = exp(-ch3) :* dchq32
			}
			else {
				Gtemp0 = exp(xb[index3,]) :/ gam
				Gtemp1 = exp(gam:* y[index3,1])
				Gtemp  = exp(-Gtemp0 :* (Gtemp1:-1))
				G[index3,sindex1] = Gtemp :* Gtemp0 :* (Gtemp1:-1) :* merlin_util_xz(gml)
				G[index3,sindex2] = Gtemp :* Gtemp0 :* (y[index3,1] :* Gtemp1 :- Gtemp1 :/ gam :+ 1:/gam)  
			}

			//entry times
			gml.survind = 5
			if (hast) {
				dchq51 	= dchq52 	= J(Nobs5,1,0)
				qw5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,2]')
				dchq5	= chq5 :* qw5
				for (q=1;q<=Ngq;q++) {
					dchq51 		= dchq51 :+ dchq5[,q] :* merlin_util_xz(gml,qp5[,q])
					dchq52 		= dchq52 :+ dchq5[,q] :* qp5[,q]
				}
				G[index5,sindex1] = G[index5,sindex1] :- exp(-ch5) :* dchq51
				G[index5,sindex2] = G[index5,sindex2] :- exp(-ch5) :* dchq52
			}
			else {
				Gtemp0 = exp(xb[index5,]) :/ gam
				Gtemp1 = exp(gam:* y[index5,tind])
				Gtemp  = exp(-Gtemp0 :* (Gtemp1:-1))
				G[index5,sindex1] = G[index5,sindex1] :- Gtemp :* Gtemp0 :* (Gtemp1:-1) :* merlin_util_xz(gml)
				G[index5,sindex2] = G[index5,sindex2] :- Gtemp :* Gtemp0 :* (y[index5,tind] :* Gtemp1 :- Gtemp1 :/ gam :+ 1:/gam)  
			}
			sindex 		= sindex1,sindex2
			G[index3,] 	= G[index3,sindex] :/ exp(logl[index3,])
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
		//all 0 

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (Nobs2) {
			//xb
			if (hast) {
				dchx = dchxg = dchg = J(Nobs2,1,0)
				for (q=1;q<=Ngq;q++) {	
					x 		= merlin_util_xz(gml,qp2[,q])
					dchx 	= dchx :+ dchq2temp[,q] :* x[,xindex1] :* x[,xindex2]
					dchxg 	= dchxg :+ dchq2temp[,q] :* qp2[,q] :* x[,xindex4]
					dchg 	= dchg :+ dchq2temp[,q] :* qp2[,q]:^2
				}
				Hx[index2,] 	= -dchx
				Hxg[index2,] 	= -dchxg
				Hg[index2,] 	= -dchg
			}
			else {
				Htemp 		 = exp(xb[index2,]) :/ gam 
				Hx[index2,]  = -Htemp :* x[index2,xindex1] :* x[index2,xindex2] :* (exp(gam:*y[index2,1]) :-1)
				Hxg[index2,] = -Htemp :* x[index2,xindex4] :* (y[index2,1] :* exp(gam:*y[index2,1]) :- (exp(gam:*y[index2,1]) :-1):/gam)
				Hg[index2] 	 = -Htemp :* (exp(gam :* y[index2,1]) :* y[index2,1]:^2 :- 2 :* y[index2,1] :* exp(gam :* y[index2,1]):/gam :+ 2:* exp(gam :* y[index2,1]):/(gam:^2) :- 2:/(gam:^2))
				if (haslt) {
					Hx[index2,]  = Hx[index2,] :+ Htemp :* x[index2,xindex1] :* x[index2,xindex2] :* (exp(gam:*y[index2,3]) :-1)
					Hxg[index2,] = Hxg[index2,] :+ Htemp :* x[index2,xindex4] :* (y[index2,3] :* exp(gam:*y[index2,3]) :- (exp(gam:*y[index2,3]) :-1):/gam)
					Hg[index2] 	 = Hg[index2] :+ Htemp :* (exp(gam :* y[index2,3]) :* y[index2,3]:^2 :- 2 :* y[index2,3] :* exp(gam :* y[index2,3]):/gam :+ 2:* exp(gam :* y[index2,3]):/(gam:^2) :- 2:/(gam:^2))
				}
			}
		}

		//build H
		Hxsum 	= quadcolsum(Hx,1)
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

`RM' merlin_logl_gompertz_ml(	`gml' gml)
{
	model 		= gml.model
	y 			= merlin_util_depvar(gml)
	Nobs		= merlin_get_nobs(gml,model)
	hast		= gml.istimedep[model,1]
	haslt		= gml.hasltrunc[model]
	hasbh		= gml.hasbh[model,1]
	logl 		= J(Nobs,gml.ndim[gml.Nrelevels],0)

	//core
	xb			= merlin_util_xzb(gml)
	gam 		= merlin_util_dap(gml,1)
	
	//===========================================================================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 			= merlin_get_surv_index(gml)
			logl[index1,] 	= xb[index1,] :+ gam :* y[index1,1]
			if (hasbh) logl[index1,] = log(exp(logl[index1,]) :+ merlin_util_bhazard(gml))
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
					chq2[,q] = exp(merlin_util_xzb(gml,qp2[,q]) :+ gam :* qp2[,q]) 
				}
				if (!haslt) logl[index2,] = logl[index2,] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
				else 	    logl[index2,] = logl[index2,] :- (y[index2,1]:-y[index2,3]):/2 :* (chq2 * gq[,2]) 
			}
			else {
				logl[index2,] = logl[index2,] :- exp(xb[index2,]) :* (1/gam) :*(exp(gam:*y[index2,1]):-1)
				if (haslt) logl[index2,] = logl[index2,] :+ exp(xb[index2,]) :* (1/gam) :*(exp(gam:*y[index2,3]):-1)
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

`RM' merlin_gompertz_logh(`gml' gml, `RC' t)
{
	gam 	= merlin_util_dap(gml,1)
	logh 	= merlin_util_xzb(gml,t) :+ gam :* t
	if (gml.hasbh[gml.model,1]) {
		logh = log(exp(logh) :+ asarray(gml.bhazards,gml.model)[asarray(gml.surv_index,(gml.model,1))])
	}
	return(logh)
}

`RM' merlin_gompertz_ch(`gml' gml,`RC' t, | `RC' t0)
{
	gam = merlin_util_dap(gml,1)
	
	if (gml.NI[gml.model]) {
		nobs 	= rows(t)
		ch 		= J(nobs,1,0)
		Ngq 	= gml.chip
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
		loghazq = gam :* qp :+ log(t:/2 :* J(nobs,1,gq[,2]'))
		
		if (args()==2) { 
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
	else return(exp(merlin_util_xzb(gml,t)) :* (1/gam) :*(exp(gam:*t):-1))
}

`RM' merlin_gompertz_cdf(`gml' gml,`RC' t)
{
	return(1:-exp(-merlin_gompertz_ch(gml,t)))
}

`RM' merlin_gompertz_s(`gml' gml, `RC' t)
{
	return(exp(-merlin_gompertz_ch(gml,t))) 	
}

`RM' merlin_gompertz_expval(`gml' gml, `RC' t)
{
	
}

end
