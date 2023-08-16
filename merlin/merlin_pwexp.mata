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

`RM' merlin_logl_pwexp(`gml' gml, `RM' G, `RM' H)
{
	model 		= gml.model
	y 		= merlin_util_depvar(gml)
	Nobs		= merlin_get_nobs(gml,model)
	hast		= gml.istimedep[model,1]
	haslt		= gml.hasltrunc[model]
	hasbh		= gml.hasbh[model,]
	logl 		= J(Nobs,1,0)

	//core
	baseb 		= asarray(gml.distancb,(gml.model,1))
	cuts  		= asarray(gml.distancb,(gml.model,2))
	Ncuts 		= asarray(gml.distancb,(gml.model,3))
	ivindex		= asarray(gml.distancb,(gml.model,4))
	xb		= merlin_util_xzb(gml)

	//===========================================================================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 			= merlin_get_surv_index(gml)
			logl[index1,] 	= xb[index1,] :+ baseb[ivindex[index1]]
			if (hasbh[1]) logl[index1,] = log(exp(logl[index1,]) :+ merlin_util_bhazard(gml))
			if (hasbh[2]) logl[index1,] = log(exp(logl[index1,]) :* merlin_util_bhazard(gml))
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
					chq2[,q] = exp(merlin_util_xzb(gml,qp2[,q]) :+ merlin_pwexp_baselogh(gml,qp2[,q],cuts,Ncuts,baseb))
				}
				if (!haslt) logl[index2] = logl[index2] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
				else 		logl[index2] = logl[index2] :- (y[index2,1]:-y[index2,3]):/2 :* (chq2 * gq[,2]) 
			}
			else {
				expxb2 			= exp(xb[index2,])
				basech2 		= merlin_pwexp_basech(gml,y[index2,1],cuts,Ncuts,baseb)
				logl[index2] 	= logl[index2] :- expxb2 :* basech2
				if (haslt) {
					gml.survind		= 4
					index4			= merlin_get_surv_index(gml)
					expxb4			= exp(xb[index4])
					basech20 		= merlin_pwexp_basech(gml,y[index4,3],cuts,Ncuts,baseb)
					logl[index2] 	= logl[index2] :+ expxb4 :* basech20
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
					chq3[,q] = exp(merlin_util_xzb(gml,qp3[,q]) :+ merlin_pwexp_baselogh(gml,qp3[,q],cuts,Ncuts,baseb))
				}
				ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
				logl[index3,] = 1:-exp(-ch3)
			}
			else {
				expxb3	= exp(xb[index3,])
				basech3 = merlin_pwexp_basech(gml,y[index3,1],cuts,Ncuts,baseb)
				logl[index3,] = -exp(-expxb3 :* basech3)
			}
			
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
					chq5[,q] = exp(merlin_util_xzb(gml,qp5[,q]) :+ merlin_pwexp_baselogh(gml,qp5[,q],cuts,Ncuts,baseb))
				}
				ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
				logl[index5,] = logl[index5,] :- (1:-exp(-ch5))
			}
			else {
				expxb5 = exp(xb[index5,])
				basech5 = merlin_pwexp_basech(gml,y[index5,tind],cuts,Ncuts,baseb)
				logl[index5,] = logl[index5,] :+ exp(-expxb5 :* basech5)
			}
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
			G[index1,sindex2] = designmatrix2(ivindex[index1],Ncuts+1)
		}  

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (merlin_get_nobs(gml,model)) {
			if (hast) {
				dch2 = dch22 = J(Nobs2,1,0)
				if (!haslt) qw2 = y[index2,1]:/2 :* J(Nobs2,1,gq[,2]')
				else		qw2 = (y[index2,1]:-y[index2,3]):/2 :* J(Nobs2,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dch2 	= dch2  :+ chq2[,q] :* merlin_util_xz(gml,qp2[,q]) :* qw2[,q]
					dch22 	= dch22 :+ exp(merlin_util_xzb(gml,qp2[,q])) :* designmatrix2(merlin_pwexp_time_interval(gml,qp2[,q],cuts,Ncuts),Ncuts+1) :* qw2[,q]
				}
				dch22 = dch22 :* J(Nobs2,1,exp(baseb)')
				G[index2,sindex1] = G[index2,sindex1] :- dch2
				G[index2,sindex2] = G[index2,sindex2] :- dch22 
			}
			else {
				G[index2,sindex1] 	= G[index2,sindex1] :- expxb2 :* basech2 :* x[index2,]
				dbch20 				= merlin_pwexp_basech_db(gml,y[index2,1],cuts,Ncuts,baseb)
				G[index2,sindex2] 	= G[index2,sindex2] :- expxb2 :* dbch20
				if (haslt) {
					G[index4,sindex1] = G[index4,sindex1] :+ expxb4 :* basech20 :* x[index4,]
					dbch200	= merlin_pwexp_basech_db(gml,y[index4,3],cuts,Ncuts,baseb)
					G[index4,sindex2] = G[index4,sindex2] :+ expxb4 :* dbch200
				}
			}
		}

		//interval censoring -> cdf function
		gml.survind = 3
		if (Nobs3) {	
			//exit times
			if (hast) {
				dchq3 	= dch32 =  J(Nobs3,1,0)
				qw3		= y[index3,1] :/ 2 :* J(Nobs3,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dchq3 = dchq3 :+ chq3[,q] :* qw3[,q] :* merlin_util_xz(gml,qp3[,q])
					dch32 = dch32 :+ exp(merlin_util_xzb(gml,qp3[,q])) :* designmatrix2(merlin_pwexp_time_interval(gml,qp3[,q],cuts,Ncuts),Ncuts+1) :* qw3[,q]
				}
				dch32 = dch32 :* J(Nobs3,1,exp(baseb)')
				G[index3,sindex1] = exp(-ch3) :* dchq3
				G[index3,sindex2] = exp(-ch3) :* dch32
			}
			else {
				G[index3,sindex1] = exp(-expxb3 :* basech3) :* expxb3 :* basech3 :* merlin_util_xz(gml)
				G[index3,sindex2] = exp(-expxb3 :* basech3) :* expxb3 :* merlin_pwexp_basech_db(gml,y[index3,1],cuts,Ncuts,baseb)
			}
			
			//entry times
			gml.survind = 5
			if (hast) {
				dchq5 	= dch52 = J(Nobs5,1,0)
				qw5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,2]')
				for (q=1;q<=Ngq;q++) {
					dchq5 = dchq5 :+ chq5[,q] :* qw5[,q] :* merlin_util_xz(gml,qp5[,q])
					dch52 = dch52 :+ exp(merlin_util_xzb(gml,qp5[,q])) :* designmatrix2(merlin_pwexp_time_interval(gml,qp5[,q],cuts,Ncuts),Ncuts+1) :* qw5[,q]
				}
				dch52 = dch52 :* J(Nobs3,1,exp(baseb)')
				G[index5,sindex1] = G[index5,sindex1] :- exp(-ch5) :* dchq5
				G[index5,sindex2] = G[index5,sindex2] :- exp(-ch5) :* dch52
			}
			else {
				G[index5,sindex1] = G[index5,sindex1] :- exp(-expxb5 :* basech5) :* expxb5 :* basech5 :* merlin_util_xz(gml)
				G[index5,sindex2] = G[index5,sindex2] :- exp(-expxb5 :* basech5) :* expxb5 :* merlin_pwexp_basech_db(gml,y[index5,tind],cuts,Ncuts,baseb)
			}
			sindex 				= sindex1,sindex2
			G[index3,sindex] 	= G[index3,sindex] :/ exp(logl[index3,])
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
		xindex5 = xindex6 = xindex7 = J(1,0,.)
		diagind = 1
		for (i=1; i<=NHbs[2]; i++) {
			refind = 1
			while (refind<=i) {
				xindex5 = xindex5,i
				xindex6 = xindex6,refind
				if (refind==i) xindex7 = xindex7,diagind
				refind++
				diagind++
			}
		}
		Hg = J(Nobs,cols(xindex5),0)

		//exactly observed events -> hazard function
		//all 0

		//exactly observed events and/or right censoring -> survival function
		gml.survind = 2
		if (Nobs2) {
			//xb
			if (hast) {
				d2ch2 = d2ch22 = J(Nobs2,1,0)
				for (q=1;q<=Ngq;q++) {
					x = merlin_util_xz(gml,qp2[,q])
					d2ch2 	= d2ch2 :+ chq2[,q] :* x[,xindex1] :* x[,xindex2] :* qw2[,q]
					x2 		= designmatrix2(merlin_pwexp_time_interval(gml,qp2[,q],cuts,Ncuts),Ncuts+1) :* J(Nobs2,1,exp(baseb)')
					d2ch22 	= d2ch22 :+ x[,xindex4] :* exp(merlin_util_xzb(gml,qp2[,q])) :* x2[,xindex3] :* qw2[,q]
				}
				Hx[index2,]  		= - d2ch2
				Hxg[index2,] 		= - d2ch22
				Hg[index2,xindex7]	= - dch22

			}
			else {
				Hx[index2,]  		= -expxb2 :* basech2 :* x[index2,xindex1] :* x[index2,xindex2]
				Hxg[index2,] 		= -expxb2 :* x[index2,xindex4] :* dbch20[,xindex3]
				Hg[index2,xindex7]	= -expxb2 :* dbch20
				if (haslt) {
					Hx[index4,]  		= Hx[index4,] :+ expxb4 :* basech20 :* x[index4,xindex1] :* x[index4,xindex2]
					Hxg[index4,] 		= Hxg[index4,] :+ expxb4 :* x[index4,xindex4] :* dbch200[,xindex3]
					Hg[index4,xindex7]	= Hg[index4,xindex7] :+ expxb4 :* dbch200
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
		Hindex2 = NHbs[1] + NHbs[2] + gml.skip[model] 
		H[|Hindex1,Hindex1\Hindex2,Hindex2|] = Hm
		
		return(logl)
}

`RM' merlin_logl_pwexp_ml(`gml' gml)
{
	model 		= gml.model
	y 		= merlin_util_depvar(gml)
	Nobs		= merlin_get_nobs(gml,model)
	hast		= gml.istimedep[model,1]
	hasbh		= gml.hasbh[model,1]
	logl 		= J(Nobs,gml.ndim[gml.Nrelevels],0)

	//core
	baseb 		= asarray(gml.distancb,(gml.model,1))
	cuts  		= asarray(gml.distancb,(gml.model,2))
	Ncuts 		= asarray(gml.distancb,(gml.model,3))
	ivindex		= asarray(gml.distancb,(gml.model,4))
	xb		= merlin_util_xzb(gml)

	//===========================================================================================================================//
	// log likelihood

		//exactly observed events -> hazard function
		gml.survind = 1
		if (merlin_get_nobs(gml,model)) {
			index1 		= merlin_get_surv_index(gml)
			logl[index1,] 	= xb[index1,] :+ baseb[ivindex[index1]]
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
				gq 	= merlin_gq(Ngq,"legendre")
				qp2	= y[index2,1] :/ 2 :* J(Nobs2,1,gq[,1]') :+ y[index2,1]:/2
				for (q=1;q<=Ngq;q++) {
					chq2[,q] = exp(merlin_util_xzb(gml,qp2[,q]) :+ merlin_pwexp_baselogh(gml,qp2[,q],cuts,Ncuts,baseb))
				}
				logl[index2,] = logl[index2,] :- y[index2,1]:/2 :* (chq2 * gq[,2]) 
			}
			else {
				expxb2 		= exp(xb[index2,])
				basech2 	= merlin_pwexp_basech(gml,y[index2,1],cuts,Ncuts,baseb)
				logl[index2,] 	= logl[index2,] :- expxb2 :* basech2
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
					chq3[,q] = exp(merlin_util_xzb(gml,qp3[,q]) :+ merlin_pwexp_baselogh(gml,qp3[,q],cuts,Ncuts,baseb))
				}
				ch3 = y[index3,1] :/ 2 :* (chq3 * gq[,2])
				logl[index3,] = 1:-exp(-ch3)
			}
			else {
				expxb3	= exp(xb[index3,])
				basech3 = merlin_pwexp_basech(gml,y[index3,1],cuts,Ncuts,baseb)
				logl[index3,] = -exp(-expxb3 :* basech3)
			}
			
			//entry times
			gml.survind = 5
			if (gml.hasltrunc[model]) 	tind = 4
			else 				tind = 3
			index5 	= merlin_get_surv_index(gml)
			if (hast) {
				Nobs5 	= merlin_get_nobs(gml)
				chq5 	= J(Nobs5,Ngq,.)
				qp5		= y[index5,tind] :/ 2 :* J(Nobs5,1,gq[,1]') :+ y[index5,tind]:/2
				for (q=1;q<=Ngq;q++) {
					chq5[,q] = exp(merlin_util_xzb(gml,qp5[,q]) :+ merlin_pwexp_baselogh(gml,qp5[,q],cuts,Ncuts,baseb))
				}
				ch5 = y[index5,tind] :/ 2 :* (chq5 * gq[,2])
				logl[index5,] = logl[index5,] :- (1:-exp(-ch5))
			}
			else {
				expxb5 = exp(xb[index5,])
				basech5 = merlin_pwexp_basech(gml,y[index5,tind],cuts,Ncuts,baseb)
				logl[index5,] = logl[index5,] :+ exp(-expxb5 :* basech5)
			}
			//logL
			logl[index3,] = log(logl[index3,])
		}

	if (gml.todo==0) return(logl)
	
}

`RM' merlin_pwexp_logh(`gml' gml,`RC' t)
{
	logh  = merlin_util_xzb(gml,t)
	baseb = asarray(gml.distancb,(gml.model,1))
	cuts  = asarray(gml.distancb,(gml.model,2))
	Ncuts = cols(cuts)
	logh  = logh :+ merlin_pwexp_baselogh(gml,t,cuts,Ncuts,baseb)
	
	if (gml.hasbh[gml.model,1]) {
		logh = log(exp(logh) :+ asarray(gml.bhazards,gml.model)[asarray(gml.surv_index,(gml.model,1))])
	}
	return(logh)
}

`RM' merlin_pwexp_time_interval(`gml' gml,`RC' t, `RR' cuts, `RS' Ncuts)
{
	baseindex = J(merlin_get_nobs(gml),1,.)
	
	//first interval
	ix 	= selectindex(t:<cuts[1])
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) baseindex[ix] = J(nrc[1],1,1)
	//mid
	if (Ncuts>1) {
		for (i=1;i<Ncuts;i++) {
			ix 	= selectindex((cuts[i]:<=t) :& (t:<cuts[i+1]))
			nrc = rows(ix)\cols(ix)
			if (nrc[1] & nrc[2]) baseindex[ix] = J(nrc[1],1,i+1)
		} 
	}
	//last interval
	ix 	= selectindex(cuts[Ncuts]:<=t)
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) baseindex[ix] = J(nrc[1],1,Ncuts+1)
	
	return(baseindex)
}

`RM' merlin_pwexp_baselogh(`gml' gml,`RC' t, `RR' cuts, `RS' Ncuts, `RC' baseb)
{
	return(baseb[merlin_pwexp_time_interval(gml,t,cuts,Ncuts)])
}

`RM' merlin_pwexp_ch(`gml' gml, `RC' t,| `RC' t0)
{
	hast0 = args()==3
	baseb = asarray(gml.distancb,(gml.model,1))
	cuts  = asarray(gml.distancb,(gml.model,2))
	Ncuts = cols(cuts)
	
	if (gml.NI[gml.model]) {
		nobs 	= rows(t)
		ch 		= J(nobs,1,0)
		Ngq 	= 30
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
		loghazq = log(t:/2 :* J(nobs,1,gq[,2]'))
		if (hast0) {
			for (q=1;q<=Ngq;q++) {
				ch = ch :+ exp(merlin_util_xzb(gml,qp[,q],t0) :+ merlin_pwexp_baselogh(gml,qp[,q],cuts,Ncuts,baseb) :+ loghazq[,q])
			}
		}
		else {
			for (q=1;q<=Ngq;q++) {
				ch = ch :+ exp(merlin_util_xzb(gml,qp[,q]) :+ merlin_pwexp_baselogh(gml,qp[,q],cuts,Ncuts,baseb) :+ loghazq[,q])
			}
		}
		return(ch)
		
	}
	else {
		if (hast0) 	return(exp(merlin_util_xzb(gml,t,t0)) :* merlin_pwexp_basech(gml,t,cuts,Ncuts,baseb))
		else 		return(exp(merlin_util_xzb(gml,t)) :* merlin_pwexp_basech(gml,t,cuts,Ncuts,baseb))
	}
}

`RM' merlin_pwexp_basech(`gml' gml,`RC' t, `RR' cuts, `RS' Ncuts, `RC' baseb)
{
	ch0 = J(merlin_get_nobs(gml),1,.)
	expbaseb = exp(baseb)

	//calc. base ch at each knot
	basech = cuts
	if (Ncuts>1) 	steps = cuts :- (0,cuts[1..(Ncuts-1)])
	else 			steps = cuts

	basech = steps :* expbaseb[1::Ncuts]'
	basech = runningsum(basech)

	//first interval
	ix 	= selectindex(t:<cuts[1])
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) ch0[ix] = expbaseb[1] :* t[ix]
	//mid
	if (Ncuts>1) {
		for (i=1;i<Ncuts;i++) {
			ix 	= selectindex(cuts[i]:<=t :& t:<cuts[i+1])
			nrc = rows(ix)\cols(ix)
			if (nrc[1] & nrc[2]) ch0[ix] = basech[i] :+ (t[ix] :- cuts[i]) :* expbaseb[i+1]
		} 
	}
	//last interval
	ix 	= selectindex(cuts[Ncuts]:<=t)
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) ch0[ix] = basech[Ncuts] :+ (t[ix] :- cuts[Ncuts]) :* expbaseb[Ncuts+1]
	
	return(ch0)
}

`RM' merlin_pwexp_basech_db(`gml' gml,`RC' t, `RR' cuts, `RS' Ncuts, `RC' baseb)
{
	dbch0 		= J(merlin_get_nobs(gml),Ncuts+1,0)
	expbaseb 	= exp(baseb)

	//calc. base ch contribution in each interval
	if (Ncuts>1) 	steps = cuts :- (0,cuts[1..(Ncuts-1)])
	else 			steps = cuts
	intbasech 	= steps :* expbaseb[1::Ncuts]'

	//first interval
	ix 	= selectindex(t:<cuts[1])
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) dbch0[ix,1] = expbaseb[1] :* t[ix]
	//mid
	if (Ncuts>1) {
		for (i=1;i<Ncuts;i++) {
			ix 	= selectindex(cuts[i]:<=t :& t:<cuts[i+1])
			nrc = rows(ix)\cols(ix)
			if (nrc[1] & nrc[2]) {
				for (j=1;j<=i;j++) dbch0[ix,j] = J(nrc[1],1,intbasech[j])
				dbch0[ix,i+1] = (t[ix] :- cuts[i]) :* expbaseb[i+1]
			}
		} 
	}
	//last interval
	ix 	= selectindex(cuts[Ncuts]:<=t)
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) {
		for (j=1;j<=Ncuts;j++) dbch0[ix,j] = J(nrc[1],1,intbasech[j])
		dbch0[ix,Ncuts+1] = (t[ix] :- cuts[Ncuts]) :* expbaseb[Ncuts+1]
	}
	
	return(dbch0)
}

`RM' merlin_pwexp_cdf(`gml' gml,`RC' t)
{
	return(1:-exp(-merlin_pwexp_ch(gml,t)))
}

`RM' merlin_pwexp_s(`gml' gml, | `RC' t)
{
	return(exp(-merlin_pwexp_ch(gml,t))) 
}

end
