*! version 2.0.0  03mar2024

local gml 	struct merlin_struct scalar
local pgml	pointer(struct merlin_struct scalar) scalar
local TR 	transmorphic
local RS 	real scalar
local RC 	real colvector
local SS 	string scalar
local PS 	pointer scalar
local RR 	real rowvector
local RM 	real matrix
local PC 	pointer colvector
local PM 	pointer matrix
local SC 	string colvector

version 17

mata:

real matrix merlin_stexcess_logh(`gml' gml, `RM' t)
{
	haz_expect = exp(merlin_util_xzb(gml,t))
        haz_excess = exp(merlin_util_xzb_mod(gml,2,t))
	return(log(haz_expect :+ 
                        gml.indicator[merlin_get_index(gml)] :* 
                        haz_excess))
}

real matrix merlin_stexcess_h(`gml' gml, `RM' t)
{
	haz_expect = exp(merlin_util_xzb(gml,t))
        haz_excess = exp(merlin_util_xzb_mod(gml,2,t))
	return(haz_expect :+ 
                        gml.indicator[merlin_get_index(gml)] :* 
                        haz_excess)
}

`RM' merlin_stexcess_ch(`gml' gml, `RC' t, | `RC' t0)
{
	nobs 	= rows(t)
	ch 	= J(nobs,1,0)
	Ngq 	= gml.chip

	gq 	= merlin_gq(Ngq,"legendre")
	qp	= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
	qw	= t :/ 2 :* J(nobs,1,gq[,2]')
	if (args()==2) {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ (exp(merlin_util_xzb(gml,qp[,q])) :+
				exp(merlin_util_xzb_mod(gml,2,qp[,q]))) :* 
				qw[,q]
		}
	}
	else {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ (exp(merlin_util_xzb(gml,qp[,q],t0)) :+
				exp(merlin_util_xzb_mod(gml,2,qp[,q],t0))) :* 
				qw[,q]
		}
	}
	return(ch)
}

real matrix merlin_stexcess_logl(`gml' gml , | `RM' G, `RM' H)
{
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml,model)
	haslt	= gml.hasltrunc[model]
	logl 	= J(Nobs,1,0)

	// exact events
	gml.survind = 1
	if (merlin_get_nobs(gml,model)) {
		index1 	= merlin_get_surv_index(gml)
		haz_expect = exp(merlin_util_xzb(gml))
		haz_excess = gml.indicator[index1] :* 	///
			exp(merlin_util_xzb_mod(gml,2))
		haz = haz_expect :+ haz_excess
		logl[index1] = log(haz)
	}

	// right censoring
	gml.survind = 2
	Nobs2 	= merlin_get_nobs(gml)
	if (Nobs2) {
		index2 	= merlin_get_surv_index(gml)
		Ngq 	= gml.chip
		chq1 = chq2 = J(Nobs2,Ngq,0)
		gq 	= merlin_gq(Ngq,"legendre")
		if (!haslt) {
			qp2 = y[index2,1] :/ 2 :* J(Nobs2,1,gq[,1]') 	///
				:+ y[index2,1]:/2
		}
		else {
			qp2 = (y[index2,1] :- y[index2,3]) :/ 2 	///
				:* J(Nobs2,1,gq[,1]') :+ (y[index2,1] 	///
				:+ y[index2,3]):/2
		}
		for (q=1;q<=Ngq;q++) {
			chq1[,q] = exp(merlin_util_xzb(gml,qp2[,q])) 
			chq2[,q] = gml.indicator[index2] :* 		///
				exp(merlin_util_xzb_mod(gml,2,qp2[,q]))
		}
		if (!haslt) {
			logl[index2] = logl[index2] :- y[index2,1]:/2 	///
				:* ((chq1 :+ chq2) * gq[,2]) 
		}
		else {
			logl[index2] = logl[index2] :- 			///
				(y[index2,1]:-y[index2,3]):/2 		///
				:* ((chq1:+chq2) * gq[,2]) 
		}
	}	
		
	if (gml.todo==0) return(logl)
	
	//indexes for eqn 1
	model   = 1
	NHbs 	= asarray(gml.NHbs,model)
	sindex1 = (1..NHbs[1]) :+ gml.skip[model]
	//indexes for eqn 2
	model   = 2
	NHbs 	= asarray(gml.NHbs,model)
	sindex2 = (1..NHbs[1]) :+ gml.skip[model]
	
	gml.model = 1
	gml.survind = 0
	x1  = merlin_util_xz(gml)
	gml.model = 2
	gml.survind = 0
	x2  = merlin_util_xz(gml)
	
	G[index1,sindex1] = haz_expect :* x1[index1,]
	G[index1,sindex2] = haz_excess :* x2[index1,]
	G[index1,] = G[index1,] :/ haz

	dchq1 = dchq2 = J(Nobs2,1,0)
	if (!haslt) qw2  = y[index2,1]:/2 :* J(Nobs2,1,gq[,2]')
	else {
		qw2  = (y[index2,1]:-y[index2,3]):/2 :* J(Nobs2,1,gq[,2]')
	}
	for (q=1;q<=Ngq;q++) {
		gml.model = 1
		dchq1 = dchq1 :+ chq1[,q] :* 	///
				merlin_util_xz(gml,qp2[,q]) :* qw2[,q]
		gml.model = 2
		dchq2 = dchq2 :+ chq2[,q] :* 	///
			merlin_util_xz(gml,qp2[,q]) :* qw2[,q]
	}
	G[index2,sindex1] = G[index2,sindex1] :- dchq1
	G[index2,sindex2] = G[index2,sindex2] :- dchq2
	
	if (gml.todo==1) return(logl)
	
	return(logl)		
}

end
