
local GMLS 		struct merlin_struct scalar
local pGMLS		pointer(struct merlin_struct scalar) scalar
local RS 		real scalar
local SS 		string scalar
local PS 		pointer scalar
local RM 		real matrix
local SM		string matrix
local PC 		pointer colvector
local SC 		string colvector
local SR		string rowvector
local TR		transmorphic
local RC		real colvector
local RR		real rowvector
local PM		pointer matrix

version 14.2

mata:

`RC' merlin_get_index(`GMLS' gml)
{
	if (gml.imputing) 	index = asarray(gml.ImputeIndex,gml.modtoind)
	else {
		if (gml.survind==-99) return(1::gml.N)
		index = asarray(gml.xbindex,gml.modtoind)
	}

	if (!gml.survind) 	return(index)
	else 		        return(index[merlin_get_surv_index(gml)])
}

//get model specific index for survival outcomes
//can't if else, as it could be 0
`RC' merlin_get_surv_index(`GMLS' gml)
{
	if (gml.imputing) {
		if 	(gml.survind==1) return(asarray(gml.surv_index_imp,(gml.modtoind,1)))
		else if (gml.survind==2) return(asarray(gml.surv_index_imp,(gml.modtoind,2)))
		else if (gml.survind==3) return(asarray(gml.surv_index_imp,(gml.modtoind,3)))
		else if (gml.survind==5) return(asarray(gml.surv_index_imp,(gml.modtoind,5)))
		else if (gml.survind==4) return(asarray(gml.surv_index_imp,(gml.modtoind,4)))
		else if (gml.survind==6) return(asarray(gml.surv_index_imp,(gml.modtoind,6)))
	}
	else {		
		if 	(gml.survind==1) return(asarray(gml.surv_index,(gml.modtoind,1)))
		else if (gml.survind==2) return(asarray(gml.surv_index,(gml.modtoind,2)))
		else if (gml.survind==3) return(asarray(gml.surv_index,(gml.modtoind,3)))
		else if (gml.survind==5) return(asarray(gml.surv_index,(gml.modtoind,5)))
		else if (gml.survind==4) return(asarray(gml.surv_index,(gml.modtoind,4)))
		else if (gml.survind==6) return(asarray(gml.surv_index,(gml.modtoind,6)))
                else if (gml.survind==7) return(asarray(gml.surv_index,(gml.modtoind,7)))
		else 	return(asarray(gml.cox_index,(gml.modtoind,2,gml.survind)))
	}
}

`RC' merlin_get_adpanelindex(`GMLS' gml, `RS' lev)
{
	index = asarray(gml.adpanelindexes,(lev,gml.modtoind))
	if (gml.survind==0) return(index)
	else 		    return(index[merlin_get_surv_index(gml)])
}

`RS' merlin_get_nobs(`GMLS' gml,| `RS' mod)
{
	if (args()==1) mod = gml.modtoind 
	
	if (gml.imputing) {
		if 	(gml.survind==0) 	return(rows(asarray(gml.ImputeIndex,mod)))		//poor
		else if (gml.survind==1) 	return(gml.NsurvImp[mod,1])
		else if (gml.survind==2) 	return(gml.NsurvImp[mod,2])
		else if (gml.survind==3) 	return(gml.NsurvImp[mod,3])
		else if (gml.survind==5) 	return(gml.NsurvImp[mod,5])
		else if (gml.survind==4) 	return(gml.NsurvImp[mod,4])
		else 					 	return(gml.NsurvImp[mod,6])
	}
	else {
		if 	(gml.survind==0) 	return(gml.Nobs[gml.Nlevels,mod])
		else if (gml.survind==1) 	return(gml.Nsurv[mod,1])
		else if (gml.survind==2) 	return(gml.Nsurv[mod,2])
		else if (gml.survind==3) 	return(gml.Nsurv[mod,3])
		else if (gml.survind==5) 	return(gml.Nsurv[mod,5])
		else if (gml.survind==4) 	return(gml.Nsurv[mod,4])
		else if (gml.survind==6) 	return(gml.Nsurv[mod,6])
                else if (gml.survind==7) 	return(gml.Nsurv[mod,7])
		else return(asarray(gml.cox_index,(gml.modtoind,1,gml.survind)))
	}
}

`RM' merlin_xz_var(`GMLS' gml, `RS' i, `RS' j, | `RC' t, `RC' t0)
{	
	if (gml.imputing)	return(gml.ImputedValue)
	else {
		hast	= args()==5
		hast0	= args()==6
		info 	= asarray(gml.elinfo,(gml.model,i,j))
		istvars	= asarray(info,1)
		
		if (hast) {
	 		if 	(istvars[1]) 	var = t
			else if (istvars[2]) 	var = merlin_util_depvar(gml)[,3]
			else 			var = asarray(info,2)[merlin_get_index(gml)]
		}
		else if (hast0) {
			if 	(istvars[1]) 	var = t
			else if (istvars[2]) 	var = t0
			else 			var = asarray(info,2)[merlin_get_index(gml)]
		}
		else {
			if 	(istvars[1]) 	var = asarray(gml.timevars,gml.model)[merlin_get_index(gml)]
			else if (istvars[2]) 	var = merlin_util_depvar(gml)[,3]
			else 			var = asarray(info,2)[merlin_get_index(gml)]
		}
	
		return(var)
	}
}

void merlin_xz_var_update(`GMLS' gml, `RS' c, `RS' el, `SS' v, `RC' val)
{
	if (asarray(gml.elvarname,(gml.model,c,el))==v) {
		if (rows(val)==1) 	asarray(gml.elinfo,(gml.model,c,el),J(gml.N,1,val))
		else 			asarray(gml.elinfo,(gml.model,c,el),val)
	}
}


`RM' merlin_xz_b(`GMLS' gml, `RS' i, `RS' j, `RS' re)
{
	elinfo 	= asarray(gml.elinfo,(gml.model,i,j))	
	lev	= asarray(elinfo,1)
	reindex	= asarray(elinfo,2)[re]

	if (gml.adapt[lev]) {
		if (gml.ltflag) {
			if (lev==gml.Nrelevels) return(asarray(gml.aghip2_lt,(lev,reindex))[merlin_get_adpanelindex(gml,lev),])
			else 			return(asarray(gml.aghip2_lt,(lev,reindex))[merlin_get_adpanelindex(gml,lev),gml.qind[,lev+1]])
		}
		else {
			if (lev==gml.Nrelevels) return(asarray(gml.aghip2,(lev,reindex))[merlin_get_adpanelindex(gml,lev),])
			else 			return(asarray(gml.aghip2,(lev,reindex))[merlin_get_adpanelindex(gml,lev),gml.qind[,lev+1]])		
		}
	}
	else {
		if (lev==gml.Nrelevels)	return(J(merlin_get_nobs(gml),1,asarray(gml.b,lev)[reindex,]))
		else 			return(J(merlin_get_nobs(gml),1,asarray(gml.b,lev)[reindex,gml.qind[,lev+1]]))
	}
}

`RM' merlin_xz_blups(`GMLS' gml, `RS' i, `RS' j, `RS' re)
{
	elinfo 	= asarray(gml.elinfo,(gml.model,i,j))	
	lev	= asarray(elinfo,1)
	reindex	= asarray(elinfo,2)[re]
	if (lev==gml.Nrelevels) return(asarray(gml.blups,lev)[merlin_get_adpanelindex(gml,lev),reindex])
// 	else 			return(asarray(gml.blups,(lev,info[2]))[merlin_get_adpanelindex(gml,lev),gml.qind[,lev+1]])		
}

`RM' merlin_xz_specblups(`GMLS' gml, `RS' i, `RS' j, `RS' re)
{
	elinfo 	= asarray(gml.elinfo,(gml.model,i,j))	
	lev	= asarray(elinfo,1)
	reindex	= asarray(elinfo,2)[re]
	if (lev==gml.Nrelevels) return(asarray(gml.blups,lev)[reindex])
// 	else 			return(asarray(gml.blups,(lev,info[2]))[merlin_get_adpanelindex(gml,lev),gml.qind[,lev+1]])		
}

`RM' merlin_xz_mf(`GMLS' gml, `RS' i, `RS' j, | `RC' t)
{
	if (args()==3)	return((*asarray(gml.elinfo,(gml.model,i,j))[1])(gml,asarray(gml.timevars,gml.model)[merlin_get_index(gml)]))
	else 		return((*asarray(gml.elinfo,(gml.model,i,j))[1])(gml,t))
}

`RM' merlin_xz_fp(`GMLS' gml, `RS' c, `RS' el, `RS' deriv, | `RC' t, `RC' t0)
{
	hast	= args()==5
	hast0	= args()==6
	fpinfo 	= asarray(gml.elinfo,(gml.model,c,el))
	istvars	= asarray(fpinfo,1)
	index	= merlin_get_index(gml)
	powers 	= asarray(fpinfo,3)		

	if (hast) {
		if 	(istvars[1]) 	fpvar = t
		else if (istvars[2]) 	fpvar = merlin_util_depvar(gml)[,3]
		else 			fpvar = asarray(fpinfo,2)[index]
	}
	else if (hast0) {
		if 	(istvars[1]) 	fpvar = t
		else if (istvars[2]) 	fpvar = t0
		else 			fpvar = asarray(fpinfo,2)[index]
	}
	else {
		if 	(istvars[1]) 	fpvar = asarray(gml.timevars,gml.model)[index]
		else if (istvars[2]) 	fpvar = merlin_util_depvar(gml)[,3]
		else 			fpvar = asarray(fpinfo,2)[index]
	}

	hasoffset 	= asarray(fpinfo,5) //offset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	fpvar = fpvar :+ t0
			else 		fpvar = fpvar :+ merlin_util_depvar(gml)[,3]
		} 
		else fpvar = fpvar :+ asarray(fpinfo,7)[index]
	}
	hasoffset 	= asarray(fpinfo,6) //moffset

	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	fpvar = fpvar :- t0
			else 		fpvar = fpvar :- merlin_util_depvar(gml)[,3]
		} 
		else fpvar = fpvar :+ asarray(fpinfo,8)[index]
	}	
	
	return(merlin_fp(fpvar,powers,deriv))
}

`RM' merlin_xz_rcs(`GMLS' gml, `RS' c, `RS' el, `RS' deriv, | `RC' t, `RC' t0)
{
	hast		= args()==5
	hast0		= args()==6
	rcsinfo 	= asarray(gml.elinfo,(gml.model,c,el))

	istvars 	= merlin_util_istimevar(gml,c,el) 	//asarray(rcsinfo,1)
	index		= merlin_get_index(gml)

	if (hast) {
		if 	(istvars[1]) 	rcsvar = t
		else if (istvars[2]) 	rcsvar = merlin_util_depvar(gml)[,3]
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}
	else if (hast0) {
		if 	(istvars[1]) 	rcsvar = t
		else if (istvars[2]) 	rcsvar = t0
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}
	else {
		if 	(istvars[1]) 	rcsvar = asarray(gml.timevars,gml.model)[index]
		else if (istvars[2]) 	rcsvar = merlin_util_depvar(gml)[,3]
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}

	hasoffset = asarray(rcsinfo,7) //offset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	rcsvar = rcsvar :+ t0
			else 		rcsvar = rcsvar :+ merlin_util_depvar(gml)[,3]
		} 
		else rcsvar = rcsvar :+ asarray(rcsinfo,8)[index]
	}
	hasoffset 	= asarray(rcsinfo,9) //moffset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	rcsvar = rcsvar :- t0
			else 		rcsvar = rcsvar :- merlin_util_depvar(gml)[,3]
		} 
		else rcsvar = rcsvar :+ asarray(rcsinfo,10)[index]
	}
	
	islog 	= asarray(rcsinfo,4)
	knots 	= asarray(rcsinfo,3)
	isorth 	= asarray(rcsinfo,5)

	if (isorth) {
		rmat = asarray(rcsinfo,6)
		if (islog)  result = merlin_rcs(log(rcsvar),knots,deriv,rmat)
		else 	    result = merlin_rcs(rcsvar,knots,deriv,rmat)		
	}
	else {
		if (islog)  result = merlin_rcs(log(rcsvar),knots,deriv)
		else 	    result = merlin_rcs(rcsvar,knots,deriv)
	}

	if (deriv & islog)  return(result:/rcsvar)
	else 		    return(result)
	
}

`RM' merlin_xz_pc(`GMLS' gml, `RS' c, `RS' el, | `RC' t, `RC' t0)
{
	hast	= args()==4
	hast0	= args()==5
	info 	= asarray(gml.elinfo,(gml.model,c,el))
	
	istvars = asarray(info,1)
	index	= merlin_get_index(gml)

	if (hast) {
		if 	(istvars[1]) 	var = t
		else if (istvars[2]) 	var = merlin_util_depvar(gml)[,3]
		else 			var = asarray(info,2)[index]
	}
	else if (hast0) {
		if 	(istvars[1]) 	var = t
		else if (istvars[2]) 	var = t0
		else 			var = asarray(info,2)[index]
	}
	else {
		if 	(istvars[1]) 	var = asarray(gml.timevars,gml.model)[index]
		else if (istvars[2]) 	var = merlin_util_depvar(gml)[,3]
		else 			var = asarray(info,2)[index]
	}

	hasoffset 	= asarray(info,4) //offset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	var = var :+ t0
			else 		var = var :+ merlin_util_depvar(gml)[,3]
		} 
		else var = var :+ asarray(info,5)[index]
	}
	hasoffset 	= asarray(info,6) //moffset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	var = var :- t0
			else 		var = var :- merlin_util_depvar(gml)[,3]
		} 
		else var = var :+ asarray(info,7)[index]
	}
	
	knots 	= asarray(info,3)

	return(merlin_pc(var,knots,asarray(info,8)))	
}

`RM' merlin_xz_rcs_EV(`GMLS' gml, `RS' c, `RS' el, | `RC' t)
{
	hast		= args()==4
	rcsinfo 	= asarray(gml.elinfo,(gml.model,c,el))
	
	istvar 		= asarray(rcsinfo,1)
	index		= merlin_get_index(gml)
	
	yindex		= asarray(rcsinfo,12)
	
	if (hast) 	ev = merlin_util_expval_mod(gml,yindex,t)
	else 		ev = merlin_util_expval_mod(gml,yindex)
	
	//now expand into spline basis
	ncols = cols(ev)		//may have random effects
	
	islog 	= asarray(rcsinfo,4)
	if (islog) ev = log(ev)
	
	knots 	= asarray(rcsinfo,3)
	isorth 	= asarray(rcsinfo,5)
	
	if (isorth) {
        
		rmat 	= asarray(rcsinfo,6)
		result 	= merlin_rcs(ev,knots,deriv,rmat)
				
	}
	else {
		result 	= merlin_rcs(ev,knots,deriv)
	}

	return(result)
	
}

`RM' merlin_xz_bs(`GMLS' gml, `RS' c, `RS' el, | `RC' t, `RC' t0)
{
	hast		= args()==4
	hast0 		= args()==5
	rcsinfo 	= asarray(gml.elinfo,(gml.model,c,el))

	istvars 	= asarray(rcsinfo,1)
	index		= merlin_get_index(gml)

	if (hast) {
		if (istvars[1]) 	rcsvar = t
		else if (istvars[2]) 	rcsvar = merlin_util_depvar(gml)[,3]
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}
	else if (hast0) {
		if (istvars[1]) 	rcsvar = t
		else if (istvars[2]) 	rcsvar = t0
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}
	else {
		if (istvars[1]) 	rcsvar = asarray(gml.timevars,gml.model)[index]
		else if (istvars[2]) 	rcsvar = merlin_util_depvar(gml)[,3]
		else 			rcsvar = asarray(rcsinfo,2)[index]
	}

	hasoffset 	= asarray(rcsinfo,6) //offset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	rcsvar = rcsvar :+ t0
			else 		rcsvar = rcsvar :+ merlin_util_depvar(gml)[,3]
		} 
		else rcsvar = rcsvar :+ asarray(rcsinfo,7)[index]
	}
	hasoffset 	= asarray(rcsinfo,11) //moffset
	if (hasoffset[1]) {
		if (hasoffset[2]) {
			if (hast0) 	rcsvar = rcsvar :- t0
			else 		rcsvar = rcsvar :- merlin_util_depvar(gml)[,3]
		} 
		else rcsvar = rcsvar :+ asarray(rcsinfo,12)[index]
	}	
	
	islog 	= asarray(rcsinfo,5)
	knots 	= asarray(rcsinfo,3)
	degree 	= asarray(rcsinfo,4)
	Nbasis 	= asarray(rcsinfo,8)
	hasint	= asarray(rcsinfo,10)
	
	if (islog) 	result = merlin_bs(log(rcsvar),knots,degree,Nbasis,hasint)
	else 		result = merlin_bs(rcsvar,knots,degree,Nbasis,hasint)
	
	return(result)
	
}

`RM' merlin_xz_ps(`GMLS' gml, `RS' c, `RS' el, | `RC' t)
{
	hast	= args()==4
	rcsinfo = asarray(gml.elinfo,(gml.model,c,el))
	
	istvar 	= asarray(rcsinfo,1)
	index	= merlin_get_index(gml)
	
	if (hast) {
		if (istvar) rcsvar = t
		else 		rcsvar = asarray(rcsinfo,2)[index]
	}
	else {
		if (istvar) rcsvar = asarray(gml.timevars,gml.model)[index]
		else 		rcsvar = asarray(rcsinfo,2)[index]
	}
	
	hasoffset 	= asarray(rcsinfo,6)
	if (hasoffset) rcsvar = rcsvar :+ asarray(rcsinfo,7)[index]
	hasmoffset 	= asarray(rcsinfo,11)
	if (hasmoffset) rcsvar = rcsvar :+ asarray(rcsinfo,12)[index]
	
	islog 	= asarray(rcsinfo,5)
	knots 	= asarray(rcsinfo,3)
	degree 	= asarray(rcsinfo,4)
	Nbasis 	= asarray(rcsinfo,8)
	hasint	= asarray(rcsinfo,10)
	
	if (islog) 	result = merlin_ps(log(rcsvar),knots,degree,Nbasis,hasint)
	else 		result = merlin_ps(rcsvar,knots,degree,Nbasis,hasint)
	
	return(result)
	
}

// `RM' merlin_xz_EV(`GMLS' gml, `RS' i, `RS' j, | `RC' t)
// {
// 	k = gml.model
//         //update index for appropriate outcome
// 	info = asarray(gml.elinfo,(gml.model,i,j))
// 	gml.model = k2 = asarray(info,1)
// 	asarray(info,1), asarray(info,2), asarray(info,3)
// 	if (asarray(info,2)) {
// 		nobs = merlin_get_nobs(gml,k)
// 		t = J(nobs,1,asarray(info,3))
// 	}
// 	if (args()==3) 	res =  (*gml.invlinks[k2])(merlin_util_xzb(gml))
// 	else 		res =  (*gml.invlinks[k2])(merlin_util_xzb(gml,t))	
// 	gml.model = k
// 	return(res)
// }

`RM' merlin_xz_iEV(`GMLS' gml, `RS' i, `RS' j, | `RC' t)
{
	k = gml.model
	ngq = 30
	gq = merlin_gq(ngq,"legendre")
	nodes = t :* J(rows(t),1,gq[,1]'):/2 :+ t:/2
	weights = t :* J(rows(t),1,gq[,2]'):/2
	
	gml.model = asarray(gml.elinfo,(gml.model,i,j))			//update index for appropriate outcome
	res = J(rows(t),1,0)
	for (q=1;q<=ngq;q++) {
		res = res :+ weights[,q] :* (*gml.invlinks[gml.model])(merlin_util_xzb(gml,nodes[,q]))	
	}
	gml.model = k
	return(res)
}

`RM' merlin_xz_dt(`GMLS' gml, `RS' i, `RS' j,| `RC' t)
{
	p = asarray(gml.elinfo,(gml.model,i,j))[1]
	
	if (args()==3)	{
		tvar = asarray(gml.timevars,gml.model)[merlin_get_index(gml)]
		h = 1e-8 * (abs(tvar) + 1e-8):/2
		lh = (*p)(tvar :+ h)
		rh = (*p)(tvar :- h)
		return((lh:-rh):/(2:*h))
	}
	else {
		h = 1e-8 * (abs(tvar) + 1e-8):/2
		lh = (*p)(t :+ h)
		rh = (*p)(t :- h)
		return((lh:-rh):/(2:*h))
	}
}

`RM' merlin_xz_fpdt(`GMLS' gml, `RS' i, `RS' j, | `RC' t)
{
	hast			= args()==4
	fpinfo 			= asarray(gml.elinfo,(gml.model,c,el))
	istvar			= asarray(fpinfo,1)
	index			= merlin_get_index(gml)
	
	if (hast) {
		if (istvar) fpvar = t
		else 		fpvar = asarray(fpinfo,2)[index]
	}
	else {
		if (istvar) fpvar = asarray(gml.timevars,gml.modtoind)
		else 		fpvar = asarray(fpinfo,2)[index]
	}
	powers 			= asarray(fpinfo,3)		
	hasanyoffset 	= asarray(fpinfo,5)
	if (hasanyoffset) fpvar = fpvar :+ asarray(fpinfo,6)[index]

	return(merlin_fp(tvar,powers,1))
}

`RM' merlin_xz_d2t(`GMLS' gml, `RS' i, `RS' j,| `RC' t)
{
	p = asarray(gml.elinfo,(gml.model,i,j))
	if (args()==3)	{
		tvar = asarray(gml.timevars,gml.model)[merlin_get_index(gml)]
		h = 1e-8 * (abs(tvar) + 1e-8):/2
		lh = (*p)(tvar :+ h)
		rh = (*p)(tvar :- h)
		mid = 2:*(*p)(tvar)
		return((lh :- mid :+ rh):/(h:^2))
	}
	else {
		h = 1e-8 * (abs(t) + 1e-8):/2
		lh = (*p)(t :+ h)
		rh = (*p)(t :- h)
		mid = 2:*(*p)(t) 
		return((lh :- mid :+ rh):/(h:^2))
	}
}

void merlin_mweights(`GMLS' gml)
{
	for (i=1;i<=1;i++) {
		gml.model = 2
		asarray(gml.weights,(1,1),1:/(*gml.Plnl[gml.model])(gml))
	}
}

`RM' merlin_identity(`RM' res)
{
	return(res)
}

`RM' merlin_exp(`RM' res)
{
	return(exp(res))
}

`RM' merlin_invlogit(`RM' res)
{
	return(invlogit(res))
}

`RM' merlin_tanh(`RM' res)
{
	return(tanh(res))
}

`RM' merlin_normal(`RM' res)
{
	return(normal(res))
}

`RM' merlin_invcll(`RM' res)
{
	return(1:-exp(-exp(res)))
}

`RM' merlin_rcs(real colvector variable,	///
                real rowvector knots,| 		///
                real scalar deriv,			///
                real matrix rmatrix			///
        )
{
	real scalar Nobs, Nknots, kmin, kmax, interior, Nparams
	real matrix splines, knots2
	
	if (args()==2) deriv = 0

	//======================================================================================================================================//
	// Extract knot locations

		Nobs 	= rows(variable)
		Nknots 	= cols(knots)
		kmin 	= knots[1,1]
		kmax 	= knots[1,Nknots]
	
		if (Nknots==2) interior = 0
		else interior = Nknots - 2
		Nparams = interior + 1
		
		splines = J(Nobs,Nparams,.)
		
	//======================================================================================================================================//
	// Calculate splines

		if (Nparams>1) {
			lambda = J(Nobs,1,(kmax:-knots[,2..Nparams]):/(kmax:-kmin))
			knots2 = J(Nobs,1,knots[,2..Nparams])
		}

		if (deriv==0) {
			splines[,1] = variable
			if (Nparams>1) {
				splines[,2..Nparams] = (variable:-knots2):^3 :* (variable:>knots2) :- lambda:*((variable:-kmin):^3):*(variable:>kmin) :- (1:-lambda):*((variable:-kmax):^3):*(variable:>kmax) 
			}
		}
		else if (deriv==1) {
			splines[,1] = J(Nobs,1,1)
			if (Nparams>1) {
				splines[,2..Nparams] = 3:*(variable:-knots2):^2 :* (variable:>knots2) :- lambda:*(3:*(variable:-kmin):^2):*(variable:>kmin) :- (1:-lambda):*(3:*(variable:-kmax):^2):*(variable:>kmax) 	
			}
		}
		else if (deriv==2) {
			splines[,1] = J(Nobs,1,0)
			if (Nparams>1) {
				splines[,2..Nparams] = 6:*(variable:-knots2) :* (variable:>knots2) :- lambda:*(6:*(variable:-kmin)):*(variable:>kmin) :- (1:-lambda):*(6:*(variable:-kmax)):*(variable:>kmax) 	
			}
		}
		else if (deriv==3) {
			splines[,1] = J(Nobs,1,0)
			if (Nparams>1) {
				splines[,2..Nparams] = 6:*(variable:>knots2) :- lambda:*6:*(variable:>kmin) :- (1:-lambda):*6:*(variable:>kmax)
			}
		}
		
		//orthog
		if (args()==4) {
			real matrix rmat
			rmat = luinv(rmatrix)
			if (deriv==0) splines = (splines,J(Nobs,1,1)) * rmat[,1..Nparams]
			else splines = splines * rmat[1..Nparams,1..Nparams]
		}
		return(splines)
}

`RM' merlin_orthog(`RM' x) 	//from rcsgen.ado SSC
{
	meanx = mean(x)
	v = x :- meanx ,J(rows(x),1,1) 
	q = J(rows(v),0,.)
	R = J(cols(v),cols(v),0)
	R[cols(v),] = (meanx,1)
	for (i=1;i<=cols(x);i++){
                r = norm(v[,i])/sqrt(rows(v))
                q = q, (v[,i]:/ r)
                R[i,i] = r
                for (j = i + 1; j<=cols(x); j++){
                        r = (q[,i]' * v[,j])/rows(v)
                        v[,j] = v[,j] - r*q[,i]
                        R[i,j] = r 
                }
	}
	return(R)
}



`RM' merlin_gq(real scalar n, string scalar inttype)
{
	i = range(1,n,1)'
	i1 = range(1,n-1,1)'
	alpha = 0
			
	if(inttype == "legendre") {
                muzero = 2
                a = J(1,n,0)
                b = i1:/sqrt(4 :* i1:^2 :- 1)
	}
	else if(inttype == "laguerre") {
                a = 2 :* i :- 1 :+ alpha
                b = sqrt(i1 :* (i1 :+ alpha))
                muzero = gamma(alpha :+ 1)
	}

	A= diag(a)
	for(j=1;j<=n-1;j++){
		A[j,j+1] = b[j]
		A[j+1,j] = b[j]
	}
	symeigensystem(A,vec,nodes)
	weights = (vec[1,]:^2:*muzero)'
	weights = weights[order(nodes',1)]
	nodes = nodes'[order(nodes',1)']
	
	return(nodes,weights)
}

`RM' merlin_drawnorm( `RC' mvec,		///
                        `RM' V,			///
                        `RS' n)	
{
	`RS' nvars
	`RM' cholV, z, res
	
	nvars = rows(mvec)
	if (nvars!=rows(V) & nvars!=cols(V)) {
		errprintf("Dimensions not valid\n")
		exit(198)
	}
	
	cholV = cholesky(V)
	
	z1 = runiform(nvars,n/2)
	z2 = z1,(1:-z1)
	z = invnormal(z2)
	//z = rnormal(nvars,n,0,1)

	//draw via x = mu + A*z
	res = J(nvars,n,.)
	for (i=1; i<=n; i++) res[.,i] = mvec + cholV*z[., i]
	return(res)
}

`RM' merlin_drawt( 				///
                        `RC' mvec,		///
                        `RM' V,			///
                        `RS' df,		///
                        `RS' n			///
                        )	
{
	`RM' res
	res = merlin_drawnorm(mvec,V,n):/sqrt(rchi2(rows(mvec),n,df):/df)
	return(res)
}

/*
Multivariate normal probability density function
*/

`RM' merlin_dmvnorm(`RM' X, `RM' V)
{
	`RM' invsigma, Vvecs, Vvals, q2, denom
	
	symeigensystem(V,Vvecs,Vvals)
	invsigma = Vvecs * (Vvecs' :/ Vvals')
	q2 = 0.5 :* rowsum((X*invsigma):*X)
	d1 = -0.5 :* (rows(V)*log(2*pi()) :+ sum(log(Vvals)) )
	return(exp(d1:-q2))
}

`RM' merlin_lndmvnorm(`RM' X, `RM' V)
{
	`RM' invsigma, Vvecs, Vvals, q2, denom
	
	symeigensystem(V,Vvecs,Vvals)
	invsigma = Vvecs * (Vvecs' :/ Vvals')
	q2 = 0.5 :* quadrowsum((X*invsigma):*X)
	d1 = -0.5 :* (rows(V)*log(2*pi()) :+ quadsum(log(Vvals)) )
	return(d1:-q2)
}

/*
Expand a matrix into all permutations
*/

`RM' merlin_expand_matrix(`RM' x,| `RS' weights )
{
	`RS' nrows, ncols, nexp, index
	`RM' newx
	
	if (rows(x)==1) return(x)
	else {
		nrows = rows(x)
		ncols = cols(x)	
		nexp  = ncols^(nrows-1)
		newx = J(1,nexp,x)
		index = 1
		for (i=1; i<=nrows-1; i++) {
			reps = ncols^(nrows-i)
			xrep = J(reps,1,x[i,])'
			xrep = rowshape(xrep,1)
			newx[i,] = J(1,index,xrep)
			index = index * ncols
		}
		if (args()>1) {
			for (i=2; i<=nrows; i++) newx[1,] = newx[1,] :* newx[i,]
			return(newx[1,])
		}
		else return(newx)
	}	
}

/*
Outer product by column
*/

`RM' merlin_outerprod_by_col(`RM' x)
{
	`RS' ncols, nrows
	`RM' newx
	nrows = rows(x)
	ncols = cols(x)
	newx = J(nrows,ncols^2,.)
	for (i=1;i<=nrows;i++) newx[i,] = x[i,] # x[i,] //rowshape(cross(x[i,],x[i,]),1)
	return(newx)
}

/*
error check for se(blups) when 0
*/

void merlin_seblup_check(`RS' x)
{
	if (x==0 | x==.) {
		errprintf("BLUP calculation failed in adaptive quadrature algorithm\n")
		errprintf("Try increasing gh()\n")
		exit(1986)
	}
}

`RM' merlin_fp(	`RC' variable,	///
				`RR' degs,|		///
				`RS' deriv		///
			)
{
	`RM' res
	if (args()==2 | deriv==0) {		
		if (cols(degs)==2) {
			if (degs[1,1]==degs[1,2]) {
				if (degs[1,1]==0) res = log(variable),log(variable):^2 
				else res = variable:^degs[1,1],log(variable):*variable:^degs[1,1]
			}
			else {
				res = J(rows(variable),2,.)
				if (degs[1,1]==0) res[,1] = log(variable)
				else res[,1] = variable:^degs[1,1]
				if (degs[1,2]==0) res[,2] = log(variable)
				else res[,2] = variable:^degs[1,2]
			}
		}
		else {
			if (degs[1,1]==0) res = log(variable) 
			else res = variable:^degs[1,1]
		}
	}
	else {
		if (cols(degs)==2) {
			if (degs[1,1]==degs[1,2]) {  //done
				if (degs[1,1]==0) res = (1:/variable),(2:*log(variable):/variable) 
				else res = (degs[1,1]:*variable:^(degs[1,1]:-1)),(variable:^(degs[1,1]:-1) :* (degs[1,1]:*log(variable) :+ 1))
			}
			else {
				res = J(rows(variable),2,.)							//done
				if (degs[1,1]==0) res[,1] = 1:/variable
				else res[,1] = degs[1,1]:*variable:^(degs[1,1]:-1)
				if (degs[1,2]==0) res[,2] = 1:/variable
				else res[,2] = degs[1,2]:*variable:^(degs[1,2]:-1)
			}
		}
		else {
			if (degs[1,1]==0) res = 1:/variable					//done
			else res = degs[1,1]:*variable:^(degs[1,1]:-1)		//done
		}
	}
	return(res)
}

void merlin_error(`SS' text)
{
	errprintf(text+"\n")
	exit(1986)
} 

`RM' head(`RM' x)
{
	return(x[1::10,])
}

`RM' merlin_bs(		`RC' x,			///
					`RR' knots, 	///
					`RS' degree,	///
					`RS' Nbasis,	///
					`RS' hasint		///
					)
{
	splines = J(rows(x),Nbasis,0)
	
	for (s=1; s<=Nbasis; s++) splines[,s] = merlin_bs_basis(x,knots,s,degree)
	
	//handle boundary knot 2
	index = selectindex(x:==knots[cols(knots)])
	Nri = rows(index)
	Nci = cols(index)
	if (Nri & Nri) splines[index,Nbasis] = J(Nri,1,1) 
	
	if (hasint) return(splines)
	else 		return(splines[,2..Nbasis])
	
}

`RC' merlin_bs_basis(`RC' x, `RR' knots, `RS' i, `RS' k)
{
	result 	= J(rows(x),1,0)

	if (k==0) {	
		index 	= selectindex((knots[i]:<=x) :* (x:<knots[i+1]))
		Ni		= rows(index)
		Ni2 	= cols(index)
		if (Ni & Ni2) result[index]	= J(Ni,1,1)
	} 
	else {	
		if ((knots[i+k]:-knots[i])!=0) 		result = result :+ merlin_bs_basis(x,knots,i,k-1) :* (x :- knots[i]) :/ (knots[i+k]:-knots[i])
		if ((knots[i+k+1]:-knots[i+1])!=0) 	result = result :+ merlin_bs_basis(x,knots,i+1,k-1) :* (knots[i+k+1]:-x) :/(knots[i+k+1]:-knots[i+1])
	}
	return(result)
	
}

`RM' merlin_ps(	`RC' x,		///
		`RR' knots, 	///
		`RS' degree,	///
		`RS' Nbasis,	///
		`RS' hasint	///
		)
{
	splines = J(rows(x),Nbasis,0)
	
	for (s=1; s<=Nbasis; s++) splines[,s] = merlin_bs_basis(x,knots,s,degree)
	
	//handle boundary knot 2
	index = selectindex(x:==knots[cols(knots)])
	Nri = rows(index)
	Nci = cols(index)
	if (Nri & Nri) splines[index,Nbasis] = J(Nri,1,1) 
	
	if (hasint) return(splines)
	else 		return(splines[,2..Nbasis])
	
}

`RM' merlin_xz_gp(`GMLS' gml, `RS' c, `RS' el, | `RC' t)
{
	hast		= args()==4
	gpinfo 		= asarray(gml.elinfo,(gml.model,c,el))
	
	istvar 		= asarray(rcsinfo,1)
	index		= merlin_get_index(gml)
	
	if (hast) {
		if (istvar) gpvar = t
		else 		gpvar = asarray(gpinfo,2)[index]
	}
	else {
		if (istvar) gpvar = asarray(gml.timevars,gml.model)[index]
		else 		gpvar = asarray(gpinfo,2)[index]
	}
	
	return(merlin_gp(gpvar))
}

`RM' merlin_pc(`RC' var, `RR' knots, `RS' noref)
{
	
	Nknots 	= cols(knots)
	Nobs 	= rows(var)
	
	if (noref) 	Nv = Nknots + 1
	else		Nv = Nknots
	
	splines = J(Nobs,Nv,.)

	if (noref) {
		//first interval
		splines[,1] = (var :< knots[1])
		//mid
		if (Nknots>1) {
			for (i=1;i<Nknots;i++) {
				splines[,i+1] = (knots[i]:<=var :& var:<knots[i+1])
			} 
		}
		//last interval
		splines[,Nknots+1] = (knots[Nknots]:<=var)
	}
	else {
		//mid
		if (Nknots>1) {
			for (i=1;i<Nknots;i++) {
				splines[,i] = (knots[i]:<=var :& var:<knots[i+1])
			} 
		}
		//last interval
		splines[,Nknots] = (knots[Nknots]:<=var)
	}
	
	return(splines)
}

`RM' designmatrix2(`RC' v, `RS' nc)
{
        scalar  i 
        matrix  res 

        if (rows(v)==0) return(J(0,0,.))
        res = J(rows(v),nc,0)
        for (i=1;i<=rows(res);i++) res[i,v[i]] = 1
        return(res)
}

`SS' merlin_get_letter(`RS' x)
{
	return(substr("abcdefghijklmnopqrstuvwxyz",x,1))
}

`RC' rows_cols(`RM' x)
{
	return(rows(x)\cols(x))
}

end
