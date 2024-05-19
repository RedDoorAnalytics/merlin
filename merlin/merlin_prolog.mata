*! version 1.0.0 MJC

local RS real scalar
local RC real colvector
local RR real rowvector
local RM real matrix
local pGML pointer(struct merlin_struct scalar) scalar
local GML struct merlin_struct scalar
local TR transmorphic
local TS transmorphic scalar

version 14.2
mata:

// prolog
void merlin_prolog(	`RR' b,		///
                        `TS' M, 	///
                        `RS' lnl) 
{
	`pGML' p
	p = &moptimize_util_userinfo(M,1)

	if ((*p).gridsearch) merlin_vcv_gridsearch(*p,b)

	if (sum((*p).adapt) & (*p).iter<1001 & (*p).iter<=(*p).adaptit) {

		`GML' 	gml
		`RS' 	oldlnl, newlnl , it
		
		gml = moptimize_util_userinfo(M,1)
		
		// need to fill up linear predictors and myb (used in _utils), 
		// first time through
		merlin_xb(gml,b)
		gml.myb = b
		oldlnl 	= quadsum(merlin_logl_panels(1,gml),1)
		
		it = 0	
		if (gml.showadapt) merlin_di_adaptiter(it,oldlnl)
		
		it = 1
		(*gml.Pupdateip)(gml)
		newlnl = quadsum(merlin_logl_panels(1,gml),1)
		if (gml.showadapt) merlin_di_adaptiter(it,newlnl)
		
		while (reldif(oldlnl,newlnl)>gml.atol & it<=1000) {
			it++
			swap(oldlnl,newlnl)
			(*gml.Pupdateip)(gml)
			newlnl = quadsum(merlin_logl_panels(1,gml),1)
			if (gml.showadapt) merlin_di_adaptiter(it,newlnl)
		}
		
		//update stored nodes -> needs work!
		merlin_gh_post_ip(gml)

		//===========================================================//
		//left truncation

		if (gml.hasmargltrunc) {	

			gml.ltflag = 1
			oldlnl = quadsum(merlin_ltrunc(gml),1)
			
			it = 0	
			if (gml.showadapt) merlin_di_adaptiter(it,oldlnl)
			
			it = 1
			merlin_ltrunc_gh_update_ip(gml)
			newlnl = quadsum(merlin_ltrunc(gml),1)
			if (gml.showadapt) merlin_di_adaptiter(it,newlnl)
			
			while (reldif(oldlnl,newlnl)>gml.atol & it<=1000) {
				it++
				swap(oldlnl,newlnl)
				merlin_ltrunc_gh_update_ip(gml)
				newlnl = quadsum(merlin_ltrunc(gml),1)		
				if (gml.showadapt) {
					merlin_di_adaptiter(it,newlnl)
				}
			}
			//update stored nodes
			merlin_ltrunc_gh_post_ip(gml)
		
		}
	}

// 	if ((*p).Npenals>0 & (*p).iter>0) {
// 		gml = moptimize_util_userinfo(M,1)
// 		merlin_xb(gml,b)
// 		merlin_update_penalties(M,gml,b)
// 	}
	
	p->iter = (*p).iter :+ 1
	
}

/*
	display adapted log-likelihood
*/
void merlin_di_adaptiter(`RS' iter, `RS' val)
{
	st_numscalar("mliter",iter)
	st_numscalar("ll",val)		
	stata(`"di as txt "-- Iteration " _col(14) mliter as txt ":" _col(19) "Adapted log likelihood = " as res %10.0g ll"')
}

void merlin_mc_update(	`GML' gml, 	///
                        `RS' lev, 	///
                        `RR' b, 	///
                        `TS' M, 	///
                        `RS' lnl) 
{
	//increase in steps of ?
	gml.Pgml->ndim[lev] = gml.ndim[lev] = gml.ndim[lev] + 2			//must update external and passed struct
	merlin_update_Zs_bs(gml,lev)	//!!re-write to add 3 new draws to existing
}

// Update adaptive quadrature locations and scales
void merlin_gh_update_ip(`GML' gml)
{
	for (i=1; i<=1; i++) {
		ndim = gml.ndim[i]
		ind1 = 1; ind2 = ndim
		stackednodes 	= J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
		detvcv 		= J(gml.Nobs[i,1],1,.)
		basenodes 	= asarray(gml.baseGHnodes,i)
		baseweights 	= asarray(gml.baseGHweights,i)
		L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
		numer 		= asarray(gml.Li_ip,gml.qind) :/ L_i
		baseweights2 	= baseweights'
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			ipij 	 = asarray(gml.aghip,(i,j))
			newblups = numer[j,] * (ipij :* baseweights2)'
			vcv_new  = (numer[j,] * 
					(merlin_outerprod_by_col(ipij') :* 
					baseweights)) :- 
					newblups # newblups
			vcv_new  = rowshape(vcv_new,gml.Nres[i])
			nodes 	 = newblups' :+ cholesky(vcv_new) * basenodes
			asarray(gml.aghip,(i,j),nodes)
			detvcv[j] = det(vcv_new)
			stackednodes[|ind1,.\ind2,.|] = nodes'
			ind1 = ind1 + ndim
			ind2 = ind2 + ndim
		}
		//update logl extra contribution
		asarray(gml.aghlogl,
			i,
			rowshape(merlin_lndmvnorm(stackednodes,I(gml.Nres[i])),
				 gml.Nobs[i,1]) 
			:+ log(detvcv):/2)
		//update stacked nodes, and RE specific stacked nodes
		asarray(gml.stackednodes,i,stackednodes')
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2,
				(i,r),
				rowshape(res[r,],gml.Nobs[i,1]))
		}	
	}
}

void merlin_gh_post_ip(`GML' gml)
{
	for (i=1; i<=1; i++) {
		asarray(gml.Pgml->stackednodes,i,asarray(gml.stackednodes,i))
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			asarray(gml.Pgml->aghip,(i,j),asarray(gml.aghip,(i,j)))
		}
		for (r=1; r<=gml.Nres[i]; r++) {
			asarray(gml.Pgml->aghip2,(i,r),asarray(gml.aghip2,(i,r)))
		}
		asarray(gml.Pgml->aghlogl,i,asarray(gml.aghlogl,i))
	}
}	

// Update adaptive quadrature locations and scales
void merlin_ltrunc_gh_update_ip(`GML' gml)
{
	for (i=1; i<=1; i++) {
		ndim = gml.ndim[i]
		Nobs = max(gml.Nsurv[,7])
		ind1 = 1
		ind2 = ndim
		stackednodes 	= J(ndim*Nobs,gml.Nres[i],.)
		detvcv 		= J(Nobs,1,.)
		basenodes 	= asarray(gml.baseGHnodes,i)
		baseweights 	= asarray(gml.baseGHweights,i)
		L_i 		= asarray(gml.Li_ip_lt,gml.qind) * baseweights
		numer 		= asarray(gml.Li_ip_lt,gml.qind) :/ L_i
		baseweights2 	= baseweights'

		for (j=1; j<=Nobs; j++) {
			ipij 		= asarray(gml.aghip_lt,(i,j))
			newblups 	=  numer[j,] * (ipij :* baseweights2)'
			vcv_new 	= (numer[j,] * (merlin_outerprod_by_col(ipij') :* baseweights)) :- newblups # newblups
			vcv_new 	= rowshape(vcv_new,gml.Nres[i])
			nodes 		= newblups' :+ cholesky(vcv_new) * basenodes
			asarray(gml.aghip_lt,(i,j),nodes)
			detvcv[j] 	= det(vcv_new)
			stackednodes[|ind1,.\ind2,.|] = nodes'
			ind1 = ind1 + ndim
			ind2 = ind2 + ndim
		}

		//update logl extra contribution
		asarray(gml.aghlogl_lt,i,rowshape(merlin_lndmvnorm(stackednodes,I(gml.Nres[i])),Nobs) :+ log(sqrt(detvcv)))
		//update stacked nodes, and RE specific stacked nodes
		asarray(gml.stackednodes_lt,i,stackednodes')
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes_lt,i)
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2_lt,(i,r),rowshape(res[r,],Nobs))
		}
	}
}

void merlin_ltrunc_gh_post_ip(`GML' gml)
{
	for (i=1; i<=1; i++) {
		Nobs = max(gml.Nsurv[,7])
		asarray(gml.Pgml->stackednodes_lt,i,asarray(gml.stackednodes_lt,i))
		for (j=1; j<=Nobs; j++) 	asarray(gml.Pgml->aghip_lt,(i,j),asarray(gml.aghip_lt,(i,j)))
		for (r=1;r<=gml.Nres[i];r++) 	asarray(gml.Pgml->aghip2_lt,(i,r),asarray(gml.aghip2_lt,(i,r)))
		asarray(gml.Pgml->aghlogl_lt,i,asarray(gml.aghlogl_lt,i))
	}
}

void merlin_mc_update_ip(`GML' gml) 
{

	ip = gml.ip
	ndim = ip

	for (i=1; i<=1; i++) {
		
		ind1 = 1; ind2 = ndim
		stackednodes = J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
		stackednodes2 = J(gml.Nobs[i,1],ndim,.)
		L_i = quadrowsum(asarray(gml.Li_ip,gml.qind))
		numer = asarray(gml.Li_ip,gml.qind) :/ L_i
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			ipij 		= asarray(gml.aghip,(i,j))
			newblups 	=  numer[j,] * (ipij')
			vcv_new 	= (numer[j,] * (merlin_outerprod_by_col(ipij'))) :- rowshape((newblups')*newblups,1)
			vcv_new 	= rowshape(vcv_new,gml.Nres[i])
			
			/*rseed(gml.seed)
			nodes = merlin_drawnorm(newblups',vcv_new,ip)
			asarray(gml.aghip,(i,j),nodes)
			nodes2 = cholesky(asarray(gml.vcvs,i)) * nodes
			asarray(gml.aghip2,(i,j),nodes2)*/

			nodes = newblups' :+ cholesky(vcv_new) * asarray(gml.bdraws,i)
			asarray(gml.aghip,(i,j),nodes)
			nodes2 = cholesky(asarray(gml.vcvs,i)) * nodes
			asarray(gml.aghip2,(i,j),nodes2)
			for (k=1;k<=ndim;k++) stackednodes2[j,k] = lnmvnormalden(newblups', vcv_new, nodes[,k])
			stackednodes[|ind1,.\ind2,.|] = nodes'
			ind1 = ind1 + ndim
			ind2 = ind2 + ndim
		}
		//update logl extra contribution
		//asarray(gml.aghlogl,i,rowshape(merlin_lndmvnorm(stackednodes,I(gml.Nres[i])),gml.Nobs[i,1]) :+ logdetChol)
		
		//update stacked nodes, and RE specific stacked nodes
		asarray(gml.stackednodes,i,stackednodes')
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2,(i,r),rowshape(res[r,],gml.Nobs[i,1]))
		}
		
		asarray(gml.aghlogl,i,rowshape(merlin_lndmvnorm(stackednodes,I(gml.Nres[i])),gml.Nobs[i,1]) :- stackednodes2)

	}	
}

void merlin_update_ip_newNip(`GML' gml, `RS' i)
{
	//update nodes
	if (gml.adapt[i]) {
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)	
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2,(i,r),rowshape(res[r,],gml.Nobs[i,1]))
		}	
	}
	else {
		if (gml.usegh[i]) {
// 			if (gml.hasRExb[i]) {
// 				vcv = asarray(gml.vcvs,i)
// 				ndim = gml.ndim[i]
// 				ind1 = 1; ind2 = ndim
// 				stackednodes = J(gml.ndim[i]*gml.Nobs[i,1],gml.Nres[i],.)
// 				for (j=1;j<=gml.Nobs[i,1];j++) {
// 					nodes = cholesky(invvech(vcv[j,])) * asarray(gml.baseGHnodes,i)
// 					stackednodes[|ind1,.\ind2,.|] = nodes'
// 					ind1 = ind1 + ndim
// 					ind2 = ind2 + ndim
// 				}
// 				stackednodes = stackednodes'			//!!fix
// 				for (r=1;r<=gml.Nres[i];r++) {
// 					asarray(gml.aghip2,(i,r),rowshape(stackednodes[r,],gml.Nobs[i,1]))
// 				}				
// 			}
// 			else {
				bnodes = cholesky(asarray(gml.vcvs,i)) * asarray(gml.baseGHnodes,i)
				asarray(gml.b,i,bnodes)	
// 			}
		}
		else {
			if (gml.renormal[i]) 	{
				asarray(gml.b,i,cholesky(asarray(gml.vcvs,i)) * asarray(gml.bdraws,i))
			}
			else {	
				rseed(gml.seed)			//reset seed
				bnodes = merlin_drawt(J(gml.Nres[i],1,0),asarray(gml.vcvs,i),gml.df[i],gml.ndim[i])
				asarray(gml.b,i,bnodes)	
			}
		}
			
	}	
}

// Update adaptive quadrature locations and scales
void merlin_gh_update_ip_alllevs(`GML' gml)
{
	
	gml.Nrelevels
	gml.adapt
	gml.Nobs
	for (i=1; i<=gml.Nrelevels; i++) {
		
		ndim = gml.ndim[i]
		ind1 = 1; ind2 = ndim
		11
		stackednodes 	= J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
		detvcv 		= J(gml.Nobs[i,1],1,.)
		basenodes 	= asarray(gml.baseGHnodes,i)
		baseweights 	= asarray(gml.baseGHweights,i)
		33
		L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
		numer 		= asarray(gml.Li_ip,gml.qind) :/ L_i
		baseweights2 	= baseweights'
		44
		numer
		for (j=1; j<=gml.Nobs[i,1]; j++) {
			j
			ipij 	 = asarray(gml.aghip,(i,j))
			asarray(gml.aghip,(i,j))
			newblups = numer[j,] * (ipij :* baseweights2)'
			j,j
			vcv_new  = (numer[j,] * 
					(merlin_outerprod_by_col(ipij') :* 
					baseweights)) :- 
					newblups # newblups
			vcv_new  = rowshape(vcv_new,gml.Nres[i])
			nodes 	 = newblups' :+ cholesky(vcv_new) * basenodes
			asarray(gml.aghip,(i,j),nodes)
			detvcv[j] = det(vcv_new)
			stackednodes[|ind1,.\ind2,.|] = nodes'
			ind1 = ind1 + ndim
			ind2 = ind2 + ndim
		}
		99
		//update logl extra contribution
		asarray(gml.aghlogl,
			i,
			rowshape(merlin_lndmvnorm(stackednodes,I(gml.Nres[i])),
				 gml.Nobs[i,1]) 
			:+ log(detvcv):/2)
		56
		//update stacked nodes, and RE specific stacked nodes
		asarray(gml.stackednodes,i,stackednodes')
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2,
				(i,r),
				rowshape(res[r,],gml.Nobs[i,1]))
		}	
	}
}

end
