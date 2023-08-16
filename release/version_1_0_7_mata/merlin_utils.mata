*! version 1.0.0 ?????2017

/*
Utility functions for use in user-defined model
	
History
?????2017: version 0.1.0
*/

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

version 14.1

mata:

/*
	merlin_util_depvar()
	
	-> extract dependent variable(s) for current model
	-> for survival outcomes, first column in time, second is event indicator
*/

`RM' merlin_util_depvar(`GMLS' gml)
{
	return(asarray(gml.y,gml.model))
}

`RM' merlin_util_depvar_mod(`GMLS' gml, `RS' mod)
{
	return(asarray(gml.y,mod))
}

/*
	merlin_util_ap()
	
	-> extract ancillary parameters for current model
*/

`RC' merlin_util_ap(`GMLS' gml, `RS' index)
{
	return(asarray(gml.apxb,(gml.model,index)))
}

`RC' merlin_util_ap_mod(`GMLS' gml, `RS' index, `RS' mod)
{
	return(asarray(gml.apxb,(mod,index)))
}

/*
	merlin_util_xzb()
	
	-> extract main complex predictor
*/

`RM' merlin_util_xzb(`GMLS' gml, | `RC' t)
{

	hast 	= args()==2
	mod 	= gml.model
	moduse	= gml.modtoind

	Ncmps 	= gml.Ncmps[mod]						//# of components
	Nels 	= asarray(gml.Nels,mod)					//# els per component

	b 		= gml.myb
	cmpbix	= asarray(gml.CmpBIndex,mod)
	xb 		= J(gml.Nobs[gml.Nlevels,moduse],1,0)

	for (c=1;c<=Ncmps;c++) {

		eltype	 = asarray(gml.elindex,(mod,c))
		cmp		 = J(gml.Nobs[gml.Nlevels,moduse],1,1)

		//build fixed component
		for (el=1; el<=Nels[c]; el++) {

			rebuild = 0
			
			if 		(eltype[el]==1)		{										//variable
				elvars 	= merlin_xz_var(gml,c,el)
				rebuild = 1
			}
			else if (eltype[el]==8) {											//rcs()
				if (hast) 	elvars = merlin_xz_rcs(gml,c,el,0,t)
				else 		elvars = merlin_xz_rcs(gml,c,el,0)
				rebuild = 1
			}
			else if (eltype[el]==9) {											//fp()
				if (hast) 	elvars = merlin_xz_fp(gml,c,el,t)
				else 		elvars = merlin_xz_fp(gml,c,el)
				rebuild = 1
			}
			else if (eltype[el]==14) {											//bs()
				if (hast) 	elvars = merlin_xz_bs(gml,c,el,t)
				else 		elvars = merlin_xz_bs(gml,c,el)
				rebuild = 1
			}
			else if (eltype[el]==15) {											//ps()
				if (hast) 	elvars = merlin_xz_ps(gml,c,el,t)
				else 		elvars = merlin_xz_ps(gml,c,el)
				rebuild = 1
			}

			//rebuild for interactions
			if (rebuild) {
				Nold = cols(cmp)
				Nnew = cols(elvars)
				copyold = cmp
				cmp = J(gml.Nobs[gml.Nlevels,moduse],0,.)
				for (j=1;j<=Nnew;j++) {
					cmp = cmp,(copyold :* elvars[,j])	//!!poor
				}
			}
		}

		//apply coefficients here
		cmp 	= cmp * b[|1,cmpbix[c,1]\1,cmpbix[c,2]|]'
		
		//now add in elements which can contain random effects -> dimensions effects above
		nocontrib 	= 1
		recmp 		= J(gml.Nobs[gml.Nlevels,moduse],1,1)

		for (el=1; el<=Nels[c]; el++) {

			if 		(eltype[el]==2) {											//random effect
				if (gml.fixedonly==0) 	recmp = recmp :* merlin_xz_b(gml,c,el)						
				else 					recmp = recmp :* 0						//-> for fixedonly predictions
				nocontrib 	= 0
				continue
			}
			else if (eltype[el]==4) {											//EV[]
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_expval_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_expval_mod(gml,elmod)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==5) {											//iEV[]
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_expval_integ_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_expval_integ_mod(gml,elmod)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==6) {											//dEV[]
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_expval_deriv_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_expval_deriv_mod(gml,elmod)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==7) {											//d2EV[]
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_expval_deriv2_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_expval_deriv2_mod(gml,elmod)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==3) {											//mf() 
				if (hast) 	recmp  = recmp :* merlin_xz_t(gml,c,el,t)
				else 		recmp  = recmp :* merlin_xz_t(gml,c,el)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==10) {											//XB[]		
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_xzb_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_xzb_mod(gml,elmod)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==11) {											//iXB[]
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_xzb_integ_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_xzb_integ_mod(gml,elmod)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==12) {											//dXB[]
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_xzb_deriv_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_xzb_deriv_mod(gml,elmod)
				nocontrib = 0
				continue
			}
			else if (eltype[el]==13) {											//d2XB[]
				elmod = asarray(gml.elinfo,(mod,c,el))
				if (hast) 	recmp  = recmp :* merlin_util_xzb_deriv2_mod(gml,elmod,t)
				else 		recmp  = recmp :* merlin_util_xzb_deriv2_mod(gml,elmod)
				nocontrib = 0
			}
			
		}

		if (nocontrib) 	xb = xb :+ cmp
		else 			xb = xb :+ cmp :* recmp

	}

	if (gml.hascons[mod]) xb = xb :+ b[cmpbix[c,1]] 		//intercept
	return(xb)	
}

/*
	merlin_util_xzb_mod()
	
	-> extract complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	hast		= args()==3
	tmpmod 		= gml.model
	gml.model 	= mod2
		if (hast) 	result = merlin_util_xzb(gml,t)
		else 		result = merlin_util_xzb(gml)
	gml.model 	= tmpmod
	return(result)
}

/*
	merlin_util_expval()
	
	-> extract expected value for current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval(`GMLS' gml, | `RC' t)
{
	hast		= args()==2
	if (hast) 	result = (*gml.expvalP[gml.model])(gml,t)
	else 		result = (*gml.expvalP[gml.model])(gml)
	return(result)
}

/*
	merlin_util_expval_mod()
	
	-> extract expected value for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	hast		= args()==3
	tmpmod 		= gml.model
	gml.model 	= mod2
		if (hast) 	result = (*gml.expvalP[mod2])(gml,t)
		else 		result = (*gml.expvalP[mod2])(gml)
	gml.model 	= tmpmod
	return(result)
}


/*
	merlin_util_xzb_deriv()
	
	-> extract d/dt of complex linear predictor for the current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_deriv(`GMLS' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[asarray(gml.xbindex,gml.modtoind)]
	
	smp  			= c("epsdouble")
	hstep 			= J(rows(tvar),1,1)
	index 			= selectindex(abs(tvar):<=1)
	hstep[index] 	= abs(tvar)[index]
	hstep 			= hstep :* smp :^(1/3)
	lh 				= merlin_util_xzb(gml,tvar :+ hstep)
	rh 				= merlin_util_xzb(gml,tvar :- hstep)
	return((lh :- rh):/(2:*hstep))
}

/*
	merlin_util_xzb_mod_deriv()
	
	-> extract d/dt of complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_deriv_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	tmpmod = gml.model
	gml.model = mod2
		if (args()==3) 	res = merlin_util_xzb_deriv(gml,t)
		else 			res = merlin_util_xzb_deriv(gml)
	gml.model = tmpmod
	return(res)
}

/*
	merlin_util_expval_deriv()
	
	-> extract d/dt of expected value of current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_deriv(`GMLS' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[asarray(gml.xbindex,gml.modtoind)]
	
	smp  			= c("epsdouble")
	hstep 			= J(rows(tvar),1,1)
	index 			= selectindex(abs(tvar):<=1)
	hstep[index] 	= abs(tvar)[index]
	hstep 			= hstep :* smp :^(1/3)
	lh 				= merlin_util_expval(gml,tvar :+ hstep)
	rh 				= merlin_util_expval(gml,tvar :- hstep)
	return((lh :- rh):/(hstep:*2))
}

/*
	merlin_util_expval_mod_deriv()
	
	-> extract d/dt of complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_deriv_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	tmpmod = gml.model
	gml.model = mod2
		if (args()==3) 	res = merlin_util_expval_deriv(gml,t)
		else 			res = merlin_util_expval_deriv(gml)
	gml.model = tmpmod
	return(res)
}

/*
	merlin_util_xzb_integ()
	
	-> extract integral wrt t of complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_integ(`GMLS' gml, | `RC' t)
{
	ngq 	= 15
	//gq 		= merlin_gq(ngq,"legendre")
	gq	= (0.991455371120813,-0.991455371120813,0.949107912342759,-0.949107912342759,0.864864423359769,-0.864864423359769,0.741531185599394,-0.741531185599394,0.586087235467691,-0.586087235467691,0.405845151377397,-0.405845151377397,0.207784955007898,-0.207784955007898,0)'
	gq = gq,( 0.022935322010529,0.022935322010529,0.063092092629979,0.063092092629979,0.104790010322250,0.104790010322250,0.140653259715525,0.140653259715525,0.169004726639267,0.169004726639267,0.190350578064785,0.190350578064785,0.204432940075298,0.204432940075298,0.209482141084728)'
	nodes 	= t :* J(rows(t),1,gq[,1]'):/2 :+ t:/2
	weights = t :* J(rows(t),1,gq[,2]'):/2
	
	res = J(rows(t),1,0)
	for (q=1;q<=ngq;q++) {
		res = res :+ weights[,q] :* (*gml.invlinks[gml.model])(merlin_util_xzb(gml,nodes[,q]))	
	}
	return(res)
}

/*
	merlin_util_xzb_mod_integ()
	
	-> extract integral of complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_integ_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	tmpmod = gml.model
	gml.model = mod2
		if (args()==3) 	res = merlin_util_xzb_integ(gml,t)
		else 			res = merlin_util_xzb_integ(gml)
	gml.model = tmpmod
	return(res)
}

/*
	merlin_util_expval_integ()
	
	-> extract integral wrt t of expected value for current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_integ(`GMLS' gml, | `RC' t)
{
	ngq 	= 15
	//gq 		= merlin_gq(ngq,"legendre")
	gq		= (0.991455371120813,-0.991455371120813,0.949107912342759,-0.949107912342759,0.864864423359769,-0.864864423359769,0.741531185599394,-0.741531185599394,0.586087235467691,-0.586087235467691,0.405845151377397,-0.405845151377397,0.207784955007898,-0.207784955007898,0)'
	gq 		= gq,( 0.022935322010529,0.022935322010529,0.063092092629979,0.063092092629979,0.104790010322250,0.104790010322250,0.140653259715525,0.140653259715525,0.169004726639267,0.169004726639267,0.190350578064785,0.190350578064785,0.204432940075298,0.204432940075298,0.209482141084728)'
	nodes 	= t :* J(rows(t),1,gq[,1]'):/2 :+ t:/2
	weights = t :* J(rows(t),1,gq[,2]'):/2
	
	res = J(rows(t),1,0)
	for (q=1;q<=ngq;q++) {
		res = res :+ weights[,q] :* merlin_util_expval(gml,nodes[,q])	
	}
	return(res)
}

/*
	merlin_util_expval_mod_integ()
	
	-> extract integral wrt t of expected value for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_integ_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	tmpmod 		= gml.model
	gml.model 	= mod2
		if (args()==3) 	res = merlin_util_expval_integ(gml,t)
		else 			res = merlin_util_expval_integ(gml)
	gml.model 	= tmpmod
	return(res)
}


//2nd deriv

`RM' merlin_util_xzb_deriv2(`GMLS' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[asarray(gml.xbindex,gml.modtoind)]
	
	smp  			= c("epsdouble")
	hstep 			= J(rows(tvar),1,1)
	index 			= selectindex(abs(tvar):<=1)
	hstep[index] 	= abs(tvar)[index]
	hstep 			= hstep :* smp :^(1/3)
	lh 				= merlin_util_xzb(gml,tvar :+ hstep)
	rh 				= merlin_util_xzb(gml,tvar :- hstep)
	mid				= merlin_util_xzb(gml,tvar)
	return((lh :- 2 :* mid :+ rh):/(hstep:^2))
}

/*
	merlin_util_xzb_deriv2_mod()
	
	-> extract d2/dt2 of complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_deriv2_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	tmpmod = gml.model
	gml.model = mod2
		if (args()==3) 	res = merlin_util_xzb_deriv2(gml,t)
		else 			res = merlin_util_xzb_deriv2(gml)
	gml.model = tmpmod
	return(res)
}

//2nd deriv

`RM' merlin_util_expval_deriv2(`GMLS' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[asarray(gml.xbindex,gml.modtoind)]
	
	smp  			= c("epsdouble")
	hstep 			= J(rows(tvar),1,1)
	index 			= selectindex(abs(tvar):<=1)
	hstep[index] 	= abs(tvar)[index]
	hstep 			= hstep :* smp :^(1/3)
	lh 				= merlin_util_expval(gml,tvar :+ hstep)
	rh 				= merlin_util_expval(gml,tvar :- hstep)
	mid				= merlin_util_expval(gml,tvar)
	return((lh :- 2 :* mid :+ rh):/(hstep:^2))
}

/*
	merlin_util_expval_deriv2_mod()
	
	-> extract d2/dt2 of expected value for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_deriv2_mod(`GMLS' gml, `RS' mod2, | `RC' t)
{
	tmpmod = gml.model
	gml.model = mod2
		if (args()==3) 	res = merlin_util_expval_deriv2(gml,t)
		else 			res = merlin_util_expval_deriv2(gml)
	gml.model = tmpmod
	return(res)
}

/*
	merlin_util_timevar()
	
	-> get timevar for specific model
*/

`RC' merlin_util_timevar(`GMLS' gml)
{
	return(asarray(gml.timevars,gml.model))
}

`RC' merlin_util_timevar_mod(`GMLS' gml, `RS' mod)
{
	return(asarray(gml.timevars,mod))
}

end
