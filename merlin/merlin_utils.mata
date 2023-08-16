*! version 1.0.0 ?????2017

/*
Utility functions
	
History
?????2017: version 0.1.0
*/

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
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

/*
	merlin_util_xzb()
	-> extract main complex predictor
*/

`RM' merlin_util_xzb(`gml' gml, | `RC' t, `RC' t0)
{
	`RS' hast, hast0, mod
	`RM' x, bindex, xb

	hast 	= args()>=2
	hast0 	= args()==3
	mod 	= gml.model

	if (hast) {
		if (hast0)	x = merlin_util_xz_simple(gml,t,t0)
		else 		x = merlin_util_xz_simple(gml,t)
	}
	else 			x = merlin_util_xz_simple(gml)

	bindex	= asarray(gml.X_bindex,(mod,1))
	xb 	= x[,bindex[1,]] * gml.myb[bindex[2,]]'

	if (gml.simple[mod]) return(xb)

	//now add in elements which contain [] i.e. random effects or EV etc. that will need updating dynamically
	bindex 	= asarray(gml.X_bindex,(mod,2))
	
	if (bindex!=J(2,0,.)) {

		Ncmps 	= gml.Ncmps[mod]		        //# of components
		Nels 	= asarray(gml.Nels,mod)			//# els per component
		cmpxix	= asarray(gml.CmpXBIndex,(mod,1))
		cmpbix	= asarray(gml.CmpXBIndex,(mod,2))

		if (hast) 	nobs = rows(t) 
		else 		nobs = merlin_get_nobs(gml)
	
		if (!gml.Nrelevels) {
			
			//if no random effects, the design matrix x can just be directly updated with any interactions
			//followed by matrix multiplication for finish the predictor
		
			for (c=1;c<=Ncmps;c++) {

				Xindex 		= cmpxix[c,1]..cmpxix[c,2]
				eltype	 	= asarray(gml.elindex,(mod,c))
				
				for (el=1; el<=Nels[c]; el++) {
					elmod = asarray(gml.elinfo,(mod,c,el))
					if (eltype[el]==4) {												//EV[]	
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_expval_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_expval_mod(gml,elmod)
						continue
					}
					else if (eltype[el]==5) {											//iEV[]
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_expval_integ_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_expval_integ_mod(gml,elmod)
						continue
					}
					else if (eltype[el]==6) {											//dEV[]
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_expval_deriv_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_expval_deriv_mod(gml,elmod)
						continue
					}
					else if (eltype[el]==7) {											//d2EV[]
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_expval_deriv2_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_expval_deriv2_mod(gml,elmod)
						continue
					}
					else if (eltype[el]==10) {											//XB[]		
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_mod(gml,elmod)
						continue
					}
					else if (eltype[el]==11) {											//iXB[]
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_integ_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_integ_mod(gml,elmod)
						continue
					}
					else if (eltype[el]==12) {											//dXB[]
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_deriv_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_deriv_mod(gml,elmod)
						continue
					}
					else if (eltype[el]==13) {											//d2XB[]
						if (hast) 	x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_deriv2_mod(gml,elmod,t)
						else 		x[,Xindex]  = x[,Xindex] :* merlin_util_xzb_deriv2_mod(gml,elmod)
					}	
				}
			}
			
			return(xb :+ x[,bindex[1,]] * gml.myb[,bindex[2,]]')
		}
		else {

			//random effects first, multiplying and adding each
			//then other elements, as they apply to the whole thing as they are a standard interaction
			//need to loop as anything could contain a random effect
			for (c=1;c<=Ncmps;c++) {
				
				Xindex 		= cmpxix[c,1]..cmpxix[c,2]
				Bindex 		= cmpbix[c,1]..cmpbix[c,2]
				NX			= cols(Xindex)
				eltype	 	= asarray(gml.elindex,(mod,c))
				contrib 	= 0
				
				//random effects must be handled first
				for (el=1; el<=Nels[c]; el++) {
					if (eltype[el]==2) {											//random effect
						xb2		= J(nobs,1,0)
						xb2base = x[,Xindex] :* gml.myb[,Bindex]
// 						x[,Xindex]
						for (re=1;re<=NX;re++) {
							if (gml.fixedonly==0) xb2 = xb2 :+ xb2base[,re] :* merlin_xz_b(gml,c,el,re)
							else if (gml.fixedonly==2) xb2 = xb2 :+ xb2base[,re] :* merlin_xz_blups(gml,c,el,re)		//-> for fitted predictions
							else if (gml.fixedonly==3) xb2 = xb2 :+ xb2base[,re] :* merlin_xz_specblups(gml,c,el,re)	//-> for fitted predictions
							//fixedonly are 0 already
							contrib = 1
						}
					}
				}

				if (!contrib) xb2 = x[,Xindex] * gml.myb[,Bindex]'

				//-> one for each colvector
				//-> other eltypes first (that are all colvectors)
				for (el=1; el<=Nels[c]; el++) {
					if (eltype[el]==4) {												//EV[]	
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_expval_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_expval_mod(gml,elmod)
						contrib = 1
						continue
					}
					else if (eltype[el]==5) {											//iEV[]
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_expval_integ_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_expval_integ_mod(gml,elmod)
						contrib = 1
						continue
					}
					else if (eltype[el]==6) {											//dEV[]
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_expval_deriv_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_expval_deriv_mod(gml,elmod)
						contrib = 1
						continue
					}
					else if (eltype[el]==7) {											//d2EV[]
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_expval_deriv2_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_expval_deriv2_mod(gml,elmod)
						contrib = 1
						continue
					}
					else if (eltype[el]==10) {											//XB[]		
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_xzb_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_xzb_mod(gml,elmod)
						contrib = 1
						continue
					}
					else if (eltype[el]==11) {											//iXB[]
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_xzb_integ_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_xzb_integ_mod(gml,elmod)
						contrib = 1
						continue
					}
					else if (eltype[el]==12) {											//dXB[]
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_xzb_deriv_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_xzb_deriv_mod(gml,elmod)
						contrib = 1
						continue
					}
					else if (eltype[el]==13) {											//d2XB[]
						elmod = asarray(gml.elinfo,(mod,c,el))
						if (hast) 	xb2 = xb2 :* merlin_util_xzb_deriv2_mod(gml,elmod,t)
						else 		xb2 = xb2 :* merlin_util_xzb_deriv2_mod(gml,elmod)
						contrib = 1
					}
					
				}
				if (contrib) xb = xb :+ xb2
			}
		}
	}
	return(xb)	
}

/*
	merlin_util_xzb_mod()
	-> extract complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_mod(`gml' gml, `RS' mod2, | `RC' t)
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
	merlin_util_depvar()
	-> extract dependent variable(s) for current model
	-> for survival outcomes, first column in time, second is event indicator
*/

`RM' merlin_util_depvar(`gml' gml)
{
	return(asarray(gml.y,gml.model)[merlin_get_index(gml),])
}

`RM' merlin_util_depvar_mod(`gml' gml, `RS' mod)
{
	return(asarray(gml.y,mod)[merlin_get_index(gml),])
}

/*
	merlin_util_ap()
	-> extract ancillary parameters for current model
*/

`RC' merlin_util_ap(`gml' gml, `RS' index)
{
	return(asarray(gml.apxb,(gml.model,index)))
}

`RC' merlin_util_ap_mod(`gml' gml, `RS' index, `RS' mod)
{
	return(asarray(gml.apxb,(mod,index)))
}

/*
	merlin_util_dap()
	-> extract distributional ancillary parameters for current model
*/

`RC' merlin_util_dap(`gml' gml, `RS' index)
{
	return(asarray(gml.distancb,(gml.model,index)))
}

`RC' merlin_util_dap_mod(`gml' gml, `RS' index, `RS' mod)
{
	return(asarray(gml.distancb,(mod,index)))
}

/*
	merlin_util_expval()
	-> extract expected value for current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval(`gml' gml, | `RC' t)
{
	hast = args()==2
	if (hast) 	result = (*gml.expvalP[gml.model])(gml,t)
	else 		result = (*gml.expvalP[gml.model])(gml)
	return(result)
}

/*
	merlin_util_expval_mod()
	-> extract expected value for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_mod(`gml' gml, `RS' mod2, | `RC' t)
{
	hast = args()==3
	if (!hast) {
		if (gml.issurv[gml.model] | gml.tvarnames[gml.model]!="") {
			t 		= merlin_util_timevar(gml)
			hast 	= 1
		}
	}
	
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

`RM' merlin_util_xzb_deriv(`gml' gml, | `RC' t)
{
	hast = args()==2
	if (gml.simple[gml.model]) {
		if (hast) 	return(merlin_util_xzb_deriv_simple(gml,t))
		else 		return(merlin_util_xzb_deriv_simple(gml))
	}
	
	if (hast) 	tvar = t
	else 		tvar = merlin_util_timevar(gml)

	smp  		= c("epsdouble")
	hstep 		= J(rows(tvar),1,1)
	index 		= selectindex(abs(tvar):<=1)
	hstep[index] 	= abs(tvar)[index]
	hstep 		= hstep :* smp :^(1/3)
	lh 		= merlin_util_xzb(gml,tvar :+ hstep)
	rh 		= merlin_util_xzb(gml,tvar :- hstep)
	return((lh :- rh):/(2:*hstep))
}

/*
	merlin_util_xzb_mod_deriv()
	-> extract d/dt of complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_deriv_mod(`gml' gml, `RS' mod2, | `RC' t)
{
	hast = args()==3
	if (!hast) {
		if (gml.issurv[gml.model] | gml.tvarnames[gml.model]!="") {
			t 		= merlin_util_timevar(gml)
			hast 	= 1
		}
	}

	tmpmod 		= gml.model
	gml.model 	= mod2
	if (args()==3) 	res = merlin_util_xzb_deriv(gml,t)
	else 			res = merlin_util_xzb_deriv(gml)
	gml.model 	= tmpmod
	return(res)
}

/*
	merlin_util_expval_deriv()
	-> extract d/dt of expected value of current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_deriv(`gml' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[merlin_get_index(gml)]
	
// gml.model,gml.modtoind,hast
// head(tvar)
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

`RM' merlin_util_expval_deriv_mod(`gml' gml, `RS' mod2, | `RC' t)
{
	hast = args()==3
	if (!hast) {
		if (gml.issurv[gml.model] | gml.tvarnames[gml.model]!="") {
			t 		= merlin_util_timevar(gml)
			hast 	= 1
		}
	}

	tmpmod 		= gml.model
	gml.model 	= mod2
	if (hast) 	res = merlin_util_expval_deriv(gml,t)
	else 		res = merlin_util_expval_deriv(gml)
	gml.model 	= tmpmod
	return(res)
}

/*
	merlin_util_xzb_integ()
	-> extract integral wrt t of complex linear predictor for a specific model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xzb_integ(`gml' gml, | `RC' t)
{
	ngq = 15
	gq 	= (0.991455371120813,-0.991455371120813,0.949107912342759,-0.949107912342759,0.864864423359769,-0.864864423359769,0.741531185599394,-0.741531185599394,0.586087235467691,-0.586087235467691,0.405845151377397,-0.405845151377397,0.207784955007898,-0.207784955007898,0)'
	gq 	= gq,(0.022935322010529,0.022935322010529,0.063092092629979,0.063092092629979,0.104790010322250,0.104790010322250,0.140653259715525,0.140653259715525,0.169004726639267,0.169004726639267,0.190350578064785,0.190350578064785,0.204432940075298,0.204432940075298,0.209482141084728)'
	//gq 		= merlin_gq(ngq,"legendre")
	if (args()==1) t = merlin_util_timevar(gml)
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

`RM' merlin_util_xzb_integ_mod(`gml' gml, `RS' mod2, | `RC' t)
{
	hast = args()==3
	if (!hast) {
		if (gml.issurv[gml.model] | gml.tvarnames[gml.model]!="") {
			t 		= merlin_util_timevar(gml)
			hast 	= 1
		}
	}
	
	tmpmod = gml.model
	gml.model = mod2
		if (hast) 	res = merlin_util_xzb_integ(gml,t)
		else 		res = merlin_util_xzb_integ(gml)
	gml.model = tmpmod
	return(res)
}

/*
	merlin_util_expval_integ()
	-> extract integral wrt t of expected value for current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_expval_integ(`gml' gml, | `RC' t)
{
	ngq 	= 15
	//gq 		= merlin_gq(ngq,"legendre")
	gq		= (0.991455371120813,-0.991455371120813,0.949107912342759,-0.949107912342759,0.864864423359769,-0.864864423359769,0.741531185599394,-0.741531185599394,0.586087235467691,-0.586087235467691,0.405845151377397,-0.405845151377397,0.207784955007898,-0.207784955007898,0)'
	gq 		= gq,( 0.022935322010529,0.022935322010529,0.063092092629979,0.063092092629979,0.104790010322250,0.104790010322250,0.140653259715525,0.140653259715525,0.169004726639267,0.169004726639267,0.190350578064785,0.190350578064785,0.204432940075298,0.204432940075298,0.209482141084728)'
	if (args()==1) t = merlin_util_timevar(gml)
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

`RM' merlin_util_expval_integ_mod(`gml' gml, `RS' mod2, | `RC' t)
{
	hast = args()==3
	if (!hast) {
		if (gml.issurv[gml.model] | gml.tvarnames[gml.model]!="") {
			t 		= merlin_util_timevar(gml)
			hast 	= 1
		}
	}

	tmpmod 		= gml.model
	gml.model 	= mod2
		if (hast) 	res = merlin_util_expval_integ(gml,t)
		else 		res = merlin_util_expval_integ(gml)
	gml.model 	= tmpmod
	return(res)
}


//2nd deriv

`RM' merlin_util_xzb_deriv2(`gml' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[merlin_get_index(gml)]
	
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

`RM' merlin_util_xzb_deriv2_mod(`gml' gml, `RS' mod2, | `RC' t)
{
	hast = args()==3
	if (!hast) {
		if (gml.issurv[gml.model] | gml.tvarnames[gml.model]!="") {
			t 		= merlin_util_timevar(gml)
			hast 	= 1
		}
	}

	tmpmod = gml.model
	gml.model = mod2
		if (hast) 	res = merlin_util_xzb_deriv2(gml,t)
		else 		res = merlin_util_xzb_deriv2(gml)
	gml.model = tmpmod
	return(res)
}

//2nd deriv

`RM' merlin_util_expval_deriv2(`gml' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[merlin_get_index(gml)]
	
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

`RM' merlin_util_expval_deriv2_mod(`gml' gml, `RS' mod2, | `RC' t)
{
	hast = args()==3
	if (!hast) {
		if (gml.issurv[gml.model] | gml.tvarnames[gml.model]!="") {
			t 		= merlin_util_timevar(gml)
			hast 	= 1
		}
	}
	
	tmpmod 		= gml.model
	gml.model 	= mod2
		if (hast) 	res = merlin_util_expval_deriv2(gml,t)
		else 		res = merlin_util_expval_deriv2(gml)
	gml.model 	= tmpmod
	return(res)
}

/*
	merlin_util_timevar()
	-> get timevar for specific model
*/

`RC' merlin_util_timevar(`gml' gml)
{
	return(asarray(gml.timevars,gml.model)[merlin_get_index(gml)])
}

`RC' merlin_util_timevar_mod(`gml' gml, `RS' mod)
{
	return(asarray(gml.timevars,mod)[merlin_get_index(gml)])
}

/*
	updates particular element values -> called from overoutcome() in predict
*/

void merlin_util_xzb_update(`gml' gml, `SR' vars, `RM' vals ,| `RS' mod)
{
	if (args()==3) mod = gml.model

	Nats 	= cols(vals)
	Ncmps 	= gml.Ncmps[mod]						//# of components
	Nels 	= asarray(gml.Nels,mod)					//# els per component

	nobs 	= merlin_get_nobs(gml)

	for (a=1; a<=Nats; a++) {
	
		for (c=1; c<=Ncmps; c++) {

			eltype	 = asarray(gml.elindex,(mod,c))

			for (el=1; el<=Nels[c]; el++) {

				if 		(eltype[el]==1)		merlin_xz_var_update(gml,c,el,vars[a],vals[,a])		//varname
// 				else if (eltype[el]==8)		merlin_xz_rcs_update(gml,c,el,vars[a],vals[a])		//rcs
// 				else if (eltype[el]==9) 	merlin_xz_fp_update(gml,c,el,vars[a],vals[a])		//fp
// 				else if (eltype[el]==14) 	merlin_xz_bs_update(gml,c,el,vars[a],vals[a])		//bs
// 				else if (eltype[el]==15) 	merlin_xz_ps_update(gml,c,el,vars[a],vals[a])		//ps
				
			}
			
		}
	}
	
}


/*
	merlin_util_xz()
	-> extract design matrix
*/

`RM' merlin_util_xz(`gml' gml, | `RC' t, `RC' t0)
{

	hast 	= args()==2
	hast0 	= args()==3
	
	if (gml.simple[gml.model]) {
		if (hast) {
			if (hast0)	return(merlin_util_xz_simple(gml,t,t0))
			else 		return(merlin_util_xz_simple(gml,t))
		}
		else 			return(merlin_util_xz_simple(gml))
	}
	
	
//!! below not called until gf1/2 for random effects	
//!! not worked on yet
	
	mod 	= gml.model
	Ncmps 	= gml.Ncmps[mod]						//# of components
	Nels 	= asarray(gml.Nels,mod)					//# els per component

	if (hast) 	nobs = rows(t) 
	else 		nobs = merlin_get_nobs(gml)
		
	X = J(nobs,gml.eqnindex[mod,2]-gml.eqnindex[mod,1]+1,.)

	xind = 1
	for (c=1;c<=Ncmps;c++) {

		eltype	 = asarray(gml.elindex,(mod,c))
		cmp		 = J(nobs,1,1)

		//build fixed component
		for (el=1; el<=Nels[c]; el++) {

			rebuild = 0
			
			if 		(eltype[el]==1)		{										//variable
				elvars 	= merlin_xz_var(gml,c,el)
				rebuild = 1
			}
			else if (eltype[el]==8) {											//rcs()
				if (hast0)		elvars = merlin_xz_rcs(gml,c,el,0,t,t0)
				else if (hast) 	elvars = merlin_xz_rcs(gml,c,el,0,t)
				else 			elvars = merlin_xz_rcs(gml,c,el,0)
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
				cmp = J(nobs,0,.)
				for (j=1;j<=Nnew;j++) {
					cmp = cmp,(copyold :* elvars[,j])	//!!poor
				}
			}

		}

		//apply coefficients here
// 		cmp 	= cmp * b[|1,cmpbix[c,1]\1,cmpbix[c,2]|]'
		
		//now add in elements which can contain random effects -> dimensions effects above
// 		nocontrib 	= 1
// 		recmp 		= J(nobs,1,1)
//
// 		for (el=1; el<=Nels[c]; el++) {
//
// 			if 		(eltype[el]==2) {											//random effect
// 				if (gml.fixedonly==0) 	recmp = recmp :* merlin_xz_b(gml,c,el)						
// 				else 					recmp = recmp :* 0						//-> for fixedonly predictions
// 				nocontrib 	= 0
// 				continue
// 			}
// 			else if (eltype[el]==4) {											//EV[]
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_expval_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_expval_mod(gml,elmod)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==5) {											//iEV[]
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_expval_integ_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_expval_integ_mod(gml,elmod)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==6) {											//dEV[]
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_expval_deriv_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_expval_deriv_mod(gml,elmod)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==7) {											//d2EV[]
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_expval_deriv2_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_expval_deriv2_mod(gml,elmod)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==3) {											//mf() 
// 				if (hast) 	recmp  = recmp :* merlin_xz_t(gml,c,el,t)
// 				else 		recmp  = recmp :* merlin_xz_t(gml,c,el)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==10) {											//XB[]		
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_xzb_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_xzb_mod(gml,elmod)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==11) {											//iXB[]
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_xzb_integ_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_xzb_integ_mod(gml,elmod)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==12) {											//dXB[]
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_xzb_deriv_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_xzb_deriv_mod(gml,elmod)
// 				nocontrib = 0
// 				continue
// 			}
// 			else if (eltype[el]==13) {											//d2XB[]
// 				elmod = asarray(gml.elinfo,(mod,c,el))
// 				if (hast) 	recmp  = recmp :* merlin_util_xzb_deriv2_mod(gml,elmod,t)
// 				else 		recmp  = recmp :* merlin_util_xzb_deriv2_mod(gml,elmod)
// 				nocontrib = 0
// 			}
//			
// 		}

// 		if (nocontrib) 	xb = xb :+ cmp
// 		else 			xb = xb :+ cmp :* recmp
		
		nnew = cols(cmp)
		if (nnew==1) 	X[,xind] = cmp
		else 			X[,xind..(xind+nnew-1)] = cmp
		
		xind = xind + nnew
	}

	if (gml.hascons[mod]) {
// 		xb = xb :+ b[cmpbix[c,1]] 		//intercept
		X[,xind] = J(nobs,1,1)
	}
	return(X)	
}

/*
	merlin_util_xz_deriv()
	
	-> extract d/dt of the design matrix of the complex linear predictor for the current model
	-> possibly dependent on new timevar
*/

`RM' merlin_util_xz_deriv(`gml' gml, | `RC' t)
{
	hast 	= args()==2
	
	if (gml.simple[gml.model]) {
		if (hast) {
			/*if (hast0)	return(merlin_util_xz_deriv_simple(gml,t,t0))
			else*/ 		return(merlin_util_xz_deriv_simple(gml,t))
		}
		else 			return(merlin_util_xz_deriv_simple(gml))
	}
	
	if (hast) 	tvar = t
	else 		tvar = asarray(gml.timevars,gml.model)[merlin_get_index(gml)]
	
	smp  		= c("epsdouble")
	hstep 		= J(rows(tvar),1,1)
	index 		= selectindex(abs(tvar):<=1)
	hstep[index] 	= abs(tvar)[index]
	hstep 		= hstep :* smp :^(1/3)
	lh 		= merlin_util_xz(gml,tvar :+ hstep)
	rh 		= merlin_util_xz(gml,tvar :- hstep)
	return((lh :- rh):/(2:*hstep))
}

`RR' merlin_util_istimevar(`gml' gml, `RS' c, `RS' el)
{
	return(asarray(gml.eltvar,(gml.model,c))[el,])
}

`RM' merlin_util_bhazard(`gml' gml)
{
	if 	(gml.hasbh[gml.model,1]) return(asarray(gml.bhazards,gml.model)[merlin_get_index(gml),])
	else if (gml.hasbh[gml.model,2]) return(asarray(gml.bhazards,(gml.model,gml.survind)))
}

`RM' merlin_util_bHazard(`gml' gml)
{
	return(asarray(gml.bhazards,(gml.model,gml.survind+2)))
}

end
