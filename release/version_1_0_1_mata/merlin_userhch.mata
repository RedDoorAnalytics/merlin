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
-> functions for family(user, hazard())
*/

`RM' merlin_logl_userhaz(`gml' gml)
{	
	mod			= gml.model
	y 			= merlin_util_depvar(gml)

	if (gml.nobhaz[mod]) 	logh = y[,2] :* (log((*gml.userhaz[mod])(gml,y[,1])))
	else 					logh = y[,2] :*  log((*gml.userhaz[mod])(gml,y[,1]) :+ asarray(gml.bhazards,mod))
	
	return(logh :- merlin_userhaz_ch_ni(gml))
}

`RM' merlin_userhaz_logh(`gml' gml)
{	
	mod			= gml.model
	y 			= merlin_util_depvar(gml)

	if (gml.nobhaz[mod]) 	logh = log((*gml.userhaz[mod])(gml,y[,1]))
	else 					logh = log((*gml.userhaz[mod])(gml,y[,1]) :+ asarray(gml.bhazards,mod))
	return(logh)
}

`RM' merlin_userhaz_ch_ni(`gml' gml)
{	
	mod	= gml.model
	ch	= J(gml.Nobs[gml.Nlevels,mod],1,0)
	hq 	= asarray(gml.haznodes,mod)
	hw 	= asarray(gml.hazweights,mod)
	
	for (k=1;k<=gml.hazNnodes[mod];k++) {
		ch = ch :+ hw[,k] :* (*gml.userhaz[mod])(gml,hq[,k])
	}
	return(ch)
}

`RM' merlin_userhaz_s_ni(`gml' gml)
{	
	return(exp(-merlin_userhaz_ch_ni(gml)))
}

/*
-> functions for family(user, loghazard())
*/

`RM' merlin_logl_userloghaz(`gml' gml)
{
	mod			= gml.model
	y 			= merlin_util_depvar(gml)

	if (gml.nobhaz[mod]) 	logh = y[,2] :* ((*gml.userloghaz[mod])(gml,y[,1]))
	else 					logh = y[,2] :* log(exp((*gml.userloghaz[mod])(gml,y[,1])) :+ asarray(gml.bhazards,mod))
	
	return(logh :- merlin_userloghaz_ch_ni(gml))
}

`RM' merlin_userloghaz_logh(`gml' gml)
{
	mod			= gml.model
	y 			= merlin_util_depvar(gml)

	if (gml.nobhaz[mod]) 	logh = (*gml.userloghaz[mod])(gml,y[,1])
	else 					logh = log(exp((*gml.userloghaz[mod])(gml,y[,1])) :+ asarray(gml.bhazards,mod))
	
	return(logh)
}

`RM' merlin_userloghaz_ch_ni(`gml' gml)
{	
	mod	= gml.model
	ch	= J(gml.Nobs[gml.Nlevels,mod],1,0)
	hq 	= asarray(gml.haznodes,mod)
	hw 	= asarray(gml.hazweights,mod)
	
	for (k=1;k<=gml.hazNnodes[mod];k++) {
		ch = ch :+ hw[,k] :* exp((*gml.userloghaz[mod])(gml,hq[,k]))
	}
	return(ch)
}

`RM' merlin_userloghaz_s_ni(`gml' gml)
{	
	return(exp(-merlin_userloghaz_ch_ni(gml)))
}

/*
-> functions for family(user, hazard() chazard())
*/

`RM' merlin_logl_userhazchaz(`gml' gml)
{
	mod	= gml.model
	y 	= merlin_util_depvar(gml)
	
	if (gml.nobhaz[gml.model]) {
		return(y[,2] :* (log((*gml.userhaz[mod])(gml,y[,1]))) :- (*gml.userchaz[mod])(gml,y[,1]))
	}
	else {
		return(y[,2] :* log((*gml.userhaz[mod])(gml,y[,1]) :+ asarray(gml.bhazards,mod)) :- (*gml.userchaz[mod])(gml,y[,1]))
	}
}

`RM' merlin_userhazchaz_h(`gml' gml)
{	
	mod	= gml.model
	y 	= merlin_util_depvar(gml)
	return((*gml.userhaz[mod])(gml,y[,1]))
}

`RM' merlin_userhazchaz_ch(`gml' gml)
{	
	mod	= gml.model
	y 	= merlin_util_depvar(gml)
	return((*gml.userchaz[mod])(gml,y[,1]))
}

`RM' merlin_userhazchaz_s(`gml' gml)
{	
	return(exp(-merlin_userhazchaz_ch(gml)))
}

/*
-> functions for family(user, loghazard() chazard())
*/

`RM' merlin_logl_userloghazchaz(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	loghaz 	= (*gml.userloghaz[gml.model])(gml,y[,1]) 
	cumhaz 	= (*gml.userchaz[gml.model])(gml,y[,1])
	if (gml.nobhaz[gml.model]) return(y[,2] :* (loghaz) :- cumhaz)
	else return(y[,2] :* log(exp(loghaz) :+ asarray(gml.bhazards,gml.model)) :- cumhaz)
}

`RM' merlin_userloghazchaz_logh(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	loghaz 	= (*gml.userloghaz[gml.model])(gml,y[,1]) 
	if (gml.nobhaz[gml.model]) return(loghaz) 
	else return(log(exp(loghaz) :+ asarray(gml.bhazards,gml.model)))
}

`RM' merlin_userloghazchaz_ch(`gml' gml)
{
	y 	= merlin_util_depvar(gml)
	ch 	= (*gml.userchaz[gml.model])(gml,y[,1])
	return(ch)
}

`RM' merlin_userloghazchaz_s(`gml' gml)
{
	return(exp(-merlin_userloghazchaz_ch(gml)))
}






end
