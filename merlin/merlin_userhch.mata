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

`RM' merlin_userhaz_logh(`gml' gml, `RC' t)
{	
	mod	= gml.model
	haz = (*gml.Puserst[mod,1])(gml,t)
	if 		(gml.hasbh[mod,1]) 	logh = log(haz :+ merlin_util_bhazard(gml))
	else if (gml.hasbh[mod,2]) 	logh = log(haz :* merlin_util_bhazard(gml))
	else 						logh = log(haz)
	return(logh)
}

`RM' merlin_userhaz_ch(`gml' gml, `RC' t, | `RC' t0)
{	
	nobs 	= rows(t)
	ch 		= J(nobs,1,0)
	Ngq 	= gml.chip
	gq 		= merlin_gq(Ngq,"legendre")
	qp		= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
	qw		= t :/ 2 :* J(nobs,1,gq[,2]')
	if (args()==2) {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ (*gml.Puserst[gml.model,1])(gml,qp[,q]) :* qw[,q]
		}	
	}
	else {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ (*gml.Puserst[gml.model,1])(gml,qp[,q],t0) :* qw[,q]
		}	
	}
	return(ch)
}

`RM' merlin_userhaz_s(`gml' gml, `RC' t)
{	
	return(exp(-merlin_userhaz_ch(gml,t)))
}

`RM' merlin_userhaz_cdf(`gml' gml, `RC' t)
{	
	return(1:-merlin_userhaz_s(gml,t))
}

/*
-> functions for family(user, loghazard())
*/

`RM' merlin_userloghaz_logh(`gml' gml, `RC' t)
{
	mod	 = gml.model
	logh = (*gml.Puserst[mod,2])(gml,t)
	if 		(gml.hasbh[mod,1]) 	logh = log(exp(logh) :+ merlin_util_bhazard(gml))
	else if (gml.hasbh[mod,2]) 	logh = log(exp(logh) :* merlin_util_bhazard(gml))
	return(logh)
}

`RM' merlin_userloghaz_ch(`gml' gml, `RC' t, | `RC' t0)
{	
	nobs 	= rows(t)
	ch 		= J(nobs,1,0)
	Ngq 	= gml.chip
	gq 		= merlin_gq(Ngq,"legendre")
	qp		= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
	qw		= t :/ 2 :* J(nobs,1,gq[,2]')
	if (args()==2) {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ exp((*gml.Puserst[gml.model,2])(gml,qp[,q])) :* qw[,q]
		}	
	}
	else {
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ exp((*gml.Puserst[gml.model,2])(gml,qp[,q],t0)) :* qw[,q]
		}		
	}
	return(ch)
}

`RM' merlin_userloghaz_s(`gml' gml, `RC' t)
{	
	return(exp(-merlin_userloghaz_ch(gml,t)))
}

`RM' merlin_userloghaz_cdf(`gml' gml, `RC' t)
{	
	return(1:-merlin_userloghaz_s(gml,t))
}

/*
-> functions for family(user, hazard() chazard())
*/

`RM' merlin_userhazchaz_logh(`gml' gml, `RC' t)
{	
	mod	= gml.model
	haz = (*gml.Puserst[mod,1])(gml,t)
	if 		(gml.hasbh[mod,1]) 	logh = log(haz :+ merlin_util_bhazard(gml))
	else if (gml.hasbh[mod,2]) 	logh = log(haz :* merlin_util_bhazard(gml))
	else						logh = log(haz)
	return(logh)
}

`RM' merlin_userhazchaz_ch(`gml' gml, `RC' t, | `RC' t0)
{	
	if (args()==2) 	return((*gml.Puserst[gml.model,3])(gml,t))
	else 			return((*gml.Puserst[gml.model,3])(gml,t,t0))
}

`RM' merlin_userhazchaz_s(`gml' gml, `RC' t)
{	
	return(exp(-merlin_userhazchaz_ch(gml,t)))
}

`RM' merlin_userhazchaz_cdf(`gml' gml, `RC' t)
{	
	return(1:-merlin_userhazchaz_s(gml,t))
}

/*
-> functions for family(user, loghazard() chazard())
*/

`RM' merlin_userloghazchaz_logh(`gml' gml, `RC' t)
{	
	mod	 = gml.model
	logh = (*gml.Puserst[mod,2])(gml,t)
	if 		(gml.hasbh[mod,1]) 	logh = log(exp(logh) :+ merlin_util_bhazard(gml))
	else if (gml.hasbh[mod,2]) 	logh = log(exp(logh) :* merlin_util_bhazard(gml))
	return(logh)
}

`RM' merlin_userloghazchaz_ch(`gml' gml, `RC' t, | `RC' t0)
{	
	if (args()==2) 	return((*gml.Puserst[mod,3])(gml,t))
	else			return((*gml.Puserst[mod,3])(gml,t,t0))
}

`RM' merlin_userloghazchaz_s(`gml' gml, `RC' t)
{	
	return(exp(-merlin_userloghazchaz_ch(gml,t)))
}

`RM' merlin_userloghazchaz_cdf(`gml' gml, `RC' t)
{	
	return(1:-merlin_userloghazchaz_s(gml,t))
}

end
