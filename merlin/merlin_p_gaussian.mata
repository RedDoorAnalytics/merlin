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

`RM' merlin_p_gaussian_mu(`gml' gml , | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return((*gml.invlinks[gml.model])(merlin_util_xzb(gml,t)))
}

`RM' merlin_p_over_gaussian_mu(`gml' gml , `PS' pf, | `RC' t)
{
	not = args()==2
	if (not) t = merlin_util_timevar(gml)
	
	mod = gml.model
	overmod = strtoreal(st_local("overmodel"))
	overatvars = tokens(st_local("overatvars"))
	overatvals = strtoreal(tokens(st_local("overatvals")))
	
	Ngh = 15
	gh 	= _gauss_hermite_nodes(Ngh)
	x	= asarray(gml.distancb,(overmod,1)) :* gh[1,]' :* sqrt(2)
	w	= gh[2,]' :/ sqrt(pi())
		
	allres = 0

	//overoutcome
	gml.model = overmod
	`gml' gml2
	gml2 = gml
	merlin_util_xzb_update(gml2,overatvars,overatvals,overmod)
	xbi = merlin_util_xzb_mod(gml2,overmod)
	gml.model = mod

	for (q=1;q<=Ngh;q++) {
		
		//outcome
		merlin_util_xzb_update(gml,st_global(" e(response"+strofreal(overmod)+")"), xbi :+ x[q])	//add the mean
		res1 = (*pf)(gml,t)
		res1 = res1 :* w[q]
		allres = allres :+ res1
	
	}
	
	return(allres)
}

end
