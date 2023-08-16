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

`RM' merlin_p_userh_h(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	return((*gml.userhaz[gml.model])(gml,t))
}

`RM' merlin_p_userh_ch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")

	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2
	result = J(N,1,0)

	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_userh_h(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_userh_logch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	return(log(merlin_p_userh_ch(gml,t)))
}

`RM' merlin_p_userh_s(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(-merlin_p_userh_ch(gml,t)))
}

`RM' merlin_p_userh_f(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(1:-merlin_p_userh_s(gml,t))
}

`RM' merlin_p_userh_rmst(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2

	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_userh_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_userh_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_userh_rmst(gml,t))
}

/*
-> functions for family(user, loghazard())
*/

`RM' merlin_p_userlogh_h(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	return(exp((*gml.userloghaz[gml.model])(gml,t)))
}

`RM' merlin_p_userlogh_ch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2
	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_userlogh_h(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_userlogh_logch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	return(merlin_p_userlogh_ch(gml,t))
}

`RM' merlin_p_userlogh_s(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(-merlin_p_userlogh_ch(gml,t)))
}

`RM' merlin_p_userlogh_f(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(1:-merlin_p_userlogh_s(gml,t))
}

`RM' merlin_p_userlogh_rmst(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2

	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_userlogh_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_userlogh_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_userlogh_rmst(gml,t))
}

/*
-> functions for family(user, hazard() chazard())
*/

`RM' merlin_p_userhch_h(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return((*gml.userhaz[gml.model])(gml,t))
}

`RM' merlin_p_userhch_ch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return((*gml.userchaz[gml.model])(gml,t))
}

`RM' merlin_p_userhch_logch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(log((*gml.userchaz[gml.model])(gml,t)))
}

`RM' merlin_p_userhch_s(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(-merlin_p_userhch_ch(gml,t)))
}

`RM' merlin_p_userhch_f(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(1:-merlin_p_userhch_s(gml,t))
}

`RM' merlin_p_userhch_rmst(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2
	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_userhch_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_userhch_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_userhch_rmst(gml,t))
}

/*
-> functions for family(user, loghazard() chazard())
*/

`RM' merlin_p_userloghch_h(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp((*gml.userloghaz[gml.model])(gml,t)))
}

`RM' merlin_p_userloghch_ch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return((*gml.userchaz[gml.model])(gml,t))
}

`RM' merlin_p_userloghch_logch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(log((*gml.userchaz[gml.model])(gml,t)))
}

`RM' merlin_p_userloghch_s(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(-merlin_p_userloghch_ch(gml,t)))
}

`RM' merlin_p_userloghch_f(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(1:-merlin_p_userloghch_s(gml,t))
}

`RM' merlin_p_userloghch_rmst(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2
	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_p_userloghch_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_userloghch_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_userloghch_rmst(gml,t))
}

end
