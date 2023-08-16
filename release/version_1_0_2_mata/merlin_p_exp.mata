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

`RM' merlin_p_exp_h(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_util_xzb(gml,t)))
}

`RM' merlin_p_exp_ch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	if (gml.NI[gml.model]) {	
		gq 	= merlin_gq(15,"legendre")
		qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
		result = J(gml.N,1,0)
		for (q=1; q<=15; q++) {
			result = result :+ merlin_p_exp_h(gml,qp[,q]) :* gq[q,2] :* t:/2
		}
		return(result)
	}
	else return(exp(merlin_util_xzb(gml,t)) :* t )
}

`RM' merlin_p_exp_s(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(exp(-merlin_p_exp_ch(gml,t)))
}

`RM' merlin_p_exp_f(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(1:-merlin_p_exp_s(gml,t))
}

`RM' merlin_p_exp_rmst(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	gq 	= merlin_gq(15,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2

	result = J(gml.N,1,0)
	for (q=1; q<=15; q++) {
		result = result :+ merlin_p_exp_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_exp_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_exp_rmst(gml,t))
}

`RM' merlin_p_exp_mu(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(1:/exp(merlin_util_xzb(gml,t)))
}

end
