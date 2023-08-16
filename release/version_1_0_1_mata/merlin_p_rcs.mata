*! version 1.0.0 ?????2016

local GMLS 		struct merlin_struct scalar
local pGMLS		pointer(struct merlin_struct scalar) scalar
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

`RM' merlin_p_rcs_h(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_p_rcs_logh(gml,t)))
}


`RM' merlin_p_rcs_logh(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	mod 		= gml.model
	brcs 		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	if (hasorthog) {
		rmat = asarray(gml.distancb,(mod,5))
		return(merlin_util_xzb(gml,t) :+ merlin_rcs(log(t),knots,0,rmat) * brcs)
	}
	else return(merlin_util_xzb(gml,t) :+ merlin_rcs(log(t),knots,0) * brcs)
}

`RM' merlin_p_rcs_ch(`GMLS' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
		
	gq 	= merlin_gq(15,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	result = J(gml.N,1,0)
	for (q=1; q<=15; q++) {
		result = result :+ merlin_p_rcs_h(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_rcs_s(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(exp(-merlin_p_rcs_ch(gml,t)))
}

`RM' merlin_p_rcs_f(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(1:-merlin_p_rcs_s(gml,t))
}

`RM' merlin_p_rcs_rmst(`GMLS' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	gq 	= merlin_gq(15,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2

	result = J(gml.N,1,0)
	for (q=1; q<=15; q++) {
		result = result :+ merlin_p_rcs_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_rcs_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_rcs_rmst(gml,t))
}

end
