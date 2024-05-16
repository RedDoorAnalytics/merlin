*! version 2.0.0  03mar2024

local gml 	struct merlin_struct scalar
local pgml	pointer(struct merlin_struct scalar) scalar
local TR 	transmorphic
local RS 	real scalar
local RC 	real colvector
local SS 	string scalar
local PS 	pointer scalar
local RR 	real rowvector
local RM 	real matrix
local PC 	pointer colvector
local PM 	pointer matrix
local SC 	string colvector

version 17

mata:

real matrix merlin_p_stexcess_h(`gml' gml, | `RC' t)
{
	if (args()==1) t = merlin_util_timevar(gml)
	haz_expect = exp(merlin_util_xzb(gml,t))
        haz_excess = exp(merlin_util_xzb_mod(gml,2,t))
	return(haz_expect :+ 
                        gml.indicator[merlin_get_index(gml)] :* 
                        haz_excess)
}

`RM' merlin_p_stexcess_ch(`gml' gml, | `RC' t)
{
	if (args()==1) t = merlin_util_timevar(gml)
	N 	= gml.Nobs[gml.Nlevels,gml.model]	
	ch 	= J(N,1,0)
	Ngq 	= gml.chip

	gq 	= merlin_gq(Ngq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2
	qw	= t :/ 2 :* J(N,1,gq[,2]')

	for (q=1;q<=Ngq;q++) {
		ch = ch :+ (exp(merlin_util_xzb(gml,qp[,q])) :+
			gml.indicator[merlin_get_index(gml)] :*
			exp(merlin_util_xzb_mod(gml,2,qp[,q]))) :* 
			qw[,q]
	}
	return(ch)
}

end
