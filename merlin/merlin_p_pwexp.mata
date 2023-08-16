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

`RM' merlin_p_pwexp_logh(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	logh  = merlin_util_xzb(gml,t)
	baseb = asarray(gml.distancb,(gml.model,1))
	cuts  = asarray(gml.distancb,(gml.model,2))
	Ncuts = cols(cuts)
	logh  = logh :+ merlin_pwexp_baselogh(gml,t,cuts,Ncuts,baseb)
	
	return(logh)
}

`RM' merlin_p_pwexp_h(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(exp(merlin_p_pwexp_logh(gml,t)))
}

`RM' merlin_p_pwexp_ch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	baseb = asarray(gml.distancb,(gml.model,1))
	cuts  = asarray(gml.distancb,(gml.model,2))
	Ncuts = cols(cuts)
	
	if (gml.NI[gml.model]) {
		nobs 	= rows(t)
		ch 		= J(nobs,1,0)
		Ngq 	= 30
		gq 		= merlin_gq(Ngq,"legendre")
		qp		= t :/ 2 :* J(nobs,1,gq[,1]') :+ t:/2
		loghazq = log(t:/2 :* J(nobs,1,gq[,2]'))
		for (q=1;q<=Ngq;q++) {
			ch = ch :+ exp(merlin_util_xzb(gml,qp[,q]) :+ merlin_pwexp_baselogh(gml,qp[,q],cuts,Ncuts,baseb) :+ loghazq[,q])
		}
		return(ch)
		
	}
	else return(exp(merlin_util_xzb(gml,t)) :* merlin_pwexp_basech(gml,t,cuts,Ncuts,baseb))
}

`RM' merlin_p_pwexp_logch(`gml' gml, | `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(log(merlin_p_pwexp_ch(gml,t)))
}

`RM' merlin_p_pwexp_s(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(exp(-merlin_p_pwexp_ch(gml,t)))
}

`RM' merlin_p_pwexp_f(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(1:-merlin_p_pwexp_s(gml,t))
}

`RM' merlin_p_pwexp_dens(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	return(merlin_p_pwexp_h(gml,t):*merlin_p_pwexp_s(gml,t))
}

end
