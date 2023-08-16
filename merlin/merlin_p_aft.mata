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

`RM' merlin_p_aft_h(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)	
	smp  			= c("epsdouble")
	hstep 			= J(rows(t),1,1)
	index 			= selectindex(abs(t):<=1)
	hstep[index] 	= abs(t)[index]
	hstep 			= hstep :* smp :^(1/3)
	lh 				= merlin_p_aft_ch(gml,t :+ hstep)
	rh 				= merlin_p_aft_ch(gml,t :- hstep)
	return((lh :- rh):/(2:*hstep))
}

`RM' merlin_p_aft_logh(`gml' gml,| `RC' t)
{	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)	
	logh = log(merlin_p_aft_h(gml,t))
	return(logh)
}

`RM' merlin_p_aft_logch(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)	

	mod		= gml.model
	xb		= -merlin_util_xzb(gml,t)
	basis 	= log(t:*exp(xb))
	nc		= cols(basis)

	brcs		= asarray(gml.distancb,(mod,1))
	knots 		= asarray(gml.distancb,(mod,3))
	hasorthog 	= asarray(gml.distancb,(mod,4))
	if (hasorthog) rmat = asarray(gml.distancb,(mod,5))

	xbs = dxbs = J(merlin_get_nobs(gml),nc,.)
	if (hasorthog) { 
		for (i=1;i<=nc;i++) {
			xbs[,i]  = merlin_rcs(basis[,i],knots,0,rmat) * brcs
// 			dxbs[,i] = merlin_rcs(basis[,i],knots,1,rmat) * brcs
		}
	}
	else {
		for (i=1;i<=nc;i++) {
			xbs[,i]  = merlin_rcs(basis[,i],knots,0) * brcs
// 			dxbs[,i] = merlin_rcs(basis[,i],knots,1) * brcs
		}	
	}
	xbs = xbs :+ asarray(gml.distancb,(mod,2))
	return(xbs)
}

`RM' merlin_p_aft_ch(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)	
	return(exp(merlin_p_aft_logch(gml,t)))
}


`RM' merlin_p_aft_f(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)	

	return(1:-merlin_aft_s(gml,t))
}

`RM' merlin_p_aft_s(`gml' gml,| `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)	

	return(exp(-merlin_aft_ch(gml,t)))
}

`RM' merlin_p_aft_rmst(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	N = gml.Nobs[gml.Nlevels,gml.model]
	Nq = 30
	gq 	= merlin_gq(Nq,"legendre")
	qp	= t :/ 2 :* J(N,1,gq[,1]') :+ t:/2

	result = J(N,1,0)
	for (q=1; q<=Nq; q++) {
		result = result :+ merlin_aft_s(gml,qp[,q]) :* gq[q,2] :* t:/2
	}
	return(result)
}

`RM' merlin_p_aft_rmft(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(t:-merlin_p_aft_rmst(gml,t))
}

`RM' merlin_p_aft_dens(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	return(merlin_p_aft_h(gml,t):*merlin_p_aft_s(gml,t))
}

end


