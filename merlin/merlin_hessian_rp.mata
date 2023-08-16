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

`RM' merlin_rp_hessian_clp(`gml' gml)
{
	y	= merlin_util_depvar(gml)
	mod = gml.model
	nb 	= gml.eqnindex[mod,2] - gml.eqnindex[mod,1] + 1
	index1 = index2 = J(1,0,.)
	for (i=1; i<=nb; i++) {
		refind = 1
		while (refind<=i) {
			index1 = index1,i
			index2 = index2,refind
			refind++
		}
	}
	
	X = merlin_util_xz(gml,y[,1])
	Xprime = merlin_util_xz_deriv(gml,y[,1])
	
	brcs 	= asarray(gml.distancb,(mod,1))
	knots 	= asarray(gml.distancb,(mod,3))
	hasorth = asarray(gml.distancb,(mod,4))
	logt    = log(y[,1])
	
	if (hasorth) {
		rmat = asarray(gml.distancb,(mod,5))
		logh = - Xprime[,index1] :* Xprime[,index2] :/ ((merlin_rcs(logt,knots,1,rmat) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1])):^2
	}
	else {
		logh = - Xprime[,index1] :* Xprime[,index2] :/ ((merlin_rcs(logt,knots,1) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]) ):^2
	}
	
	return(y[,2] :* logh :-  merlin_rp_ch(gml,y[,1]) :* X[,index1] :* X[,index2])
	
}

`RM' merlin_rp_hessian_clp_dap(`gml' gml)
{
	y	= merlin_util_depvar(gml)
	mod = gml.model
	n1 	= gml.eqnindex[mod,2] - gml.eqnindex[mod,1] + 1
	n2 	= gml.Ndistancp[mod]
	
	index1 = index2 = J(1,0,.)
	for (i=1; i<=n2; i++) {
		for (j=1;j<=n1;j++) {
			index1 = index1,i
			index2 = index2,j
		}
	}
	
	brcs 	= asarray(gml.distancb,(mod,1))
	knots 	= asarray(gml.distancb,(mod,3))
	hasorth = asarray(gml.distancb,(mod,4))
	logt    = log(y[,1])
	
	X 		= merlin_util_xz(gml,y[,1])
	dX 		= merlin_util_xz_deriv(gml,y[,1])
	
	if (hasorth) {
		rmat = asarray(gml.distancb,(mod,5))
		SX	 = merlin_rcs(logt,knots,0,rmat)
		dSX	 = merlin_rcs(logt,knots,1,rmat)
		logh = - dSX[,index1] :* dX[,index2] :/ y[,1] :/ ((dSX * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1])):^2
	}
	else {
		SX = merlin_rcs(logt,knots,0)
		dSX = merlin_rcs(logt,knots,1)
		logh = - dSX[,index1] :* dX[,index2] :/ y[,1] :/ ((dSX * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]) ):^2
	}
	
	return(y[,2] :* logh :-  merlin_rp_ch(gml,y[,1]) :* X[,index2] :* SX[,index1])
}

`RM' merlin_rp_hessian_dap(`gml' gml)
{
	y	= merlin_util_depvar(gml)
	mod	= gml.model
	index1 = index2 = J(1,0,.)
	for (i=1; i<=gml.Ndistancp[mod]; i++) {
		refind = 1
		while (refind<=i) {
			index1 = index1,i
			index2 = index2,refind
			refind++
		}
	}
	
	brcs 	= asarray(gml.distancb,(mod,1))
	knots 	= asarray(gml.distancb,(mod,3))
	hasorth = asarray(gml.distancb,(mod,4))
	logt    = log(y[,1])
	
	if (hasorth) {
		rmat = asarray(gml.distancb,(mod,5))
		X = merlin_rcs(logt,knots,0,rmat)
		dX = merlin_rcs(logt,knots,1,rmat)
		logh = - (dX[,index1] :* dX[,index2] :/ y[,1]:^2) :/ ((merlin_rcs(logt,knots,1,rmat) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1])) :^ 2
	}
	else {
		X = merlin_rcs(logt,knots,0)
		dX = merlin_rcs(logt,knots,1)	
		logh = - (dX[,index1] :* dX[,index2] :/ y[,1]:^2) :/ ((merlin_rcs(logt,knots,1) * brcs) :/ y[,1] :+ merlin_util_xzb_deriv(gml,y[,1]) ) :^ 2
	}	
	return(y[,2] :* logh :- merlin_rp_ch(gml,y[,1]) :* X[,index1] :* X[,index2])
}


end
