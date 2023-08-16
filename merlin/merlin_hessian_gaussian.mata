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

`RM' merlin_gaussian_hessian_clp(`gml' gml)
{
	nb = asarray(gml.NHbs,gml.model)[1]
	index1 = index2 = J(1,0,.)
	for (i=1; i<=nb; i++) {
		refind = 1
		while (refind<=i) {
			index1 = index1,i
			index2 = index2,refind
			refind++
		}
	}

	y		= merlin_util_depvar(gml)
	X 		= merlin_util_xz(gml)
	expval 	= (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
	sigma	= merlin_util_dap(gml,1)
	
	return(-2 :* X[,index1] :* X[,index2] :/ (2 :* sigma:^2))
}

`RM' merlin_gaussian_hessian_clp_lsig(`gml' gml)
{
	y		= merlin_util_depvar(gml)
	expval 	= (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
	sigma	= merlin_util_dap(gml,1)

	return(-2 :* merlin_util_xz(gml) :* (y :- expval)  :/ (sigma:^2))
}


`RM' merlin_gaussian_hessian_lsigma(`gml' gml)
{
	y		= merlin_util_depvar(gml)
	expval 	= (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
	sigma	= merlin_util_dap(gml,1)

	return(-2 :* (y :- expval):^ 2 :/ (sigma :^ 2))
}

end
