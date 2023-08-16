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

version 14.2

mata:

`RM' metamerlin_ad_gaussian(gml)
{
	y 	= merlin_util_depvar(gml)
	xb 	= merlin_util_xzb(gml)
	sd	= merlin_util_xzb_mod(gml,2)
	return(lnnormalden(y,xb,sd))
}

end
