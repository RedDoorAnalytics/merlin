
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

`RM' merlin_bernoulli_score_clp(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	expval 	= (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
	
	return(y :* log(expval) :+ (1:-y) :* log(1:-expval))
	
	return(merlin_util_xz(gml) :* (y :/ expval :- (1:-y) :/ (1:-expval)))
}

end
