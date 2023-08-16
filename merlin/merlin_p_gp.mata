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

`RM' merlin_p_gp_mu(`gml' gml , | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	y		= merlin_util_depvar(gml)
	x 		= merlin_util_xzb(gml)
	
	l 		= exp(asarray(gml.distancb,(gml.model,1)))
	sd1 	= exp(asarray(gml.distancb,(gml.model,2)))
	sde		= exp(asarray(gml.distancb,(gml.model,3)))
	
	//K matrix
	nobs	= rows(y)
	K 		= sd1:^2 :* I(nobs)
	Kre		= sde:^2 :* I(nobs)
	
	for (i=1;i<=nobs;i++) {
		j=1
		while (j<i) {
			K[i,j] = K[j,i] = sd1:^2 :* exp(-((x[i] :- x[j]):^2):/(2:*l:^2))
			j++
		}
	}
	Kfull 	= K :+ Kre
	Kinv 	= invsym(Kfull)

	
}

end
