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

`RM' merlin_logl_gp(`gml' gml)
{
	y		= merlin_util_depvar(gml)
	x 		= merlin_util_xzb(gml)
	len		= exp(asarray(gml.distancb,(gml.model,1)))
	sd1 	= exp(asarray(gml.distancb,(gml.model,2)))
	sde		= exp(asarray(gml.distancb,(gml.model,3)))
	Nobs	= gml.Nobs[gml.Nlevels,gml.model]
	K 		= sd1:^2 :* I(Nobs)
	
	for (i=1;i<=Nobs;i++) {
		j=1
		while (j<i) {
			K[i,j] = K[j,i] = sd1:^2 :* exp(-(abs(x[i] :- x[j]):^2):/(2:*len:^2))
			j++
		}
	}
	K 		= K :+ sde:^2 :* I(Nobs)
	return(lnmvnormalden(J(Nobs,1,0),K,y))
}

`RM' merlin_logl_gp_noresid(`gml' gml)
{
	y		= merlin_util_depvar(gml)
	x 		= merlin_util_xzb(gml)
	len		= exp(asarray(gml.distancb,(gml.model,1)))
	sd1 	= exp(asarray(gml.distancb,(gml.model,2)))
	Nobs	= gml.Nobs[gml.Nlevels,gml.model]
	K 		= sd1:^2 :* I(Nobs)
	
	for (i=1;i<=Nobs;i++) {
		j=1
		while (j<i) {
			K[i,j] = K[j,i] = sd1:^2 :* exp(-(abs(x[i] :- x[j]):^2):/(2:*len:^2))
			j++
		}
	}
	return(lnmvnormalden(J(Nobs,1,0),K,y))
}

`RM' merlin_gp_expval(`gml' gml , | `RC' t)
{
	has2 = args()==1
	if (has2) 	t = merlin_util_xzb(gml)
	
	y		= merlin_util_depvar(gml)
	x 		= merlin_util_xzb(gml)
	
	l 		= exp(asarray(gml.distancb,(gml.model,1)))
	sd1 	= exp(asarray(gml.distancb,(gml.model,2)))
	sde		= exp(asarray(gml.distancb,(gml.model,3)))
	
	//K matrix
	nobs	= rows(y)
	K 		= sd1:^2 :* I(nobs)
	
	for (i=1;i<=nobs;i++) {
		j=1
		while (j<i) {
			K[i,j] = K[j,i] = sd1:^2 :* exp(-(abs(x[i] :- x[j]):^2):/(2:*l:^2))
			j++
		}
	}
	
	//if out of sample -> add on to K matrix
	if (has2) {							
		nobs2 		= rows(t)
		Kstar 		= J(nobs2,nobs,.)
		Kstarstar 	= I(nobs2)
		
		for (i=1;i<=nobs2;i++) {
			for (j=1;j<=nobs;j++) {
				Kstar[i,j] = sd1:^2 :* exp(-(abs(t[i] :- x[j]):^2):/(2:*l:^2))
			}
		}	
		
		for (i=1;i<=nobs2;i++) {
			j=1
			while (j<i) {
				Kstarstar[i,j] = Kstarstar[i,j] = sd1:^2 :* exp(-(abs(t[i] :- t[j]):^2):/(2:*l:^2))
				j++
			}
		}	
		
		K = K , Kstar' \ Kstar , Kstarstar
	}
	else nobs2 = 0
	
	Kre		= sde:^2 :* I(nobs+nobs2)
	
	Kfull 	= K :+ Kre
	Kinv 	= invsym(Kfull)
	
	return(K * Kinv * y)
}

end
