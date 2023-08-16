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

//no random effects in any model
`RC' merlin_logl_ob(`gml' gml)
{
	result = 0
	for (j=1;j<=gml.Nmodels;j++) {
		gml.model = gml.modtoind = j
		if (gml.NotNull[j,1]) {
			if (gml.hasweights[1]) result = result :+ quadsum((*gml.Plnl[j])(gml) :* asarray(gml.weights,(1,j)),1) 
			else result = result :+ quadsum((*gml.Plnl[j])(gml),1)
		}
	}
	return(result)
}

end
