*! version 1.0.0 ?????2016

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
local Egml		struct merlin_ereturn_struct scalar
local RS 		real scalar
local SS 		string scalar
local PS 		pointer scalar
local RM 		real matrix
local SM		string matrix
local PC 		pointer colvector
local SC 		string colvector
local SR		string rowvector
local TR		transmorphic
local RC		real colvector
local RR		real rowvector
local PM		pointer matrix

version 14.2

mata:

void merlin_setup_error_checks(`gml' gml)
{
	for (i=1; i<=gml.Nmodels; i++) {
		
		if (st_local("timevar"+strofreal(i))!="") {
			if (gml.familys[i]=="gamma") {
				merlin_error("timevar() not currently supported with family(gamma)")
			}
			if (gml.familys[i]=="ggamma") {
				merlin_error("timevar() not currently supported with family(ggamma)")
			}
			if (gml.familys[i]=="lognormal") {
				merlin_error("timevar() not currently supported with family(lognormal)")
			}
			if (gml.familys[i]=="loglogistic") {
				merlin_error("timevar() not currently supported with family(loglogistic)")
			}
		}
	
	}
	anycox 		= sum(gml.familys:=="cox")
	anynotcox 	= sum(gml.familys:!="cox")
	if (anycox) {
		if (anynotcox)			merlin_error("family(cox) can't be combined with other families")
		if (gml.hasanylint) 	merlin_error("linterval() not supported with family(cox)")
		if (gml.Nmodels>1)		merlin_error("family(cox) only supported as a univariate model")
		if (gml.Nlevels>1) 		merlin_error("family(cox) doesn't currently allow random effects")
	}
	anyadd		= sum(gml.familys:=="addhazard")
	if (anyadd) {
		if (gml.Nlevels>1) 		merlin_error("family(addhazard) doesn't currently allow random effects")
	}
	
}

end
