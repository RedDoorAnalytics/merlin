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

void merlin_setup_any_impute_models(`gml' gml)
{
	//model specific flag for imputation model
	gml.IsImputed 		= J(gml.Nmodels,1,0)
	gml.ImputeVarname 	= J(gml.Nmodels,1,"")
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.familys[i]=="ibernoulli" | gml.familys[i]=="igaussian") {
			gml.IsImputed[i] 		= 1
			gml.ImputeVarname[i] 	= st_local("response"+strofreal(i))
		}
	}
	gml.hasImputed = sum(gml.IsImputed)
}

//setup functions for any imputation models
void merlin_setup_impute(`gml' gml)
{
	if (gml.hasImputed) {
		gml.ImputeIndex = asarray_create("real",2)
		gml.ImputeIP = strtoreal(st_local("iintpoints"))
		merlin_impute_index(gml)
	}
}

void merlin_impute_index(`gml' gml)
{
	gml.ImputeIndex = asarray_create("real",1)
	for (i=1;i<=gml.Nmodels;i++) {
// 		if (gml.IsImputed[i]) {
			asarray(gml.ImputeIndex,i,st_data(.,st_local("coreindex"),st_local("imputetouse")))
// 		}
	}
}

end
