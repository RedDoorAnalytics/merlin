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

// logl contributions for imputation models
`RS' merlin_logl_impute(`gml' gml)
{
	if (gml.familys[2]=="ibernoulli") {
	
		//for binary missing data model, need to evaluate for 0s and 1s

		//outcome model
		//-> for the outcome model x needs to be filled in with appropriate value being integrated over
		gml.ImputedValue = 0
		gml.model = gml.modtoind = 1
		res1 = (*gml.Plnl[1])(gml)
		//missing data model
		gml.model = gml.modtoind = 2
		res1 = res1 :+ merlin_logl_ibernoulli(gml)
		res1 = exp(res1)
		
		gml.ImputedValue = 1
		gml.model = gml.modtoind = 1
		res2 = (*gml.Plnl[1])(gml)
		//missing data model
		gml.model = gml.modtoind = 2
		res2 = res2 :+ merlin_logl_ibernoulli(gml)
		res2 = exp(res2)
		
		return(quadsum(log(res1 :+ res2),1))
	}
	else {
		
		//gaussian
		Ngh = gml.ImputeIP
		gh 	= _gauss_hermite_nodes(Ngh)
		x	= asarray(gml.distancb,(2,1)) :* gh[1,]' :* sqrt(2)
		w	= gh[2,]' :/ sqrt(pi())
		
		allres = 0
		xbi = merlin_util_xzb_mod(gml,2)
		for (q=1;q<=Ngh;q++) {
			gml.ImputedValue = xbi :+ x[q]	//add the mean
			gml.model = gml.modtoind = 1
			res1 = (*gml.Plnl[1])(gml)

			//missing data model
			gml.model = gml.modtoind = 2
// 			res1 = res1 :+ merlin_logl_igaussian(gml)
			res1 = exp(res1) :* w[q]
			allres = allres :+ res1
		}
		return(quadsum(log(allres),1))
	}
	
}

end
