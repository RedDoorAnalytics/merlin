*! version 1.0.0 MJC

local TS 	transmorphic scalar
local RS 	real scalar
local RC 	real colvector
local RM 	real matrix
local Pgml 	pointer(struct merlin_struct scalar) scalar
local gml 	struct merlin_struct scalar

version 14.2
mata:

`RM' merlin_neuralnet_score_ob(`gml' gml)
{
	G = J(gml.N,cols(gml.myb),0)

	ix = 1
	
	//output models
	
	for (j=1;j<=gml.Nmodels;j++) {
		gml.model = gml.modtoind = j
		if (gml.NotNull[j,1]) {
		
			xb = merlin_util_xzb(gml)
// 			expval 	= (*gml.invlinks[gml.model])(xb)	
// 			dL = y :* log(expval) :- y :* log(1:-expval) :+ log(y:*expval:^(-1):*(1:-expval) :- 1 :+ y) 
			dL = (1:+exp(-xb)):^(-2) :* exp(-xb) 

// 			dl		= exp(dL)
			//all elements are EV[] + intercept
		
			Ncmps 	= gml.Ncmps[j]						//# of components
			
			for (c=1;c<=Ncmps;c++) {
				
				G[,ix++] = asarray(gml.Score,(j,c)) :* dL
				
			}
		
			if (gml.hascons[j]) {
				G[,ix++] = dL
			}
		
		}
	}
	
	
	
	
	
	
	return(G)
	
}

end
