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

`RS' merlin_get_penalty(`gml' gml) 
{

	cmpbix	= asarray(gml.CmpXBIndex,(gml.model,2))
	bindex 	= 1..(max(cmpbix)-1)
	if (gml.islasso) 	return(gml.lambda * sum(abs(gml.myb[bindex])))		//lasso
	else 				return(gml.lambda * sum(gml.myb[bindex]:^2))		//ridge
}

`RR' merlin_get_deriv_penalty(`gml' gml) 
{
	cmpbix	= asarray(gml.CmpXBIndex,(gml.model,2))
	bindex 	= 1..(max(cmpbix)-1)	
	if (gml.islasso) {														//lasso
			deriv = gml.lambda :* gml.myb[bindex] :/ abs(gml.myb[bindex])
			_editmissing(deriv,0)
			return(deriv)							
	}
	else 	return(gml.lambda :* gml.myb[bindex]:*2)						//ridge
}

`RS' merlin_get_deriv2_penalty(`gml' gml, `RR' b) 
{
	if (gml.islasso) 	return(0)											//lasso
	else 				return(gml.lambda * 2)								//ridge
}



// `RS' merlin_penalise_ob(`gml' gml,`RR' b) 
// {
// 	p 	= 0
// 	pix = 1
	
// 	for (mod=1; mod<=gml.Nmodels; mod++) {
		
// 		Ncmps 	= gml.Ncmps[mod]						//# of components
// 		Nels 	= asarray(gml.Nels,mod)					//# els per component
// 		cmpbix	= asarray(gml.CmpBIndex,mod)
		
// 		for (c=1;c<=Ncmps;c++) {
// 			eltype	 = asarray(gml.elindex,(mod,c))
// 			if (Nels[c]==1 & eltype[1,1]==15) {
// 				avec = b[|1,cmpbix[c,1]\1,cmpbix[c,2]|] 

// 				p = p + gml.penals[pix++] :* merlin_get_penalty(avec) :/2
// 			} 	
// 		}	
		
// 	}
// 	return(p)
// }

// `RS' merlin_get_penalty(`RR' avec)
// {
// 	N = cols(avec)
// 	p = 0

// 	for (j=3;j<=N;j++) p = p + (avec[j] - 2 * avec[j-1] + avec[j-2])^2
// 	return(p)
// }

// void merlin_update_penalties(`TR' M, `gml' gml, `RR' b)
// {
	
// 	gml.H = merlin_get_H(M)
	
// 	pix = 1

// 	for (mod=1; mod<=gml.Nmodels; mod++) {

// 		Ncmps 	= gml.Ncmps[mod]						//# of components
// 		Nels 	= asarray(gml.Nels,mod)					//# els per component
// 		cmpbix	= asarray(gml.CmpBIndex,mod)

// 		for (c=1;c<=Ncmps;c++) {
// 			eltype	 = asarray(gml.elindex,(mod,c))
// 			if (Nels[c]==1 & eltype[1,1]==15) {
// 				avec = b[|1,cmpbix[c,1]\1,cmpbix[c,2]|]

// 				gml.penals[pix] :* merlin_get_penalty(avec) :/2
// 				Hll		= gml.H :+ gml.penals[pix] :* merlin_get_penalty(avec) :/2
				
// 				//optimize problem
// 				S = optimize_init()
// 				optimize_init_which(S,"min")
// 				optimize_init_argument(S, 1, gml) 
// 				optimize_init_argument(S, 2, avec) 
// 				optimize_init_argument(S, 3, Hll) 
// 				optimize_init_evaluator(S, &merlin_lcv())
// 				optimize_init_evaluatortype(S, "d0")
// 				optimize_init_params(S,log(gml.penals[pix]))	//initial value
// 				gml.penals[pix] = exp(optimize(S))				//estimated on log scale
// 				pix++
// 			} 	
// 		}	
		
// 	}
	
// }

// void merlin_lcv(todo, p, `gml' gml, `RR' avec, Hll, v, g ,H)
// {
	
// 	//LCV = -logl + pen/2 * trace(-H_pl,H_l)
	
// 	Hpl 	= Hll :- exp(p) :* merlin_get_penalty(avec) :/ 2
// 	v 		= -gml.lnfi1 :+ trace(invsym(Hpl),Hll)

// }

end


















