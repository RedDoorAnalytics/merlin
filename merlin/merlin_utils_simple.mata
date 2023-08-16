/*
Utility functions
*/

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
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

// merlin_util_xzb()	
// -> extract main complex predictor
`RM' merlin_util_xzb_simple(`gml' gml, | `RC' t, `RC' t0)
{
	hast 	= args()==2
	hast0 	= args()==3
	bindex	= asarray(gml.X_bindex,(gml.model,1))

	if 	(hast0) x = merlin_util_update_xz_simple(gml,t,t0)
	else if (hast) 	x = merlin_util_update_xz_simple(gml,t)
	else 		x = asarray(gml.X,gml.model)[merlin_get_index(gml),]

	return(x[,bindex[1,]] * gml.myb[bindex[2,]]')
}

// merlin_util_xz()	
// -> extract main design matrix
`RM' merlin_util_xz_simple(`gml' gml, | `RC' t, `RC' t0)
{
	if 	(args()==3)	return(merlin_util_update_xz_simple(gml,t,t0))
	else if (args()==2)     return(merlin_util_update_xz_simple(gml,t))
	else 		        return(asarray(gml.X,gml.model)[merlin_get_index(gml),])
}

// merlin_util_xz()	
// -> extract main design matrix
`RM' merlin_util_xz_col_simple(`gml' gml, `RS' index)
{
	return(asarray(gml.X,gml.model)[merlin_get_index(gml),index])
}

// merlin_util_xz()	
// -> extract main design matrix
`RM' merlin_util_xz_simple_mod(`gml' gml, `RS' mod, | `RC' t, `RC' t0)
{	
	tmpmod 		= gml.model
	gml.model 	= mod
	if 	(args()==3)	result = merlin_util_update_xz_simple(gml,t,t0)
	else if (args()==2)     result = merlin_util_update_xz_simple(gml,t)
	else 		        result = asarray(gml.X,gml.model)[merlin_get_index(gml),]
	gml.model 	= tmpmod
	return(result)
}

// merlin_util_xzb_deriv()	
// -> extract d/dt of complex linear predictor for the current model
`RM' merlin_util_xzb_deriv_simple(`gml' gml, | `RC' t, `RC' t0)
{
	hast 	= args()==2
	hast0 	= args()==3
	bindex	= asarray(gml.X_bindex,(gml.model,1))

	if 	(hast0)	x = merlin_util_update_dxz_simple(gml,t,t0)
	else if (hast) 	x = merlin_util_update_dxz_simple(gml,t)
	else 		x = asarray(gml.dX,gml.model)[merlin_get_index(gml),]
        
	return(x[,bindex[1,]] * gml.myb[bindex[2,]]')
}

`RM' merlin_util_xz_deriv_simple(`gml' gml, | `RC' t)
{
	if 	(args()==3)	return(merlin_util_update_dxz_simple(gml,t,t0))
	else if (args()==2)     return(merlin_util_update_dxz_simple(gml,t))
	else 		        return(asarray(gml.dX,gml.model)[merlin_get_index(gml),])
}

`RM' merlin_util_update_xz_simple(`gml' gml, `RC' t, | `RC' t0)
{
	`RM' x
	
	hast0	= args()==3
	mod 	= gml.model
	cmpxix	= asarray(gml.CmpXBIndex,(mod,1))
	Nels 	= asarray(gml.Nels,mod)
	Nobs	= merlin_get_nobs(gml)
	x 	= asarray(gml.X,mod)[merlin_get_index(gml),]	

	//update for any time-dependency
	for (c=1;c<=gml.Ncmps[mod];c++) {
	
		istvar 	 = asarray(gml.eltvar,(mod,c))
		anytvars = sum(istvar)

		if (anytvars) {

			eltype	= asarray(gml.elindex,(mod,c))
			cmp	= J(Nobs,1,1)

			//rebuild
			for (el=1; el<=Nels[c]; el++) {

				contrib = 0
				if 	(eltype[el]==1)	{					        //variable
					if 	(istvar[el,1]) 	        elvars = t
					else if (hast0 & istvar[el,2])	elvars = t0
					else 				elvars = merlin_xz_var(gml,c,el)
					contrib = 1
				}
				else if (eltype[el]==8) {							//rcs()
					if (hast0) 	elvars = merlin_xz_rcs(gml,c,el,0,t,t0)
					else		elvars = merlin_xz_rcs(gml,c,el,0,t)
					contrib = 1
				}
				else if (eltype[el]==9) {							//fp()
					if (hast0)	elvars = merlin_xz_fp(gml,c,el,0,t,t0)
					else 		elvars = merlin_xz_fp(gml,c,el,0,t)
					contrib = 1
				}
				else if (eltype[el]==14) {							//bs()
					if (hast0)	elvars = merlin_xz_bs(gml,c,el,t,t0)
					else 		elvars = merlin_xz_bs(gml,c,el,t)
					contrib = 1
				}
				else if (eltype[el]==17) {							//pc()
					if (hast0)	elvars = merlin_xz_pc(gml,c,el,t,t0)
					else 		elvars = merlin_xz_pc(gml,c,el,t)
					contrib = 1
				}
				else if (eltype[el]==3) {							//mf()
					elvars = merlin_xz_mf(gml,c,el,t)
					contrib = 1
				}

				//rebuild for interactions
				if (contrib) { 
					Nnew = cols(elvars)
					copyold = cmp
					cmp = J(Nobs,0,.)
					for (j=1;j<=Nnew;j++) {
						cmp = cmp,(copyold :* elvars[,j])	//!!poor
					}
				}
// 				else cmp = elvars
			}

			x[|.,cmpxix[c,1]\.,cmpxix[c,2]|] = cmp
		}
	}
	
	return(x)	
}


`RM' merlin_util_update_dxz_simple(`gml' gml, `RC' t)
{	
	mod 	= gml.model
	cmpxix	= asarray(gml.CmpXBIndex,(mod,1))
	Nels 	= asarray(gml.Nels,mod)
	Nobs	= merlin_get_nobs(gml)
	x 	= asarray(gml.dX,mod)[merlin_get_index(gml),]	
	Nels	= asarray(gml.Nels,mod)
	
	for (c=1;c<=gml.Ncmps[mod];c++) {
		
		eltype 		= asarray(gml.elindex,(mod,c))
		istvar 		= asarray(gml.eltvar,(mod,c))[,1]
		anytvars 	= sum(istvar)

		if (anytvars) { 
			
			//if more than one then need to use the chain rule
			if (anytvars>1) {
				
				newvars = J(gml.Nobs,1,0)
				for (el=1;el<=Nels[c];el++) {
				
					newvars2 = J(gml.Nobs,1,1)
					for (el2=1;el2<=Nels[c];el2++) {

						if (el!=el2) {
							if 	(eltype[el2]==1) nextvars = merlin_xz_var(gml,c,el2)
							else if (eltype[el2]==8) nextvars = merlin_xz_rcs(gml,c,el2,0)
							else if (eltype[el2]==9) nextvars = merlin_xz_fp(gml,c,el2,0)
						}
						else {
							if (istvar[el]) {
								if 	(eltype[el2]==1) nextvars = J(gml.Nobs,1,1)
								else if (eltype[el2]==8) nextvars = merlin_xz_rcs(gml,c,el2,1,t)
								else if (eltype[el2]==9) nextvars = merlin_xz_fp(gml,c,el2,1,t)
							}
							else {
								if 	(eltype[el]==1)	nextvars = J(gml.Nobs,1,0)
								else if (eltype[el]==8) {
									rcsinfo  = asarray(gml.elinfo,(gml.model,c,el))
									nextvars = J(gml.Nobs,cols(asarray(rcsinfo,3))-1,0)
								}
								else if (eltype[el]==9) {
									fpinfo 	 = asarray(gml.elinfo,(gml.model,c,el))
									nextvars = J(gml.Nobs,cols(asarray(fpinfo,3)),0)
								}
							}
						}
						
						//rebuild
						if (Nels[c]>1) {
							Nnew = cols(nextvars)
							copyold = newvars2
							newvars2 = J(gml.Nobs,0,.)
							for (j=1;j<=Nnew;j++) {
								newvars2 = newvars2,(copyold :* nextvars[,j])
							}	
						}
						else newvars2 = nextvars
					}
					newvars = newvars :+ newvars2
				}
					
			}
			else {
				newvars = J(gml.Nobs,1,1)
				for (el=1;el<=Nels[c];el++) {
					if (istvar[el]) {
						if 		(eltype[el]==1)		nextvars = J(gml.Nobs,1,1)
						else if (eltype[el]==8) 	nextvars = merlin_xz_rcs(gml,c,el,1,t)
						else if (eltype[el]==9) 	nextvars = merlin_xz_fp(gml,c,el,1,t)
					}
					else {
						if 		(eltype[el]==1)		nextvars = merlin_xz_var(gml,c,el) 
						else if (eltype[el]==8) 	nextvars = merlin_xz_rcs(gml,c,el,0)
						else if (eltype[el]==9) 	nextvars = merlin_xz_fp(gml,c,el,0)
					}
					//rebuild
					if (Nels[c]) {
						Nnew = cols(nextvars)
						copyold = newvars
						newvars = J(gml.Nobs,0,.)
						for (j=1;j<=Nnew;j++) {
							newvars = newvars,(copyold :* nextvars[,j])
						}
					}
					else newvars = nextvars
				}	
			}
			x[|.,cmpxix[c,1]\.,cmpxix[c,2]|] = newvars
		}
		
	}
	
	return(x)
}

end
