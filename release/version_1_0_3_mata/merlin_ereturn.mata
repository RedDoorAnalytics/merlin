*! version 1.0.0 ?????2016

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

version 14.1

mata:

void merlin_ereturn(`SS' GML)
{
	`gml' gml
	gml = *findexternal(GML)
	
	if (gml.predict) {
		stata("ereturn local struct struct")
	}
	stata("ereturn local allvars "+gml.allvars)
	
	//model info
	stata("ereturn scalar Nmodels = "+strofreal(gml.Nmodels))
	stata("ereturn local levelvars "+invtokens(gml.E.levelvars))
	
	ind = 1
	for (k=1;k<=gml.Nmodels;k++) {
		gml.model 	= k
		f 			= gml.familys[k]
		strk 		= strofreal(k)
		stata("ereturn local family"+strk+" "+gml.familys[k])
		if (f!="null") {
			stata("ereturn local response"+strk+" "+gml.responses[ind++])
		}
		stata("ereturn local failure"+strk+" "+gml.failures[k])
		stata("ereturn local latents"+strk+" ")
		if (gml.hasltrunc[k]) stata("ereturn local ltruncated"+strk+" "+gml.ltruncated[1,k])
		
		if (f=="rcs" | f=="rp") {
			stata("ereturn local rcs"+strk+" "+invtokens(strofreal(asarray(gml.distancb,(k,3)))))
			if (asarray(gml.distancb,(k,4))) {
				stata("ereturn local orthog"+strk+" orthog")
				stata("tempname rcsrmat_"+strk)
				st_matrix(st_local("rcsrmat_"+strk),asarray(gml.distancb,(k,5)))
				stata("ereturn matrix rcsrmat_"+strk+"="+st_local("rcsrmat_"+strk))
			}
		}
		
		if (f=="user") {
			uf = gml.E.userfunctions
			if (uf[k,1]!="") stata("ereturn local llfunction"+strk+" "+uf[k,1])
			if (uf[k,2]!="") stata("ereturn local hfunction"+strk+" "+uf[k,2])
			if (uf[k,3]!="") stata("ereturn local chfunction"+strk+" "+uf[k,3])
			if (uf[k,4]!="") stata("ereturn local loghfunction"+strk+" "+uf[k,4])
			
		}
		stata("ereturn local Nvars_"+strk+" "+gml.E.Nvars[k])
		stata("ereturn local cmplabels"+strk+" "+gml.cmplabels[k])
		merlin_ereturn_els(gml)			//spline knots and rmats for predictions
		stata("ereturn local constant"+strk+" "+strofreal(gml.hascons[k]))
		
		stata("ereturn local ndistap"+strk+" "+strofreal(gml.Ndistancp[k]))
		stata("ereturn local nap"+strk+" "+strofreal(gml.Nap[k]))
	}
	
	//number of levels, including ob level
	stata("ereturn scalar Nlevels = "+strofreal(gml.Nlevels))
	
	//random effect names, in sorted order (which matches eqns), at each level
	for (i=1;i<=gml.Nrelevels;i++) {
		index = strofreal(i)
		stata("ereturn local latents"+index+" "+invtokens(asarray(gml.latlevs,i)'))
		stata("ereturn local Nres"+index+" "+strofreal(gml.Nres[i]))
		stata("ereturn local Nreparams"+index+" "+strofreal(gml.E.Nreparams[i]))
		stata("ereturn local re_eqns"+index+" "+gml.E.reeqns[i])
		stata("ereturn local re_ivscale"+index+" "+gml.E.reivscale[i])
		stata("ereturn local re_label"+index+" "+gml.E.relabel[i])
	}
	
	//integration
	for (i=1;i<=gml.Nrelevels;i++) {
		index = strofreal(i)
		stata("ereturn local intpoints"+index+" "+invtokens(strofreal(gml.ip[i])))
		if (gml.usegh[i]) {
			if (gml.adapt[i]) 	stata("ereturn local intmethod"+index+" mvaghermite")
			else 				stata("ereturn local intmethod"+index+" ghermite")
		}
		else {
			if (gml.adapt[i]) 	stata("ereturn local intmethod"+index+" mvamcarlo")
			else 				stata("ereturn local intmethod"+index+" mcarlo")
		}	
	}
	
	
	
}

void merlin_ereturn_els(`gml' gml)
{
	mod		= gml.model
	Nels 	= asarray(gml.Nels,mod)
	
	if (gml.Ncmps[mod]) {
		for (i=1;i<=gml.Ncmps[mod];i++) {
			elindex = asarray(gml.elindex,(mod,i))
// 			stata("ereturn local elindex_m"+strofreal(mod)+"_c"+strofreal(i)+" "+invtokens(strofreal(elindex')))
			for (j=1;j<=Nels[i];j++) {
				//post knots and rmat for spline components
				if (elindex[j]==8) {
					index 		= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(j)
					rcsinfo 	= asarray(gml.elinfo,(mod,i,j))
					knots		= asarray(rcsinfo,3)				
					stata("ereturn local knots_"+index+" "+invtokens(strofreal(knots)))
					orth 		= asarray(rcsinfo,5)				
					if (orth) 	{
						stata("tempname rmat_"+index)
						st_matrix(st_local("rmat_"+index),asarray(rcsinfo,6))
						stata("ereturn matrix rmat_"+index+"="+st_local("rmat_"+index))
					}
				}
			}
		}
	}
}


end
