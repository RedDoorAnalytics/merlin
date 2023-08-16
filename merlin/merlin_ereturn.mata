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
	
	if 	(gml.E.merlin==4) title = "Survival model"
	else if (gml.E.merlin==3) title = "Artificial neural network"
	else if (gml.E.merlin==2) title = "Mixed effects survival model"
	else if (gml.E.merlin==5) title = "Fixed effects regression model"
	else if (gml.E.merlin==6) title = "Joint longitudinal-survival model"
        else if (gml.E.merlin==7) title = "Excess hazard model"
	else 			  title = "Mixed effects regression model"
	stata("ereturn local title "+title)
	
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
		stata("ereturn local bhazard"+strk+" "+gml.bhvarnames[k])
		stata("ereturn local latents"+strk+" ")
		if (gml.hasltrunc[k]) 	stata("ereturn local ltruncated"+strk+" "+gml.ltruncated[1,k])
		if (gml.haslint[k]) 	stata("ereturn local linterval"+strk+" "+gml.linterval[1,k])
		stata("ereturn local timevar"+strk+" "+gml.tvarnames[k])
		if (f=="pwexponential") {
			stata("ereturn local knots"+strk+" "+invtokens(strofreal(asarray(gml.distancb,(k,2)))))
		}
		if (f=="aft" | f=="rp" | f=="prp") {
			stata("ereturn local knots"+strk+" "+invtokens(strofreal(asarray(gml.distancb,(k,3)))))
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
		if (gml.renormal[i]) 	stata("ereturn local re_dist"+index+" normal")
		else 					stata("ereturn local re_dist"+index+" t")
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
	
        if (gml.indicatorvar!="") {
                stata("ereturn local indicator "+gml.indicatorvar)
        }
        
	if (gml.haspenalty) {
		if (gml.islasso) 	stata("ereturn local penalty Lasso")
		else 				stata("ereturn local penalty Ridge")
		stata("ereturn local lambda "+strofreal(gml.lambda))
	}
	
	stata("ereturn local chintpoints "+strofreal(gml.chip))
	
	if (gml.hastmat) {
		st_matrix("tmat",gml.tmat)
		stata("ereturn matrix transmatrix = tmat, copy")
	}
	
}

void merlin_ereturn_els(`gml' gml)
{
	mod		= gml.model
	Nels 	= asarray(gml.Nels,mod)
	
	if (gml.Ncmps[mod]) {
		for (i=1;i<=gml.Ncmps[mod];i++) {
			elindex = asarray(gml.elindex,(mod,i))
			for (j=1;j<=Nels[i];j++) {
				if (elindex[j]==8) {	//rcs()
					index 		= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(j)
					rcsinfo 	= asarray(gml.elinfo,(mod,i,j))
					knots		= asarray(rcsinfo,3)				
					stata("ereturn local knots_"+index+" "+invtokens(strofreal(knots,"%10.0g")))
					orth 		= asarray(rcsinfo,5)				
					if (orth) 	{
						stata("tempname rmat_"+index)
						st_matrix(st_local("rmat_"+index),asarray(rcsinfo,6))
						stata("ereturn matrix rmat_"+index+"="+st_local("rmat_"+index))
					}
				}
				else if (elindex[j]==17) {	//pc()
					index 		= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(j)
					info 		= asarray(gml.elinfo,(mod,i,j))
					knots		= asarray(info,3)
					stata("ereturn local knots_"+index+" "+invtokens(strofreal(knots,"%10.0g")))
				}
				else if (elindex[j]==14 | elindex[j]==15) {	//bs()/ps()
					index 			= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(j)
					rcsinfo 		= asarray(gml.elinfo,(mod,i,j))
					knots			= asarray(rcsinfo,3)			
					stata("ereturn local knots_"+index+" "+invtokens(strofreal(knots,asarray(rcsinfo,9))))
					degree			= asarray(rcsinfo,4)
					stata("ereturn local degree_"+index+" "+strofreal(degree))
					Nbasis			= asarray(rcsinfo,8)
					stata("ereturn local Nbasis_"+index+" "+strofreal(Nbasis))
				}
				
			}
		}
	}
}


end
