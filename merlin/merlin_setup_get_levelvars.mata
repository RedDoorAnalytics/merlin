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

version 15.1

mata:

void merlin_get_levelvars(`gml' gml)
{
	//build new idvars
	if (gml.Nlevels>1) {
		if (!gml.predict) stata("qui sort "+invtokens(gml.levelvars[1,])+",stable")
		//now create temps of them to account for gaps in id indexes
		newidvars = J(1,0,"")
		for (i=1;i<gml.Nlevels;i++) {
			tempid = "_tempidindex"+strofreal(i)
			stata("tempvar "+tempid)
			newidvars = newidvars,st_local(tempid)
			stata("qui egen "+st_local(tempid)+
                                "= group("+invtokens(gml.levelvars[1,1..i])+
                                ") if "+gml.touse)
		}
		gml.levelvars[1,] = newidvars
	}
	
	//ob level, now full touse is built
	stata("tempvar coreindex")
	coreindex = st_local("coreindex")
	stata("qui egen "+coreindex+" = seq() if "+gml.touse)
	
	//need overall ids built with | model touses, which are then applied to the ids and stored
	gml.xbindex = asarray_create("real",1)
	for (i=1;i<=gml.Nmodels;i++) {
		asarray(gml.xbindex,i,
                        st_data(.,coreindex,gml.modeltouses[1,i]))
		gml.Nobs[gml.Nlevels,i] = rows(asarray(gml.xbindex,i))	
	}

	gml.levelvars = gml.levelvars,(st_local("coreindex")\st_local("coreindex"))
	
	//id vars setup
	gml.panelindexes = asarray_create("real",2)
	if (gml.Nlevels>1) {
		for (j=1;j<=gml.Nmodels;j++) {
			ids = st_data(.,gml.levelvars[1,],gml.modeltouses[1,j])
			for (i=1;i<gml.Nlevels;i++) {
				psetup 		= panelsetup(uniqrows(ids[,1..(i+1)]),i)
				gml.Nobs[i,j]	= rows(psetup)
				asarray(gml.panelindexes,(i,j),psetup)
			}
		}
		
		//error check for at least one observation per outcome, per 2nd lowest level
		if (gml.Nmodels>1) {
			for (j=1;j<gml.Nmodels;j++) {
				if (gml.Nobs[gml.Nrelevels,j] != gml.Nobs[gml.Nrelevels,j+1]) {
			errprintf("There must be at least one observation per outcome, per cluster\n")
			exit(198)
				}	
			}
		}

		//left truncation panel indices
// 		if (gml.hasanyltrunc) {
// 			gml.lt_panelindexes = asarray_create("real",2)
// 			for (j=1;j<=gml.Nmodels;j++) {
// 				ids = st_data(.,levelvars[1,],gml.modeltouses[1,j])
// 				ltids = select(ids,_t0:>0)
// 				asarray(gml.surv_index,(i,4))
				
				
// 				for (i=1;i<gml.Nlevels;i++) {
// 					psetup 				= panelsetup(uniqrows(ids[,1..(i+1)]),i)
// 					gml.Nobs[i,j]	 	= rows(psetup)
// 					asarray(gml.panelindexes,(i,j),psetup)
// 				}
// 			}			
// 		}
		
	}


}

end
