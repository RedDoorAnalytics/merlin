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

void merlin_get_levelvars(`gml' gml)
{
	//build new overall touse for id variables
	//which needs to cover unbalanced obs across models
	stata("tempvar idtouse")
	idtouse = st_local("idtouse")
	stata("qui gen byte "+idtouse+"= 0")
	for (i=1;i<=gml.Nmodels;i++) {
		stata("qui replace "+idtouse+"= 1 if "+st_local("touse"+strofreal(i))+"==1")
	}

	//build new idvars
	if (gml.Nlevels>1) {
		if (!gml.predict) stata("qui sort "+invtokens(gml.levelvars[1,])+",stable")
		//now create temps of them to account for gaps in id indexes
		newidvars = J(1,0,"")
		for (i=1;i<gml.Nlevels;i++) {
			tempid = "_tempidindex"+strofreal(i)
			stata("tempvar "+tempid)
			newidvars = newidvars,st_local(tempid)
			stata("qui egen "+st_local(tempid)+"= group("+invtokens(gml.levelvars[1,1..i])+") if "+idtouse)
		}
		gml.levelvars[1,] = newidvars
	}
	
	//ob level, now full touse is built
	stata("tempvar coreindex")
	coreindex = st_local("coreindex")
	stata("qui egen "+coreindex+" = seq() if "+idtouse)
	
	//need overall ids built with | model touses, which are then applied to the ids and stored
	gml.xbindex = asarray_create("real",1)
	for (i=1;i<=gml.Nmodels;i++) {
		asarray(gml.xbindex,i,st_data(.,coreindex,gml.modeltouses[1,i]))
		gml.Nobs[gml.Nlevels,i] = rows(asarray(gml.xbindex,i))	
	}

}

end
