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

void merlin_setup_levels(`gml' gml)
{
	merlin_get_cluster_varnames(gml)
	gml.E.levelvars = gml.levelvars[1,]											//to post
	gml.Nrelevels	= cols(gml.levelvars)
	gml.Nlevels 	= gml.Nrelevels + 1											//includes ob level
	gml.Nobs		= J(gml.Nlevels,gml.Nmodels,.)
}

void merlin_get_cluster_varnames(`gml' gml)
{
	//get unique cluster spec's within []
	gml.levelvars = J(2,0,"")		//first row contains level name
						//second row contains [] as well
	//get any level specs within [] - includes pipes
	idspec = J(0,1,"")
	for (i=1;i<=gml.Nmodels;i++) {
		vars 	= merlin_get_indepvars(i)
		Ndv 	= cols(vars)
		for (j=1;j<=Ndv;j++) {

			for (k=1;k<=Ndv;k++) {
				dv = strtrim(vars[1,k])
				pos = 1	
				while (pos) {
					pos = strpos(dv,"#")
					if (pos) {
						dv2 = substr(dv,1,pos-1)
						dv = substr(dv,pos+1,.)
					}
					else dv2 = dv
					
					pos1 	= strpos(dv2,"[") ; pos2 = strpos(dv2,"]") ;
					check 	= strpos(dv2,"XB") |strpos(dv2,"EV") | strpos(dv2,"rcs(") | strpos(dv2,"fp(")
					if (pos1 & !check) {
						idspec = idspec\substr(dv2,pos1+1,pos2-pos1-1)
					}
				
				}	
			}
						
		}
	}
	
	if (idspec!=J(0,1,"")) {

		//find unique ones to get # levels
		idspec = uniqrows(idspec)
		Nlevs  = rows(idspec)

		//get id vars, stripping out pipes
		idspecind = J(Nlevs,1,1)
		for (i=1;i<=Nlevs;i++) {
			ind 	= 1
			tempid 	= idspec[i]
			pos 	= strpos(tempid,">")
			while (pos) {
				idspecind[i] 	= idspecind[i] + 1
				tempid 			= substr(tempid,pos+1,.)
				pos 			= strpos(tempid,">")
			}
			pos 	= strpos(tempid,"<")
			while (pos) {
				idspecind[i] 	= idspecind[i] + 1
				tempid 			= substr(tempid,pos+1,.)
				pos 			= strpos(tempid,"<")
			}
			gml.levelvars 		= gml.levelvars,(tempid\("["+idspec[i]+"]"))
		}
		//make sure in order
		gml.levelvars = gml.levelvars[,idspecind']
		
	}
	
}
end
