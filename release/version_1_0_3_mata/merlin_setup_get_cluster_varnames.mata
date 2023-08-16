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

void merlin_get_cluster_varnames(`gml' gml)
{

	//get unique cluster spec's within []
	gml.levelvars = J(2,0,"")													//first row contains level name
																				//second row contains M/MV#[] as well
	//get any level specs within [] - includes pipes
	idspec = J(0,1,"")
	for (i=1;i<=gml.Nmodels;i++) {
		vars 	= merlin_get_indepvars(i)
		Ndv 	= cols(vars)
		for (j=1;j<=Ndv;j++) {
			dv 		= strtrim(vars[j])
			pos1 	= strpos(dv,"[") ; pos2 = strpos(dv,"]") ;
			check 	= strpos(dv,"XB") |strpos(dv,"EV") | strpos(dv,"rcs(") | strpos(dv,"fp(")
			if (pos1 & !check) {
				idspec = idspec\substr(dv,pos1+1,pos2-pos1-1)
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

		/*
		//if family(re), need level specific one row indexes
		if (sum(gml.familys:=="re")) {
			alloneidvars = J(1,0,"")
			for (i=1;i<=Nlevs;i++) {
				varoneid = "_tempid"+strofreal(i)
				stata("tempvar "+varoneid)
				stata("qui bys "+invtokens(gml.levelvars[1,(1..i)])+": gen "+st_local(varoneid)+"=_n==1 if "+st_local("touse"))
				alloneidvars = alloneidvars,st_local(varoneid)
			}
			gml.oneidvars = st_data(.,alloneidvars,st_local("touse"))
		}*/
		
		//need to sort by id vars
		//stata("qui sort "+invtokens(gml.levelvars[1,])+",stable")
		
	}
	
}
end
