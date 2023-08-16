*! version 1.0.0 ?????2016

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

void merlin_starting_values(`gml' gml, `TR' lats)
{

	//starting values
	binit = J(1,gml.Nb,0)

	//distap defaults -> needed for rp first spline term and ordinal cuts
	if (st_local("from")=="") binit[gml.initdistapindex[1,]] = gml.initdistapindex[2,]
	
	if (gml.Nrelevels & st_local("from")=="" & st_local("zeros")=="") {
		
		cmd = st_local("ZERO")

		//get rid of all components which have a [ in them and put back together again
		
		stata("_parse expand EQ GL : ZERO")
		newcmd = ""
		for (mod=1; mod<=gml.Nmodels; mod++) {
			clp = st_local("EQ_"+strofreal(mod))
			stata("local 0 "+clp)
			stata("syntax anything , [*]")
			opts = st_local("options")
			clp = st_local("anything")
			newclp = J(1,0,"")
			st_local("clp",clp)
			stata(`"gettoken cmp clp : clp, parse(" ") bind"')
			if (!strpos(st_local("cmp"),"[")) newclp = newclp,st_local("cmp")
			while (st_local("clp")!="") {
				stata(`"gettoken cmp clp : clp, parse(" ") bind"')
				if (!strpos(st_local("cmp"),"[")) newclp = newclp,st_local("cmp")
			}
			newcmd = newcmd +" ("+invtokens(newclp)+","+opts+" )"
		} 
		
		//devcode
		dc = ""
		if (st_local("devcode")!="") {
			dc = " devcode("+st_local("devcode")+")"
		}
		
		stata("di ")
		stata(`"di as text "Fitting fixed effects model:""')
		stata("qui merlin "+ newcmd + ", nogen"+dc)
		
		//update with fixed effects
		binit[gml.initbindex] = st_matrix("e(b)")

	}
	
	if (st_local("from")!="") {
		binit = st_matrix(st_local("from"))
	}
	
	//pass in ap start values
	if (st_local("apstartvalues")!="") binit[gml.initapindex[1,]] = gml.initapindex[2,]
	
	//pass in vcv inits, in case restartvalues was specified
	if (st_local("restartvalues")!="") binit[gml.initvarindex[1,]] = gml.initvarindex[2,]

	stata("tempname from")
	st_matrix(st_local("from"),binit)
	
}

end
