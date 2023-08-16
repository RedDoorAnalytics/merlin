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

void merlin_build_touses(`gml' gml)
{
	//get any variable names for touses
	idvars 			= gml.levelvars[1,]
	gml.modeltouses = J(1,0,"")
	allvars			= J(gml.Nmodels,1,st_local("ifinvars"))

	for (i=1;i<=gml.Nmodels;i++) {
		
		cmps 	= merlin_get_indepvars(i)
		Ncmps 	= cols(cmps)

		//start with any level vars if there are some
		modelvars 	= idvars
		
		//extract all variable names from complex predictor
		
		for (j=1;j<=Ncmps;j++) {
			
			//prep
			dv = strtrim(cmps[j])
			_merlin_remove_at(dv)
			
			//one element, and its a variable
			pos = 1	
			if (!strpos(dv,"[") & !strpos(dv,"(") & !strpos(dv,"#")) {
				pos = 0
				if (!(gml.predict & dv==st_global("e(timevar"+strofreal(i)+")"))) {
					modelvars = modelvars,dv											
				}
			}
			
			//hash it up
			while (pos) {
				pos = strpos(dv,"#")
				if (pos) {
					dv2 = substr(dv,1,pos-1)
					dv 	= substr(dv,pos+1,.)
				}
				else dv2 = dv
				if (!strpos(dv2,"[") & !strpos(dv2,"(")) {							//variable name
					if (!(gml.predict & dv2==st_global("e(timevar"+strofreal(i)+")"))) {
						modelvars = modelvars,dv2								
					}
				}
				//rcs()/fp()/bs()/ps()/gp() varname
				if (strpos(dv2,"pc(") | strpos(dv2,"fp(") | strpos(dv2,"rcs(") | strpos(dv2,"bs(") | strpos(dv2,"ps(") | strpos(dv2,"gp(")) {	
					hasEV = strpos(dv2,"EV[")				//avoid fn(EV[],...)
					if (!hasEV) {		
						p1 	= strpos(dv2,"(")+1
						len = strpos(dv2,",") - p1
						tpvar = strtrim(substr(dv2,p1,len))
						if (!(gml.predict & tpvar==st_global("e(timevar"+strofreal(i)+")"))) {
							//if predicting and the var is equal to timevar then skip
							modelvars = modelvars,tpvar
						}
						//get offset/moffset
						len = strlen(dv2) - p1
						synt = "local 0 "+substr(dv2,p1,len)
						stata(synt)
						rc = _stata("syntax varname, [offset(varname) moffset(varname) *]")
						if (rc) exit(rc)
						
						if (st_local("offset")!="") {
							modelvars = modelvars,st_local("offset")
						}
						if (st_local("moffset")!="") {
							modelvars = modelvars,st_local("moffset")
						}
					}
				}
			}
		}

		//timevar
		//-> if not predicting as it can overide predicts timevar
		if (!gml.predict) 			modelvars = modelvars,st_local("timevar"+strofreal(i))
		//markout on predicts timevar, if there
		//unless obtaining a standardised prediction
		if (st_local("ptvar")!="" & st_local("standardise")=="") {
			modelvars = modelvars,st_local("ptvar")
		}

		//response vars
		//-> if survival model, can override predicts timevar
		//-> but needed if reffects specified
		predandsurv = gml.predict & st_local("failure"+strofreal(i))!="" & st_local("reffects")==""
		if (!gml.predict) {
			modelvars = modelvars,tokens(st_local("response"+strofreal(i)))
			//remove linterval if there, as missings are allowed
			if (gml.haslint[i]) {
				allvars[i] = allvars[i]+" "+modelvars[1,cols(modelvars)]	//still noted
				modelvars[1,cols(modelvars)] = ""
			}
		}

		//collect all variables
		allvars[i] = allvars[i]+" "+invtokens(modelvars)

		//post complete case model specific touse var
	
		modeltousei = "touse"+strofreal(i)
		stata("tempname "+modeltousei)
		
		//avoid model specific ifs/ins if being called from predict
		if (gml.predict) 	stata("mark "+st_local(modeltousei))
		else 				stata("mark "+st_local(modeltousei)+" "+st_local("if"+strofreal(i))+" "+st_local("in"+strofreal(i)))
		
		//markout vars, and update based on main touse (markout only marks out missing obs of vars)
		stata("markout "+st_local(modeltousei)+" "+invtokens(modelvars))
		stata("qui replace "+st_local(modeltousei)+"=0 if !"+gml.touse)
		
		gml.modeltouses = gml.modeltouses,st_local(modeltousei)
	}

	//to post for predict(ms)
	gml.allvars = invtokens(uniqrows(tokens(invtokens(allvars'))')')
	
	//if being called from predict(ms), update on N
	
		if (gml.predict) {
			for (i=1;i<=gml.Nmodels;i++) {		
				if (st_local("npredict")!="") {
					stata("qui replace "+gml.modeltouses[i]+" = _n <= "+st_local("npredict"))
				}
			}
		}

	//update global touse based on modeltouses, as models can have different # of obs, but everything must get read in
	//for possible indexing from other models
	
		for (i=1;i<=gml.Nmodels;i++) {		
			modeltousei = "touse"+strofreal(i)
			if (i==gml.Nmodels) st_local("touseup",st_local("touseup")+st_local(modeltousei)+"==1 ")
			else st_local("touseup",st_local("touseup")+st_local(modeltousei)+"==1 | ")
		}
		stata("qui replace "+gml.touse+" = ("+st_local("touseup")+")")
	
	stata("qui count if "+gml.touse+"==1")
	gml.N = st_numscalar("r(N)")
	
	//============================================================================================================//
	// note; 
	// everything above is based on complete case markouts
	//
	// imputed variables
		
		if (gml.hasImputed) {
			
			stata("tempvar imputetouse")
			stata("mark "+st_local("imputetouse"))

			//for each outcome model strip out imputed variable and markout
			vars = tokens(allvars[1])'
			vars = select(vars,vars:!=st_local("response"+strofreal(2)))
			stata("markout "+st_local("imputetouse")+" "+invtokens(vars'))
			//then markout on any variable in missing x predictor
			vars = tokens(allvars[2])'
			vars = select(vars,vars:!=st_local("response"+strofreal(2)))

			if (!(vars==J(0,0,""))) {
				stata("markout "+st_local("imputetouse")+" "+invtokens(vars'))
			}
			//then markout on missing x 
			stata("qui replace "+st_local("imputetouse")+" = 0 if "+st_local("response"+strofreal(2))+"!=.")
// 			stata("list y x "+gml.modeltouses[1]+" "+gml.modeltouses[2]+" "+st_local("imputetouse"))

			//overall touse index needs updating now	
			stata("qui replace "+gml.touse+" = 1 if "+st_local("imputetouse")+"==1")
			
		}
		
		
		
}

/*
- removes @ and subsequent text from a string
- if not there it's left unchanged
*/
void _merlin_remove_at(`SS' dv)
{
	if (strpos(dv,"@")) dv = substr(dv,1,strpos(dv,"@")-1)
}

end
