*! version 1.0.0 

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

void merlin_setup_check_clp(`gml' gml)
{
	//need to use gettoken with bind to handle brackets within rcs() etc which covers spaces
	for (i=1;i<=gml.Nmodels;i++) {
		gml.model = i
		st_local("rest",st_local("indepvars"+strofreal(i)))
		if (strpos(st_local("rest"),"##")) merlin_error("Stata's ## notation not supported")
		if (strpos(st_local("rest"),"*"))  merlin_error("* operator not supported")
		stata("gettoken lhs rest : rest, bind")						
		merlin_cmp_errorchecks(gml,strtrim(st_local("lhs")))
		while (st_local("rest")!="") {
			stata("gettoken lhs rest : rest, bind")
			merlin_cmp_errorchecks(gml,strtrim(st_local("lhs")))	
		}
	}
}

void merlin_cmp_errorchecks(`gml' gml,`SS' dv)
{
	//note; not effected by @'s
	
	//check each component
	pos 	= 1	
	varind 	= 0
	reind 	= 0
	while (pos) {
		pos = strpos(dv,"#")
		if (pos) {
			dv2 = substr(dv,1,pos-1)
			dv 	= substr(dv,pos+1,.)
		}
		else dv2 = dv
		merlin_dvcheck(gml,dv2,varind,reind)
	}

	if (reind>1) {
		errprintf("Only one random effect element is allowed in each component\n")
		exit(1986)
	}
	
}

void merlin_dvcheck(`gml' gml,`SS' dv2, `RS' varind, reind)
{

	if (strpos(dv2,"c.") | strpos(dv2,"i.")) merlin_error("Stata's . notation (e.g. c. or i.) in variable definitions is not currently supported")

	//check parentheses
	pos1 = strpos(dv2,"[")
	pos2 = strpos(dv2,"]")	
	if ((!pos1 & pos2) | (pos1 & !pos2)) {
		errprintf("Missing [ or ]\n")
		exit(1986)
	}

	pos5 = strpos(dv2,"(")
	pos6 = strpos(dv2,")")	
	if ((!pos5 & pos6) | (pos5 & !pos6)) {
		errprintf("Missing ( or )\n")
		exit(1986)
	}
	
	if (strpos(dv2,"@")) {														//check at
		at = substr(dv2,strpos(dv2,"@")+1,.)
		if (strtoreal(at)==.) merlin_error("@"+at+" must be a real number")
		dv2 = substr(dv2,1,strpos(dv2,"@")-1)									//strip of @ to continue with
	}
	
	hasrcs	= substr(dv2,1,4)=="rcs(" 
	
	if (!hasrcs & pos1 & pos2) { 												//random effect or ?EV[outcome]

		noEV 	= substr(dv2,1,2)!="EV" & substr(dv2,1,3)!="dEV" & substr(dv2,1,4)!="d2EV" & substr(dv2,1,3)!="iEV"
		noEV	= noEV & substr(dv2,1,2)!="XB" & substr(dv2,1,3)!="dXB" & substr(dv2,1,4)!="d2XB" & substr(dv2,1,3)!="iXB"
		if (noEV) {																//random effect
			merlin_check_name(dv2)												//checks RE name
			merlin_check_pipes(dv2,pos1,pos2)		
			merlin_confirm_idvars(dv2,pos1,pos2)	
			reind++
		}
		else {																	//EV,dEV,d2EV,iEV
			merlin_confirm_depvar(gml,dv2,pos1,pos2)								//check response variable
		}
		
	}
	else if (pos5 & pos6) {														//rcs()/fp()/mf()
		
		hasfp	= substr(dv2,1,3)=="fp(" 
		hasmf	= substr(dv2,1,3)=="mf(" 
		hasbs	= substr(dv2,1,3)=="bs(" 
		hasps	= substr(dv2,1,3)=="ps(" 
		haspc	= substr(dv2,1,3)=="pc(" 
		hasgp	= substr(dv2,1,3)=="gp(" 
		hasfn 	= hasrcs | hasfp | hasmf | hasbs | hasps | hasgp | haspc
		if (!hasfn) merlin_error("Invalid function "+dv2)
		
		pos1 	= strpos(dv2,"(") + 1
		len 	= strlen(dv2) - pos1											//can't use ) as it's in rcs(...df()... )
		
		if (hasfp) {															//fp syntax check
			
			synt = substr(dv2,pos1,len)
			if (strpos(synt,",")) 	stata("local 0 "+synt)
			else 					stata("local 0 , "+synt)
			stata("syntax varlist(min=1 max=1) , POWers(numlist max=2) [OFFset(varname) MOFFset(varname)]")
			
			powers 	= strtoreal(tokens(st_local("powers")))						//check fp powers
			nfp 	= cols(powers)
			if (nfp>2) merlin_error("Invalid "+dv2)
			
			if (st_local("offset")!="" | st_local("moffset")!="") {
				gml.offsetflag 	= gml.offsetflag + 1
			}
			
			allpows = (-2,-1,-0.5,0,0.5,1,2,3)
			for (i=1;i<=nfp;i++) {
				if (!sum(powers[1,i]:==allpows)) merlin_error("Invalid "+dv2)
			}
			
		}
		else if (hasrcs) {															//rcs syntax check
			synt = substr(dv2,pos1,len)
			if (strpos(synt,",")) stata("local 0 "+synt)
			else stata("local 0 , "+synt)
			stata("syntax anything , [DF(string) KNOTS(numlist asc) LOG ORTHog EVent OFFset(varname) MOFFset(varname)]")
			
			anything = st_local("anything")
			hasrcsEV = strpos(anything,"EV[")
			
			if (st_local("offset")!="" | st_local("moffset")!="") {
				gml.offsetflag 	= gml.offsetflag + 1
			}
			
			if (st_local("df")!="") {
				if (st_local("knots")!="") {
					merlin_error("Cannot specify both df() and knots() within rcs()")
				}
				df = strtoreal(st_local("df"))
				if (df<1 | df>10) 	merlin_error("df() must be between 1 and 10")
				if (df!=round(df)) 	merlin_error("df() must be an integer")
			}
			else {
				if (st_local("knots")=="") {
					merlin_error("One of df() or knots() needed in rcs()")
				}
			}
		
		}
		else if (haspc) {															//pc syntax check
			synt = substr(dv2,pos1,len)
			if (strpos(synt,",")) stata("local 0 "+synt)
			else stata("local 0 , "+synt)
			stata("syntax anything , KNOTS(numlist asc min=1) [OFFset(varname) MOFFset(varname) NOREFerence]")
			
			if (st_local("offset")!="" | st_local("moffset")!="") {
				gml.offsetflag 	= gml.offsetflag + 1
			}
			
			if (gml.familys[gml.model]=="rp" | gml.familys[gml.model]=="logchazard") {
				merlin_error("pc() elements can't currently be used with log cumulative hazard scale models")
			}
		}
		else if (hasbs) {															//rcs syntax check
			synt = substr(dv2,pos1,len)
			if (strpos(synt,",")) stata("local 0 "+synt)
			else stata("local 0 , "+synt)
			stata("syntax varlist(min=1 max=1) , [Order(string) Degree(string) DF(string) Knots(numlist asc) BKnots(numlist asc min=2 max=2) EVent LOG INTercept OFFset(varname) MOFFset(varname)]")
			
			if (st_local("offset")!="" | st_local("moffset")!="") {
				gml.offsetflag 	= gml.offsetflag + 1
			}
			
			if (st_local("order")!="" & st_local("degree")!="") {
				merlin_error("order() and degree() can't both be specified")
			}
			
			if (st_local("order")!="") {
				o = strtoreal(st_local("order"))
				if (o<2 | trunc(o)!=o) merlin_error("order() must be an integer > 1")
			}
		
			if (st_local("degree")!="") {
				o = strtoreal(st_local("degree"))
				if (o<1 | trunc(o)!=o) merlin_error("degree() must be an integer >= 1")
			}
		
			if (st_local("df")!="") {
				if (st_local("knots")!="") {
					merlin_error("Cannot specify both df() and knots() within rcs()")
				}
				df = strtoreal(st_local("df"))
				if (df<1 | df>10) 	merlin_error("df() must be between 1 and 10")
				if (df!=round(df)) 	merlin_error("df() must be an integer")
			}
			
			if (gml.familys[gml.model]=="rp" | gml.familys[gml.model]=="logchazard") {
				merlin_error("bs() elements can't currently be used with log cumulative hazard scale models")
			}
			
		}
		else if (hasps) {															//rcs syntax check
			synt = substr(dv2,pos1,len)
			if (strpos(synt,",")) stata("local 0 "+synt)
			else stata("local 0 , "+synt)
			stata("syntax varlist(min=1 max=1) , [Order(string) Degree(string) DF(string) LOG INTercept OFFset(varname) MOFFset(varname)]")
			
			if (st_local("offset")!="" | st_local("moffset")!="") {
				gml.offsetflag 	= gml.offsetflag + 1
			}
			
			if (st_local("order")!="" & st_local("degree")!="") {
				merlin_error("order() and degree() can't both be specified")
			}
			
			if (st_local("order")!="") {
				o = strtoreal(st_local("order"))
				if (o<2 | trunc(o)!=o) merlin_error("order() must be an integer > 1")
			}
		
			if (st_local("degree")!="") {
				o = strtoreal(st_local("degree"))
				if (o<1 | trunc(o)!=o) merlin_error("degree() must be an integer >= 1")
			}
		
			if (st_local("df")!="") {
				df = strtoreal(st_local("df"))
				if (df<1) 			merlin_error("df() must be >= 1")
				if (df!=round(df)) 	merlin_error("df() must be an integer")
			}
			
		}
		else if (hasgp) {														//gaussian process
			synt = substr(dv2,pos1,len)
			if (strpos(synt,",")) 	stata("local 0 "+synt)
			else 					stata("local 0 , "+synt)
			stata("syntax varlist(min=1 max=1) , []")
		}
		
	}
	else {																		//variable
		merlin_confirmvars(dv2)				
		varind++
	}
}

void merlin_check_varname(`SS' dv)
{
	pos = strpos(dv,"#")					
	if (pos) {
		merlin_confirmvars(substr(dv,1,pos-1))	
		dv = substr(dv,pos+1,.)
	}	
}

void merlin_check_name(`SS' dv)
{	
	`SS' rest
	`RS' hasm
	
	hasm 	= substr(dv,1,1)=="M"
	
	if (!hasm) {
		errprintf("invalid latent variable specification;\n")
		errprintf("random effect names must start with M\n")
		exit(1986)	
	}
	num = substr(dv,2,1)
	if (strtoreal(num)==.) {
		errprintf("invalid latent variable specification;\n")
		errprintf("M must be followed by an integer\n")
		exit(1986)	
	}
}

void merlin_check_pipes(`SS' dv, `RS' pos1, `RS' pos2)
{
	`SS' newvar
	
	newvar = substr(dv,pos1+1,pos2-pos1-1)
	if (strpos(newvar,"<") & strpos(newvar,">")) {
		errprintf("invalid latent variable specification;\n")
		errprintf("both '<' and '>' detected within the same specification\n")
		exit(1986)
	}
}

void merlin_confirm_idvars(`SS' dv, `RS' pos1, `RS' pos2)
{
	`SS' newvar,copyvar, tempidvar
	`RS' pipepos
	
	ind1 = ind2 = 0
	newvar = substr(dv,pos1+1,pos2-pos1-1)
	copyvar = newvar
	pipepos = strpos(copyvar,"<")
	while (pipepos) {
		ind1 = 1
		tempidvar = substr(copyvar,1,pipepos-1)
		merlin_confirmvars(tempidvar)
		copyvar = substr(copyvar,pipepos+1,.)
		pipepos = strpos(copyvar,"<")
	}
	copyvar = newvar
	pipepos = strpos(copyvar,">")
	while (pipepos) {
		ind2 = 1
		tempidvar = substr(copyvar,1,pipepos-1)
		merlin_confirmvars(tempidvar)
		copyvar = substr(copyvar,pipepos+1,.)
		pipepos = strpos(copyvar,">")
	}	
	if (ind1 & ind2) {
		errprintf("Can't have both > and < in level specifier\n")
		exit(1986)
	}
}

void merlin_confirm_depvar(`gml' gml,dv,pos1,pos2)	
{

	`SS' yvar
	yvar = substr(dv,pos1+1,pos2-pos1-1)
	//strip of , options if there
	if (strpos(yvar,",")) {
		op = strtrim(substr(yvar,strpos(yvar,",")+1,.))
		if (op!="impute") {
			merlin_error("Invalid option in EV[]")
		}
		yvar = substr(yvar,1,strpos(yvar,",")-1)
	}

	//check yvar is a response variable for a model which isn't the current, stored in gml.model
	if (gml.familys[gml.model]=="null") 	resp = "_merlin_null"
	else 									resp = tokens(st_local("response"+strofreal(gml.model)))[1]		//accounts for survival

	if (yvar==resp | yvar==strofreal(gml.model)) {
		errprintf("You can't include the expected value of an outcome as a dependent variable in its own model\n")
		exit(1986)
	}
	flag = 1

	for (i=1;i<=gml.Nmodels;i++) {
		if (i!=gml.model) {
			if (gml.familys[i]=="null") 	resp = "_merlin_null"
			else							resp = tokens(st_local("response"+strofreal(i)))[1]		//accounts for survival
			if (yvar==resp | yvar==strofreal(i) | gml.familys[i]=="null") {
				flag 	= 0
				k 		= i
			}
		}
	}

	if (flag) merlin_error("Error; "+yvar+" not a response variable/model")
	
	merlin_confirm_expval(gml,k)

}

void merlin_confirmvars(`SR' var1)
{
	//strip off @ if there
	if (strpos(var1,"@")) {
		var1 = substr(var1,1,strpos(var1,"@")-1)
	}
	
	for (i=1;i<=cols(var1);i++) {
		stata("capture confirm numeric variable "+var1[1,i]+", exact")
		stata("local err = _rc")
		er = strtoreal(st_local("err"))
		if (er) merlin_error("variable "+var1[1,i]+" not found")
	}
}

void merlin_confirm_expval(`gml' gml, `RS' k)
{
	f = gml.familys[k]
	if (f=="gompertz" | f=="lquantile" | f=="loghazard" | f=="user" | f=="ordinal") {
		merlin_error("EV[]/XB[] of a family("+f+") not currently supported")
	}
}

end
