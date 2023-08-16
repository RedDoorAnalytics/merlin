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

void merlin_build_cmps(`gml' gml)
{
	/*
	(response cmp1 cmp2 ..., )

	-> Each cmpi can be

	el1#el2#...
	M1[id]@eqnname
	MV1#[id]@eqnname

	--> Each elj can be
	
	varname						- variable name											//done
	M?[] 						- random effect											//done
	fp()						- fractional polynomial									//done
	rcs()						- rcs function of time									//done
	mf()						- mata function of time									//done
	EV[] 						- expected value of outcome								//done
	dEV[] 						- first derivative of expected value of outcome			//done
	d2EV[] 						- second derivative of expected value of outcome		//done
	iEV[]						- integral of expected value of outcome over t			//done
	XB[] 						- expected value of CP									//done
	dXB[] 						- first derivative of expected value of CP				//done
	d2XB[] 						- second derivative of expected value of CP				//done
	iXB[]						- integral of expected value of CP over t				//done
	*/
	
	gml.Ncmps 		= J(gml.Nmodels,1,0)
	gml.elindex 	= asarray_create("real",2)
	gml.elinfo 		= asarray_create("real",3)
	gml.Nels 		= asarray_create("real",1)
	
	gml.eqnindex = J(gml.Nmodels,2,1)
	
	gml.initbindex	= J(1,0,.)
	
	eqn = 1
	
	for (mod=1;mod<=gml.Nmodels;mod++) {
		gml.model = gml.modtoind = mod
		gml.eqnindex[mod,1] = eqn
		
		depvars = merlin_get_indepvars(mod)
		if (depvars!="") gml.Ncmps[mod] = cols(depvars)
		Nelsmat = todospec = hascoef = J(0,1,.)

		for (i=1;i<=gml.Ncmps[mod];i++) {
			eltype 	= J(0,1,.)
			Nels 	= 0														//get Nels per cmp
			
			dv 		= strtrim(depvars[i])
			
			//strip off and store @
			atpos	= strpos(dv,"@")
			if (atpos) {
				at = substr(dv,atpos+1,.)														//get constraint
				dv = substr(dv,1,atpos-1)														//strip off @ if there
			}
			
			//hash up the component
			pos 	= 1	
			while (pos) {
				pos = strpos(dv,"#")
				if (pos) {
					dv2 = substr(dv,1,pos-1)
					dv 	= substr(dv,pos+1,.)
				}
				else dv2 = dv
				if (dv2!="") {  
					Nels++
					eltype = eltype\merlin_build_els(gml,mod,i,Nels,dv2)
				}
			}
			Nelsmat = Nelsmat\Nels
			asarray(gml.elindex,(mod,i),eltype)

			if (atpos) 	merlin_build_variables(gml,mod,i,eltype,eqn,at)
			else		merlin_build_variables(gml,mod,i,eltype,eqn)
			
		}
		//constant
		if (gml.hascons[mod]) {
			gml.initbindex = gml.initbindex,eqn
			st_local("xb"+strofreal(mod),st_local("xb"+strofreal(mod)) + " /cons"+strofreal(mod))		
			eqn++
		}
		
		gml.eqnindex[mod,2] = eqn-1
		
		//distributional ancillary parameters
		if (gml.familys[mod]=="ordinal") {
			if (gml.links[mod]=="logit") 	omod = "qui ologit"
			else 							omod = "qui oprobit"
			stata(omod+" "+st_local("response"+strofreal(mod)))						//extra starting values for ordinal repsonses
			ob = st_matrix("e(b)")	
			for (dap=1;dap<=gml.Ndistancp[mod];dap++) {
				gml.initbindex = gml.initbindex,eqn
				gml.initdistapindex = gml.initdistapindex,(eqn\ob[dap])
				eqn++
			}		
		}
		else {
			for (dap=1;dap<=gml.Ndistancp[mod];dap++) {
				gml.initbindex = gml.initbindex,eqn
				
				//special case for rp first spline term starting value
				if (dap==1 & gml.familys[mod]=="rp") {
					gml.initdistapindex = gml.initdistapindex,(eqn\1)
				}
				
				eqn++
			}		
		}
		
		//ancillary parameters
		for (dap=1;dap<=gml.Nap[mod];dap++) {
			gml.initbindex = gml.initbindex,eqn
			gml.initapindex = gml.initapindex,(eqn\strtoreal(st_local("apstartvalues")))
			eqn++
		}
		
		//store
		asarray(gml.Nels,mod,Nelsmat)
		
	}
	gml.Nb = eqn-1
}

/*
merlin_build_els()

gives return code for type of element:
	1 - variable
	2 - random effect
	3 - time/variable function
	4 - EV
	5 - iEV
	6 - dEV
	7 - d2EV
*/

`RS' merlin_build_els(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' dv2)
{
	el 			= dv2
	hassquareb	= strpos(el,"[")

	if (hassquareb) {
	
		if (strpos(el,"EV")) 			return(merlin_setup_expval(gml,mod,i,Nels,el))	//?EV
		else if (strpos(el,"XB")) {														//?XB
			
			kpos 	= strpos(el,"[")													//get response varname
			kpos2 	= strpos(el,"]")
			y 		= substr(el,kpos+1,kpos2-kpos-1)
			if (strpos(y,",")) y = substr(y,1,strpos(y,",")-1)							//strip off , options if there

			for (k=1;k<=gml.Nmodels;k++) {												//get response varname index and store
				if (st_local("response"+strofreal(k))==y | strofreal(k)==y) {
					asarray(gml.elinfo,(mod,i,Nels),k)									//store info
				}
			}
			
			if (substr(el,1,2)=="XB")		return(10)
			else if (strpos(el,"iXB")) 		return(11)
			else if (strpos(el,"dXB")) 		return(12)
			else if (strpos(el,"d2XB")) 	return(13)
			
		}
		else 						return(merlin_setup_re(gml,mod,i,Nels,el))			//random effect	
			
	}
	else {	
		
		if 		(strpos(el,"mf("))  return(merlin_setup_mf(gml,mod,i,Nels,el))			//user-defined function
		else if (strpos(el,"rcs(")) return(merlin_setup_rcs(gml,mod,i,Nels,el))			//rcs function
		else if (strpos(el,"fp(")) 	return(merlin_setup_fp(gml,mod,i,Nels,el))			//fp function
		else 						return(merlin_setup_var(gml,mod,i,Nels,el))			//varname

	}
}

void merlin_build_variables(`gml' gml, `RS' mod, `RS' i, `RC' eltype, `RS' eqn, | `SS' at)
{
	
	//build core variables to post to ml equations
	
	Nels 	= rows(eltype)
	Nobs	= gml.Nobs[gml.Nlevels,mod]

	newvars = J(Nobs,1,1)
	
	initcns = 0
	
	for (el=1;el<=Nels; el++) {

		if 		(eltype[el]==1)	{
		
			nextvars = asarray(gml.elinfo,(mod,i,el))[asarray(gml.xbindex,mod)]
			
		}
		else if (eltype[el]==2) {												//random effects
			
			nextvars = J(Nobs,1,1)
			initcns = 1
			
		}
		else if (eltype[el]==3) {												//mf()
			
			nextvars = J(Nobs,1,1)
			
		}
		else if (eltype[el]==4 | eltype[el]==5 | eltype[el]==6 | eltype[el]==7) {	//?EV
			
			nextvars = J(Nobs,1,1)
			initcns = 1
			
		}
		else if (eltype[el]==10 | eltype[el]==11 | eltype[el]==12 | eltype[el]==13) {	//?XB
			
			nextvars = J(Nobs,1,1)
			initcns = 1
			
		}
		else if (eltype[el]==8) {												//rcs()
			
			nextvars = merlin_xz_rcs(gml,i,el,0)
			
		}
		else if (eltype[el]==9) {												//fp()
			
			nextvars = merlin_xz_fp(gml,i,el)

		}
		
		//rebuild
		Nold = cols(newvars)
		Nnew = cols(nextvars)

		copyold = newvars
		newvars = J(Nobs,0,.)
		for (j=1;j<=Nnew;j++) {
			newvars = newvars,(copyold :* nextvars[,j])
		}
	
	}

	//post to Stata and ml equation local
	stub = "_cmp_"+strofreal(mod)+"_"+strofreal(i)
	if (!gml.nogen) stata("cap drop "+stub+"_*")
	Nvars = cols(newvars)
	
	gml.E.Nvars[mod] = gml.E.Nvars[mod]+" "+strofreal(Nvars)
	
	names = J(1,0,"")
	eqnames = J(1,0,"")
	for (r=1;r<=Nvars;r++) {
		if (!initcns) gml.initbindex = gml.initbindex,eqn
		eqn++
		names = names,(stub+"_"+strofreal(r))
		eqnames = eqnames,("/"+stub+"_"+strofreal(r))
		//at constraints
		if (args()==6) {
			stata("constraint free")
			stata("constraint "+st_global("r(free)")+"["+stub+"_"+strofreal(r)+"][_cons] = "+at)
			stata("local constraints "+st_local("constraints")+" "+st_global("r(free)"))			
		}
	}

	if (!gml.nogen) {
		
		build = sum(eltype:==3) | sum(eltype:==8) | sum(eltype:==9)
		if (build) {
			id = st_addvar("double",names)
			st_store(.,id,gml.modeltouses[mod],newvars)
			printf("variables created for model "+strofreal(mod)+", component "+strofreal(i)+": "+stub+"_1 to "+stub+"_"+strofreal(Nvars)+"\n")
		}
		
	}
	
	//post for xb equations
	st_local("xb"+strofreal(mod),st_local("xb"+strofreal(mod)) + " "+invtokens(eqnames))		
}

void merlin_get_cmps_labels(`gml' gml)
{
	gml.cmplabels = J(gml.Nmodels,1,"")
	
	for (mod=1;mod<=gml.Nmodels;mod++) {
		
		depvars = merlin_get_indepvars(mod)
		
		for (i=1;i<=gml.Ncmps[mod];i++) {
			
			dv 		= strtrim(depvars[i])
			
			//strip off and store @
			atpos	= strpos(dv,"@")
			if (atpos) {
				at = substr(dv,atpos+1,.)														//get constraint
				dv = substr(dv,1,atpos-1)														//strip off @ if there
			}
			
			//hash up the component
			pos 	= 1	
			first 	= 1
			while (pos) {
				pos = strpos(dv,"#")
				if (pos) {
					dv2 = substr(dv,1,pos-1)
					dv 	= substr(dv,pos+1,.)
				}
				else dv2 = dv
				
				if (dv2!="") {  
					if (first) {
							gml.cmplabels[mod] = gml.cmplabels[mod] + " " + merlin_cmp_label(dv2)
							first = 0
					}
					else	gml.cmplabels[mod] = gml.cmplabels[mod] + "#" + merlin_cmp_label(dv2)
				}
			}
			
		}
		
	}

}

`SS' merlin_cmp_label(`SS' el)
{
	hassquareb	= strpos(el,"[")

	if (hassquareb) {
	
		if (strpos(el,"EV")) 	{
			if (substr(el,1,2)=="EV")		return("EV[]")
			else if (strpos(el,"iEV")) 		return("iEV[]")
			else if (strpos(el,"dEV")) 		return("dEV[]")
			else if (strpos(el,"d2EV")) 	return("d2EV[]")
		}
		else if (strpos(el,"XB")) {	
			if (substr(el,1,2)=="XB")		return("XB[]")
			else if (strpos(el,"iXB")) 		return("iXB[]")
			else if (strpos(el,"dXB")) 		return("dXB[]")
			else if (strpos(el,"d2XB")) 	return("d2XB[]")
		}
		else {
			return(el)												//re
		}
			
	}
	else {	
		
		if 		(strpos(el,"mf("))  return("mf()")			//user-defined function
		else if (strpos(el,"rcs(")) return("rcs()")			//rcs function
		else if (strpos(el,"fp(")) 	return("fp()")			//fp function
		else 						return(el)				//varname

	}
	
}

end
