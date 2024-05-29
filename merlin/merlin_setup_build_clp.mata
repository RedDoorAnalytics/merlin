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

void merlin_build_clp(`gml' gml)
{
	merlin_build_xz(gml)			
	if (sum(gml.familys:=="rp") | sum(gml.familys:=="logchazard")) {
		merlin_build_dxz(gml)
	}
}

void merlin_build_xz(`gml' gml)
{
	/*
	(response cmp1 cmp2 ..., )

	-> Each cmpi can be
	el1#el2#...@#

	-> Each elj can be
	varname				- variable name										
	M?[] 				- random effect				
	fp()				- fractional polynomial
	rcs()				- rcs function of time	
	?EV[]/?XB[]			- other outcome function
	mf()				- mata function of time
	gp()				- gaussian process
	
	-> design matrix gets built and stored for all elements except random effects and ?EV[]/?XB[]
	*/
	
	gml.X 			= asarray_create("real",1)		//store model specific design matrices
	gml.X_bindex 		= asarray_create("real",2)		//indexes for design matrix and parameter vector, split by not [] or yes
	
	gml.E.Nvars 		= J(gml.Nmodels,1,"")
	gml.Ncmps 		= J(gml.Nmodels,1,0)			//RC - # of components for each model
	gml.CmpXBIndex 		= asarray_create("real",2)		//min and max of design matrix and coef vector, for each component
	gml.eqnindex 		= J(gml.Nmodels,2,1)			//equation min and max index for clp, for each model
	
	gml.elindex 		= asarray_create("real",2)		//RC of element types, for each model,cmp
	gml.elinfo 		= asarray_create("real",3)
	gml.eltvar 		= asarray_create("real",2)
	gml.Nels 		= asarray_create("real",1)
	gml.elvarname		= asarray_create("real",3)		
	
	gml.hasconstraint 	= asarray_create("real",1)		//is each coeff of each comp/el constrained - for gf1/2
	
	gml.initbindex		= J(1,0,.)				//index for coefficients to estimate in initial values fit
	
	beqn = 1
	for (mod=1;mod<=gml.Nmodels;mod++) {
		
		gml.model = gml.modtoind 	= mod
		Nobs				= gml.N 			//obs[gml.Nlevels,mod]
		depvars 			= merlin_get_indepvars(mod)
		if (depvars!="") {
			gml.Ncmps[mod] 		= cols(depvars)
		}
		hasconstraint 			= J(1,0,.)
		
		xeqn 				= 1
		X 				= J(Nobs,0,.)
		X1index = X2index 		= J(2,0,.)			//first row for X, second for betas (needed for multi-outcome models)
		
		gml.eqnindex[mod,1] 		= beqn
		cmpxindex = cmpbindex 		= J(0,2,.)
		
		Nelsmat = todospec = hascoef = J(0,1,.)

		for (c=1;c<=gml.Ncmps[mod];c++) {

			include_X 	= 0
			cmpxindex 	= cmpxindex\xeqn,.
			cmpbindex 	= cmpbindex\beqn,.
			eltype 		= J(0,1,.)
			eltvar		= J(0,2,.)
			elsyntax 	= J(0,1,"")
			Nels 		= 0			                //get Nels per cmp
			dv 		= strtrim(depvars[c])
			
			//strip off and store @
			atpos	= strpos(dv,"@")
			if (atpos) {
				at = substr(dv,atpos+1,.)			//get constraint
				dv = substr(dv,1,atpos-1)			//strip off @ if there
			}
			else at = ""
			
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
					eltype 		= eltype\merlin_get_element_codes(gml,mod,dv2)
					eltvar		= eltvar\merlin_get_element_istimevar(gml,mod,dv2)
					include_X  	= include_X + merlin_get_X_index_flag(dv2)
					elsyntax 	= elsyntax\dv2
				}
			}
			
			asarray(gml.elindex,(mod,c),eltype)
			asarray(gml.eltvar,(mod,c),eltvar)
			gml.istimedep[mod,] = gml.istimedep[mod,] + colsum(eltvar)
			Nelsmat = Nelsmat\Nels
			
			//update design matrix
			refbeqn = beqn
			X 		= X,merlin_build_variables(gml,mod,c,eltype,elsyntax,xeqn,beqn,at)

			cmpxindex[c,2] 	= xeqn-1
			cmpbindex[c,2] 	= beqn-1
			if (!include_X) X1index = X1index,((cmpxindex[c,1]..cmpxindex[c,2])\(cmpbindex[c,1]..cmpbindex[c,2]))
			else		X2index = X2index,((cmpxindex[c,1]..cmpxindex[c,2])\(cmpbindex[c,1]..cmpbindex[c,2]))

			//update constraint flag for new variables
			atconstraint 	= 0
			if (atpos) atconstraint = 1
			hasconstraint 	= hasconstraint,J(1,(beqn-refbeqn),atconstraint)
		}

		//constant
		if (gml.hascons[mod]) {
			cmpxindex 	= cmpxindex\xeqn,xeqn
			cmpbindex 	= cmpbindex\beqn,beqn
			X 		= X,J(Nobs,1,1)
			X1index 	= X1index,(xeqn\beqn)
			gml.initbindex 	= gml.initbindex,beqn
			hasconstraint 	= hasconstraint,0
			xeqn++
			beqn++
			st_local("xb"+strofreal(mod),st_local("xb"+strofreal(mod)) + " (cons"+strofreal(mod)+":)")	
		}

		asarray(gml.X,mod,X)				//store design matrix
		asarray(gml.X_bindex,(mod,1),X1index)		//index for parameters to match design matrix (not re's, ?[], etc)
		asarray(gml.X_bindex,(mod,2),X2index)		//index for parameters to match design matrix (re's, ?[], etc)
		asarray(gml.CmpXBIndex,(mod,1),cmpxindex)	//store indices for clp design matrix
		asarray(gml.CmpXBIndex,(mod,2),cmpbindex)	//store indices for clp betas
		asarray(gml.hasconstraint,mod,hasconstraint)	//store constraint flags
		gml.eqnindex[mod,2] = beqn-1

		//handle distributional ancillary parameters
		if (gml.familys[mod]=="ordinal" & !gml.predict) {
			if (gml.links[mod]=="logit") 	omod = "qui ologit"
			else 				omod = "qui oprobit"
			stata(omod+" "+st_local("response"+strofreal(mod)))	//extra starting values for ordinal repsonses
			ob = st_matrix("e(b)")	
			for (dap=1;dap<=gml.Ndistancp[mod];dap++) {
				gml.initbindex          = gml.initbindex,beqn
				gml.initdistapindex     = gml.initdistapindex,(beqn\ob[dap])
				beqn++
			}
		}
		else {
			for (dap=1;dap<=gml.Ndistancp[mod];dap++) {
				gml.initbindex = gml.initbindex,beqn
				//special case for rp first spline term starting value
				if (gml.familys[mod]=="rp" | gml.familys[mod]=="aft") {
					if (dap==1)     gml.initdistapindex = gml.initdistapindex,(beqn\1)
					else 	        gml.initdistapindex = gml.initdistapindex,(beqn\0)
				}
				//special case for gompertz shape parameter
				if (gml.familys[mod]=="gompertz") {
					gml.initdistapindex = gml.initdistapindex,(beqn\0.1)
				}
				//special case for ggamma
				if (gml.familys[mod]=="ggamma") {
					gml.initdistapindex = gml.initdistapindex,(beqn\1)
				}
				beqn++
			}
		}
		
		//ancillary parameters
		for (dap=1;dap<=gml.Nap[mod];dap++) {
			gml.initbindex 	= gml.initbindex,beqn
			gml.initapindex = gml.initapindex,(beqn\strtoreal(st_local("apstartvalues")))
			beqn++
		}
		
		//store
		asarray(gml.Nels,mod,Nelsmat)
		
	}
	gml.Nb = beqn-1
}

/*
merlin_build_els()
-> get return code for type of element:
	1 - variable
	2 - random effect
	3 - time/variable function
	4 - EV
	5 - iEV
	6 - dEV
	7 - d2EV
	8 - rcs()
	9 - fp()
	10 - XB
	11 - iXB
	12 - dXB
	13 - d2XB
	14 - bs()
	15 - ps()
	16 - gp()
	17 - pc()
*/

`RS' merlin_get_element_codes(`gml' gml, `RS' mod, `SS' dv2)
{
	el 			= dv2
	hassquareb	= strpos(el,"[")
	hasroundb	= strpos(el,"(")
	
	if (hassquareb & !hasroundb) {
		if (strpos(el,"EV")) {
			if (substr(el,1,2)=="EV")	return(4)
			else if (strpos(el,"iEV")) 	return(5)
			else if (strpos(el,"dEV")) 	return(6)
			else if (strpos(el,"d2EV")) 	return(7)
		}
		else if (strpos(el,"XB")) {
			if (substr(el,1,2)=="XB")	return(10)
			else if (strpos(el,"iXB")) 	return(11)
			else if (strpos(el,"dXB")) 	return(12)
			else if (strpos(el,"d2XB")) 	return(13)
		}
		else 					return(2)	//random effect
	}
	else {	
		if 	(strpos(el,"mf("))  	        return(3)	//user-defined function
		else if (strpos(el,"rcs(")) 		return(8)	//rcs function
		else if (strpos(el,"fp(")) 		return(9)	//fp function
		else if (strpos(el,"bs(")) 		return(14)	//bs function
		else if (strpos(el,"pc(")) 		return(17)	//pc function
		else 					return(1)	//varname
	}
}

//flag whether varname, or function input var is timevar()
`RR' merlin_get_element_istimevar(`gml' gml, `RS' mod, `SS' dv2)
{
	el 		= dv2
	hassquareb	= strpos(el,"[")
	hasroundb	= strpos(el,"(")
	
	if (hassquareb & !hasroundb) 	varname = substr(el,hassquareb+1,strpos(el,"]")-hassquareb+1)
	else if (hasroundb) 		varname = substr(el,hasroundb+1,strpos(el,",")-hasroundb-1)
	else 				varname = el
	if (strpos(el,"mf("))		varname = gml.tvarnames[mod] 		//override for mf()

	return(varname==gml.tvarnames[mod],varname==gml.ltruncated[mod])
}

//whether to include an element stored design matrix
`RS' merlin_get_X_index_flag(`SS' dv2)
{
	hassquareb	= strpos(dv2,"[")
	hasroundb	= strpos(dv2,"(")
	if (hassquareb & !hasroundb) 	return(1)
	else							return(0)
}

`RM' merlin_build_variables(`gml' gml, `RS' mod, `RS' c, `RC' eltype, `SC' elsyntax, `RS' xeqn, `RS' beqn, `SS' at)
{
	
	//build core variables to post to ml equations
	
	Nels 	= rows(eltype)
	Nobs	= gml.N //obs[gml.Nlevels,mod]
	newvars = J(Nobs,1,1)
	noinit	= 0

	for (el=1;el<=Nels; el++) {
		if 	(eltype[el]==1)		nextvars = merlin_setup_var(gml,mod,c,el,elsyntax[el])
		else if (eltype[el]==8) 	nextvars = merlin_setup_rcs(gml,mod,c,el,elsyntax[el])
		else if (eltype[el]==9) 	nextvars = merlin_setup_fp(gml,mod,c,el,elsyntax[el])
		else if (eltype[el]==14) 	nextvars = merlin_setup_bs(gml,mod,c,el,elsyntax[el])
		else if (eltype[el]==17) 	nextvars = merlin_setup_pc(gml,mod,c,el,elsyntax[el])
		else if (eltype[el]==3) 	nextvars = merlin_setup_mf(gml,mod,c,el,elsyntax[el])
		else if (eltype[el]==2) {
			nextvars = merlin_setup_re(gml,mod,c,el,elsyntax[el])
			noinit = 1		//don't include in fixed effects starting value fit
		}
		else if (eltype[el]==4 | eltype[el]==5 | eltype[el]==6 | eltype[el]==7 ) {
			nextvars = merlin_setup_expval(gml,mod,c,el,elsyntax[el])
			noinit = 1
		}

		//rebuild
		if (Nels>1) {
			Nnew 	= cols(nextvars)
			copyold = newvars
			newvars = J(Nobs,0,.)
			for (j=1;j<=Nnew;j++) newvars = newvars,(copyold :* nextvars[,j])
		}
		else newvars = nextvars
	}

	//post to Stata and ml equation local
	stub = "_cmp_"+strofreal(mod)+"_"+strofreal(c)
	if (!gml.nogen) stata("cap drop "+stub+"_*")
	Nvars = cols(newvars)
	
	gml.E.Nvars[mod] = gml.E.Nvars[mod]+" "+strofreal(Nvars)
	
	names 	= J(1,0,"")
	eqnames = J(1,0,"")

	for (r=1;r<=Nvars;r++) {
		if (!noinit) gml.initbindex = gml.initbindex,beqn
		names 	= names,(stub+"_"+strofreal(r))
		eqnames = eqnames,("("+stub+"_"+strofreal(r)+": )")
		//at constraints
		if (at!="") {
			stata("constraint free")
			stata("constraint "+st_global("r(free)")+"["+stub+"_"+strofreal(r)+"][_cons] = "+at)
			stata("local constraints "+st_local("constraints")+" "+st_global("r(free)"))			
		}
		xeqn++	//keep here
		beqn++	//keep here
	}

	if (!gml.nogen) {
		build = sum(eltype:==3) | sum(eltype:==8) | sum(eltype:==9) | sum(eltype:==14) | sum(eltype:==15) | sum(eltype:==17)
		if (build) {
			id = st_addvar("double",names)
			st_store(.,id,gml.touse,newvars)
			printf("variables created for model "+strofreal(mod)+", component "+strofreal(c)+": "+stub+"_1 to "+stub+"_"+strofreal(Nvars)+"\n")
		}
	}

	//post for xb equations
	st_local("xb"+strofreal(mod),st_local("xb"+strofreal(mod)) + " "+invtokens(eqnames))
	return(newvars)
}

void merlin_build_dxz(`gml' gml)
{
	gml.dX 		= asarray_create("real",1)
	gml.survind     = -99
	
	for (mod=1;mod<=gml.Nmodels;mod++) {
		
		if (gml.familys[mod]=="rp" | gml.familys[mod]=="logchazard") {

			gml.model 	= mod
			gml.modtoind 	= mod
			N		= gml.N
			Nels		= asarray(gml.Nels,mod)
			X 		= J(N,0,.)
                        Nvars_per_c     = strtoreal(tokens(gml.E.Nvars[mod]))
                        
			for (c=1;c<=gml.Ncmps[mod];c++) {
				
				eltype 		= asarray(gml.elindex,(mod,c))
				istvar 		= asarray(gml.eltvar,(mod,c))[,1]
				anytvars 	= sum(istvar)

				if (anytvars) { 
					
					//if more than one then need to use the chain rule
					if (anytvars>1) {
						
						newvars = J(N,1,0)
						for (el=1;el<=Nels[c];el++) {
						
							newvars2 = J(N,1,1)
							for (el2=1;el2<=Nels[c];el2++) {

								if (el!=el2) {
									if 	(eltype[el2]==1) nextvars = st_data(.,asarray(gml.elvarname,(mod,c,el2)),gml.touse)
									else if (eltype[el2]==8) nextvars = merlin_xz_rcs(gml,c,el2,0)
									else if (eltype[el2]==9) nextvars = merlin_xz_fp(gml,c,el2,0)
								}
								else {
									if (istvar[el]) {
										if 	(eltype[el2]==1) nextvars = J(N,1,1)
										else if (eltype[el2]==8) nextvars = merlin_xz_rcs(gml,c,el2,1)
										else if (eltype[el2]==9) nextvars = merlin_xz_fp(gml,c,el2,1)
									}
									else {
										if 	(eltype[el]==1)	nextvars = J(N,1,0)
										else if (eltype[el]==8) {
											rcsinfo  = asarray(gml.elinfo,(gml.model,c,el))
											nextvars = J(N,cols(asarray(rcsinfo,3))-1,0)
										}
										else if (eltype[el]==9) {
											fpinfo 	 = asarray(gml.elinfo,(gml.model,c,el))
											nextvars = J(N,cols(asarray(fpinfo,3)),0)
										}
									}
								}
								
								//rebuild
								Nnew 		= cols(nextvars)
								copyold 	= newvars2
								newvars2 	= J(N,0,.)
								for (j=1;j<=Nnew;j++) {
									newvars2 = newvars2,(copyold :* nextvars[,j])
								}	
							}
							newvars = newvars :+ newvars2
						}
							
					}
					else {
						newvars = J(N,1,1)
						for (el=1;el<=Nels[c];el++) {
							if (istvar[el]) {
								if 	(eltype[el]==1)	nextvars = J(N,1,1)
								else if (eltype[el]==8) nextvars = merlin_xz_rcs(gml,c,el,1)
								else if (eltype[el]==9) nextvars = merlin_xz_fp(gml,c,el,1)
							}
							else {
								if 	(eltype[el]==1)	nextvars = st_data(.,asarray(gml.elvarname,(mod,c,el)),gml.touse)
								else if (eltype[el]==8) nextvars = merlin_xz_rcs(gml,c,el,0)
								else if (eltype[el]==9) nextvars = merlin_xz_fp(gml,c,el,0)
							}
							//rebuild
							Nnew = cols(nextvars)
							copyold = newvars
							newvars = J(N,0,.)
							for (j=1;j<=Nnew;j++) {
								newvars = newvars,(copyold :* nextvars[,j])
							}
						}	
					}
					X = X,newvars
				}
				else    X = X,J(N,Nvars_per_c[c],0)
                                
			}

			//constant
			if (gml.hascons[mod]) X = X,J(N,1,0)
			//store design matrix
			asarray(gml.dX,mod,X)

		}
	}
	gml.survind = 0
}



void merlin_get_cmp_labels(`gml' gml)
{
	gml.cmplabels = J(gml.Nmodels,1,"")
	
	for (mod=1;mod<=gml.Nmodels;mod++) {
		
		depvars = merlin_get_indepvars(mod)
		
		for (i=1;i<=gml.Ncmps[mod];i++) {
			
			dv 	= strtrim(depvars[i])
			
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
			if (substr(el,1,2)=="EV")	return("EV[]")
			else if (strpos(el,"iEV")) 	return("iEV[]")
			else if (strpos(el,"dEV")) 	return("dEV[]")
			else if (strpos(el,"d2EV")) 	return("d2EV[]")
		}
		else if (strpos(el,"XB")) {	
			if (substr(el,1,2)=="XB")	return("XB[]")
			else if (strpos(el,"iXB")) 	return("iXB[]")
			else if (strpos(el,"dXB")) 	return("dXB[]")
			else if (strpos(el,"d2XB")) 	return("d2XB[]")
		}
		else {
			return(el)	//re
		}
			
	}
	else {	
		
		if 	(strpos(el,"mf("))      return("mf()")	//user-defined function
		else if (strpos(el,"rcs("))     return("rcs()")	//rcs function
		else if (strpos(el,"fp("))      return("fp()")	//fp function
		else if (strpos(el,"bs(")) 	return("bs()")	//bs function
		else if (strpos(el,"pc(")) 	return("pc()")	//pc function
		else 				return(el)	//varname

	}
	
}

end
