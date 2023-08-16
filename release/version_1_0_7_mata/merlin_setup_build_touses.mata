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
	allvars			= J(1,0,"")

	for (i=1;i<=gml.Nmodels;i++) {
		
		vars 	= merlin_get_indepvars(i)
		Nvars 	= cols(vars)

		modelvars 	= idvars
		for (j=1;j<=Nvars;j++) {
			
			if (strpos(vars[j],"@")) {
				vars[j] = substr(vars[j],1,strpos(vars[j],"@")-1)
			}
		
			dv 	= strtrim(vars[j])
			if (!strpos(dv,"[") & !strpos(dv,"(")) {
				modelvars = modelvars,dv											//one element, and its a variable
			}
			pos = 1	
			while (pos) {
				pos = strpos(dv,"#")
				if (pos) {
					dv2 = substr(dv,1,pos-1)
					dv 	= substr(dv,pos+1,.)
				}
				else dv2 = dv
				if (!strpos(dv2,"[") & !strpos(dv2,"(")) {							//variable name
					modelvars = modelvars,dv2								
				}
				//rcs()/fp()/bs()/ps()/gp() varname
				if (strpos(dv2,"fp(") | strpos(dv2,"rcs(") | strpos(dv2,"bs(") | strpos(dv2,"ps(") | strpos(dv2,"gp(")) {	
					hasEV = strpos(dv2,"EV[")				//avoid fn(EV[],...)
					if (!hasEV) {		
						p1 	= strpos(dv2,"(")+1
						len = strpos(dv2,",") - p1
						modelvars = modelvars,strtrim(substr(dv2,p1,len))
					}
				}
			}
		}

		//post into Stata fixed variable for mleqn
		//st_local("xb"+strofreal(i),invtokens(modelvars))
		
		//add timevar and response vars
		//->need response var here as likelihood contribution is only from observed responses
		
		//can override predicts timevar
		if (!gml.predict) modelvars = modelvars,st_local("timevar"+strofreal(i))					

		//if survival model, can override predicts timevar
		predandsurv = gml.predict & st_local("failure"+strofreal(i))!=""
		if (!predandsurv) {
			modelvars = modelvars,tokens(st_local("response"+strofreal(i)))
		}

		//post model specific touse var
		modeltousei = "touse"+strofreal(i)
		stata("tempname "+modeltousei)
		stata("qui gen byte "+st_local(modeltousei)+"="+gml.touse+"==1")
		stata("markout "+st_local(modeltousei)+" "+invtokens(modelvars))

		//if being called from predict(ms), update on N
		if (gml.predict) {
			if (st_local("ptvar")!="") {
				stata("markout "+st_local(modeltousei)+" "+st_local("ptvar"))
			}
			if (st_local("npredict")!="") {
				stata("qui replace "+st_local(modeltousei)+" = _n<="+st_local("npredict"))
			}
		}

		//update global touse based on modeltouses -> finished below
		if (i==gml.Nmodels) st_local("touseup",st_local("touseup")+st_local(modeltousei)+"==1 ")
		else st_local("touseup",st_local("touseup")+st_local(modeltousei)+"==1 | ")
		
		gml.modeltouses = gml.modeltouses,st_local(modeltousei)
		allvars = allvars,modelvars
	}

	//to post for predict(ms)
	gml.allvars = invtokens(uniqrows(allvars')')
	
	//update global touse based on modeltouses
	stata("qui replace "+gml.touse+" = ("+st_local("touseup")+")")
}

end
