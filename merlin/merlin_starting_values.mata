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

void merlin_starting_values(`gml' gml)
{

	from 	= st_local("from")
	hasfrom = from!=""

	if (hasfrom) {
	
		stata("confirm matrix "+from)
		gml.myb = st_matrix(from)
		merlin_xb(gml,gml.myb)
		
		if (st_local("blups")!="") {
			
			i=1
			st_view(newmu_i,.,tokens(st_local("blups")),"_merlin_flag")
			st_view(newtau_i,.,tokens(st_local("seblups")),"_merlin_flag")

			ndim 			= gml.ndim[i]
			logdetChol 		= J(gml.Nobs[i,1],1,.)
			stackednodes 	= J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
			
			ind1 = 1; ind2 = ndim
			for (j=1; j<=gml.Nobs[i,1]; j++) {
				shift 			= newmu_i[j,]'
				vcv_new 		= diag(newtau_i[j,]):^2
				nodes 			= shift :+ vcv_new * asarray(gml.baseGHnodes,i)
				logdetChol[j] 	= log(sqrt(det(vcv_new*vcv_new)))
				asarray(gml.aghip,(i,j),nodes)
				stackednodes[|ind1,.\ind2,.|] = nodes'
				ind1 = ind1 + ndim
				ind2 = ind2 + ndim
			}

			//update logl extra contribution
			asarray(gml.aghlogl,i,rowshape(merlin_lndmvnorm(stackednodes,I(gml.Nres[i])),gml.Nobs[i,1]) :+ logdetChol)
			
			//update stacked nodes, and RE specific stacked nodes
			asarray(gml.stackednodes,i, stackednodes')

			res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)
			for (r=1;r<=gml.Nres[i];r++) {
				asarray(gml.aghip2,(i,r),rowshape(res[r,],gml.Nobs[i,1]))
			}
			
		}
		
	}
	
	if (!gml.predict) {

		//starting values
		if (st_local("random")=="")     binit = J(1,gml.Nb,0)
		else 			        binit = runiform(1,gml.Nb)
		
		if (!hasfrom) {
			//intercept for addhazard
			for (mod=1;mod<=gml.Nmodels;mod++) {
                                if (gml.familys[mod]=="addhazard") binit[gml.eqnindex[mod,2]] = 0.1
                        }

			//distap defaults -> needed for rp first spline term and ordinal cuts
			binit[gml.initdistapindex[1,]] = gml.initdistapindex[2,]
		}

		if (gml.Nrelevels & !hasfrom & st_local("zeros")=="" & st_local("random")=="") {

			cmd = st_local("ZERO")

			//get rid of all components which have a [ in them and put back together again
			
			stata("_parse expand EQ GL : ZERO")
			newcmd = ""
			for (mod=1; mod<=gml.Nmodels; mod++) {
				clp = st_local("EQ_"+strofreal(mod))
				stata("local 0 "+clp)
				stata("syntax anything [if] [in], [*]")
				opts 	= st_local("options")
				clp 	= st_local("anything")
				newclp 	= J(1,0,"")
				st_local("clp",clp)
				stata(`"gettoken cmp clp : clp, parse(" ") bind"')
				if (!strpos(st_local("cmp"),"[")) newclp = newclp,st_local("cmp")
				while (st_local("clp")!="") {
					stata(`"gettoken cmp clp : clp, parse(" ") bind"')
					if (!strpos(st_local("cmp"),"[")) newclp = newclp,st_local("cmp")
				}
				newcmd = newcmd +" ("+invtokens(newclp)+" "+st_local("if")+st_local("in")+","+opts+" )"
			} 

			//devcodes
			dc = " devcode1("+st_local("devcode1")+")"
			dc = dc + " devcode2("+st_local("devcode2")+")"
			dc = dc + " devcode3("+st_local("devcode3")+")"
			dc = dc + " devcode4("+st_local("devcode4")+")"
			dc = dc + " devcode5("+st_local("devcode5")+")"
			dc = dc + " devcode6("+st_local("devcode6")+")"
			dc = dc + " devcode7("+st_local("devcode7")+")"
			
			stata("di ")
			stata(`"di as text "Fitting fixed effects model:""')
			rc = _stata("qui merlin "+ newcmd + st_local("GL_if")+st_local("GL_in") + ", nogen"+dc)
			if (rc) exit(rc)
			
			//update with fixed effects
			binit[gml.initbindex] = st_matrix("e(b)")

		}

		if (hasfrom) {
			gml.gridsearch = 0
			binit = st_matrix(st_local("from"))
		}
		else {
			if (gml.Nrelevels) 	gml.gridsearch = 1
			else			gml.gridsearch = 0
		}

		//pass in ap start values
		if (st_local("apstartvalues")!="") binit[gml.initapindex[1,]] = gml.initapindex[2,]

		//pass in vcv inits, in case restartvalues was specified
		if (st_local("restartvalues")!="") binit[gml.initvarindex[1,]] = gml.initvarindex[2,]

		stata("tempname from")
		st_matrix(st_local("from"),binit)
				
	}
}

void merlin_vcv_gridsearch(`gml' gml, `RR' binit)
{
	if (gml.Nrelevels) {

		reind 		= 1
		gml.myb 	= binit
		gml.survind = 0
		merlin_xb(gml,binit)
		logl1		= sum(merlin_logl_panels(1,gml),1)
		covariances = gml.covariances

		for (i=1;i<=gml.Nrelevels;i++) {
			if (covariances[1,i]) {
				for (j=1;j<=gml.Nres[i];j++) {
					if (!gml.initvarindex[2,reind]) {
						copyb = binit
						binit[gml.initvarindex[1,reind]] = log(sqrt(0.1))
						gml.myb 	= binit
						gml.survind = 0
						merlin_xb(gml,binit)
						logl2		= sum(merlin_logl_panels(1,gml),1)
						if (!(logl2>logl1 & (logl2:-logl1)>5 & !missing(logl2))) 	{
								binit = copyb
						}
						else 	logl1 = logl2
					}
					reind++
				}
			}
			else if (covariances[2,i]) {
				if (!gml.initvarindex[2,reind]) {
					copyb = binit
					binit[gml.initvarindex[1,reind]] = log(sqrt(0.1))
					gml.myb 	= binit
					gml.survind = 0
					merlin_xb(gml,binit)
					logl2		= sum(merlin_logl_panels(1,gml),1)
					if (!(logl2>logl1 & (logl2:-logl1)>5 & !missing(logl2))) 	{
							binit = copyb
					}
					else 	logl1 = logl2
				}
				reind++
				reind++
			}
			else if (covariances[3,i]) {
				for (j=1;j<=gml.Nres[i];j++) {
					if (!gml.initvarindex[2,reind]) {
						copyb = binit
						binit[gml.initvarindex[1,reind]] = log(sqrt(0.1))
						gml.myb 	= binit
						gml.survind = 0
						merlin_xb(gml,binit)
						logl2		= sum(merlin_logl_panels(1,gml),1)
						if (!(logl2>logl1 & (logl2:-logl1)>5 & !missing(logl2))) 	{
								binit = copyb
						}
						else 	logl1 = logl2
					}
					reind++
				}
				for (j=1;j<=gml.Nres[i];j++) {
					ind = 1
					while (ind<j) {
						reind++
						ind++
					}
				}
			}
			else {
				if (!gml.initvarindex[2,reind]) {
					copyb = binit
					binit[gml.initvarindex[1,reind]] = log(sqrt(0.1))
					gml.myb 	= binit
					gml.survind = 0
					merlin_xb(gml,binit)
					logl2		= sum(merlin_logl_panels(1,gml),1)
					if (!(logl2>logl1 & (logl2:-logl1)>5 & !missing(logl2))) 	{
							binit = copyb
					}
					else 	logl1 = logl2
				}
				reind++
			}
		}

	}
	gml.gridsearch = 0
}

end
