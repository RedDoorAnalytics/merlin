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

void merlin_init_integration(`gml' gml, levelvars, modeltouses)
{	
	//setup level specific integration techniques
	
	if (st_local("intmethod")=="") {
		gml.usegh = J(gml.Nrelevels,1,1)
		gml.adapt = J(gml.Nrelevels,1,1)
	}
	else {
		intmethod 	= tokens(st_local("intmethod"))'
		if (rows(intmethod)!=gml.Nrelevels) intmethod = J(gml.Nrelevels,1,intmethod)
		strlen 		= strlen(intmethod)
		gml.usegh	= gml.adapt = J(gml.Nrelevels,1,.)
		for (i=1;i<=gml.Nrelevels;i++) {
			if (intmethod[i]==substr("ghermite",1,max((strlen[i],2)))) {
				gml.usegh[i] = 1
				gml.adapt[i] = 0
			}
			else if (intmethod[i]==substr("mvaghermite",1,max((strlen[i],5)))) {
				gml.usegh[i] = 1
				gml.adapt[i] = 1
			} 
			else if (intmethod[i]==substr("mcarlo",1,max((strlen[i],2)))) {
				gml.usegh[i] = 0
				gml.adapt[i] = 0
			}
			else {
				errprintf("Invalid intmethod()\n")
				exit(198)
			}
		}
		if (rows(intmethod)>1 & rows(intmethod)!=gml.Nrelevels) {
			errprintf("Invalid intmethod()\n")
			exit(198)
		}
	}
	
	//setup integration points
	
	gml.ip = J(gml.Nrelevels,1,0)
	
	if (st_local("intpoints")=="") {							//defaults
		
		for (i=1;i<=gml.Nrelevels;i++) {
			if (gml.usegh[i]) 	gml.ip[i] = 7
			else 				gml.ip[i] = gml.Nres[i] * 50
		}
		
	} 
	else {

		ip = strtoreal(tokens(st_local("intpoints")))'
		if (rows(ip)>1 & rows(ip)!=gml.Nrelevels) {
			errprintf("Invalid intpoints()\n")
			exit(198)
		}
		if (rows(ip)==1) {
			for (i=1;i<=gml.Nrelevels;i++) {
				gml.ip[i] = ip
			}
		}
		else gml.ip = ip

		for (i=1;i<=gml.Nrelevels;i++) {
			if (gml.usegh[i]) {
				if (gml.ip[i]<3& gml.ip[i]>150) {
					errprintf("intpoints() must be between 3 and 150 with Gauss-Hermite quadrature\n")
					exit(198)
				}
			}
			else {
				if (gml.ip[i]<10) {
					errprintf("intpoints() must be >=10 with Monte Carlo integration\n")
					exit(198)
				}
			}
		}
	}
	
	//setup level-specific random effect distributions
	
	if (st_local("redistribution")=="") {
		gml.renormal = J(gml.Nrelevels,1,1)
	}
	else {
		redist = tokens(st_local("redistribution"))'
		if (min((redist:=="normal"):+(redist:=="t"))!=1) {
			errprintf("Invalid redistribution()\n")
			exit(198)
		}
		if (rows(redist)>1 & rows(redist)!=gml.Nrelevels) {
			errprintf("Invalid number of distributions in redistribution()\n")
			exit(198)
		}
		
		if (rows(redist)==1) {
			if (redist=="normal") gml.renormal = J(gml.Nrelevels,1,1)
			else gml.renormal = J(gml.Nrelevels,1,0)
		}
		else gml.renormal = (redist:=="normal")
		
		//update usegh for t dist default
		if (st_local("intmethod")=="") {
			for (i=1;i<=rows(gml.usegh);i++) {
				if (gml.renormal[i]==0) gml.usegh[i] = 0
			}
		}
		
		if (max((redist:=="t") :+ gml.usegh)==2) {
			errprintf("Cannot use intmethod(gh/mvagh) with redistribution(t)\n")
			exit(198)
		}
		
		if (min(gml.renormal)==0 & st_local("df")=="") {
			errprintf("df() must be specified with redistribution(t)\n")
			exit(198)
		} 
		
		df = strtoreal(tokens(st_local("df")))'
		if (rows(df)==1) gml.df = J(gml.Nrelevels,1,df)
		else {
			if (rows(df)!=sum(gml.renormal:==0)) {
				errprintf("Invalid number of df()\n")
				exit(198)			
			}
			gml.df = J(gml.Nrelevels,1,0)
			ind = 1
			for (i=1;i<=gml.Nrelevels;i++) {
				if (gml.renormal[i]==0) {
					gml.df[i] = df[ind]
					ind++
				}
			}
		}		
	}
	
	//core index that gets updated through the recursive function
	
	gml.qind = 1,J(1,gml.Nrelevels,0)	
	
	st_local("mlprolog","derivprolog(merlin_prolog())")
	
	gml.iter = 0
	if (sum(gml.adapt) | gml.todo) gml.Li_ip = asarray_create("real",cols(gml.qind))
	
	if (max(gml.usegh)==1) {
		if (max(gml.adapt :* gml.usegh)==1) {
			gml.iter = 0
			gml.atol = 1e-8
			gml.Pupdateip = &merlin_gh_update_ip()
			gml.showadapt = st_local("showadapt")!=""
		}
	}
	if (min(gml.usegh)==0) {	
		gml.seed = rseed()
		gml.bdraws = asarray_create("real",1)
		gml.atol = 1e-8
		if (max((gml.usegh:==0) :* gml.adapt)==1) {
			gml.Pupdateip = &merlin_mc_update_ip()
			gml.showadapt = st_local("showadapt")!=""
		}
	}
	
	//setup NI
	merlin_init_ip(gml, levelvars, modeltouses)
}

void merlin_init_ip(`gml' gml, levelvars, modeltouses)
{
	
	//initialise stuff needed if at least one level's integ is adaptive
	if (sum(gml.adapt)) {
		gml.aghip 				= asarray_create("real",2)
		gml.aghip2 				= asarray_create("real",2)
		gml.aghlogl				= asarray_create("real",1)	
		gml.adpanelindexes 		= asarray_create("real",2)
		gml.adNmeas 			= asarray_create("real",2)
		gml.stackednodes 		= asarray_create("real",1)
		for (j=1;j<=gml.Nmodels;j++) {
			ids = st_data(.,levelvars[1,],modeltouses[1,j])
			for (i=1;i<gml.Nlevels;i++) {
				if (gml.adapt[i]) {
					asarray(gml.adpanelindexes,(i,j),ids[,i])
				}
			}
		}
	}
	else {
		gml.b = asarray_create("real",1)
// 		if (sum(gml.hasRExb)) {
// 			gml.aghip2 = asarray_create("real",2)
// 			gml.adpanelindexes = asarray_create("real",2)
// 			for (j=1;j<=gml.Nmodels;j++) {
// 				ids = st_data(.,levelvars[1,],modeltouses[1,j])
// 				for (i=1;i<gml.Nlevels;i++) {
// 					asarray(gml.adpanelindexes,(i,j),ids[,i])
// 				}
// 			}
// 		} 
		
	}
	
	//initialise stuff for GHQ
	if (sum(gml.usegh)) {
		gml.baseGHnodes = gml.baseGHweights = asarray_create("real",1)
	}
	
	gml.ndim = J(gml.Nrelevels,1,0)

	for (i=1;i<gml.Nlevels;i++) {
		
		if (gml.usegh[i]) {
						
			if (gml.Nres[i]) {
			
				qmat 		= _gauss_hermite_nodes(gml.ip[i])
				gml.ndim[i] = gml.ip[i]:^gml.Nres[i]	
				
				x 			= J(gml.Nres[i],1,qmat[1,]):*sqrt(2)
				asarray(gml.baseGHnodes,i,merlin_expand_matrix(x))
				
				w	 		= J(gml.Nres[i],1,qmat[2,]):/sqrt(pi())
				asarray(gml.baseGHweights,i,merlin_expand_matrix(w,1)')
				
				if (gml.adapt[i]) {
					
					if (st_local("blups")=="") {
					
						asarray(gml.stackednodes,i,J(1,gml.Nobs[i,1],asarray(gml.vcvs,i) * asarray(gml.baseGHnodes,i)))
						for (r=1;r<=gml.Nres[i];r++) {
							asarray(gml.aghip2,(i,r),rowshape(asarray(gml.stackednodes,i)[r,],gml.Nobs[i,1]))
						}	
						for (j=1; j<=gml.Nobs[i,1]; j++) {
							asarray(gml.aghip,(i,j),asarray(gml.vcvs,i) * asarray(gml.baseGHnodes,i))
						}
						asarray(gml.baseGHweights,i,((2:*pi()):^(gml.Nres[i]:/2):*exp(quadcolsum(asarray(gml.baseGHnodes,i):^2):/2) :* asarray(gml.baseGHweights,i)')')
						asarray(gml.aghlogl,i,J(gml.Nobs[i,1],1,merlin_lndmvnorm(asarray(gml.baseGHnodes,i)',I(gml.Nres[i]))') :+ log(sqrt(det(asarray(gml.vcvs,i)))))
					
					}
					else {
						
// 						id = tokens(st_local("idvars"))[1,i]
// 						id
// 						var = "_tempid"+strofreal(i)
					
// 						st_view(newmu_i,.,tokens(st_local("blups")),st_local(var))
// 						st_view(newtau_i,.,tokens(st_local("seblups")),st_local(var))

// 						ndim = gml.ndim[i]
// 						logdetChol = J(gml.Nobs[i,1],1,.)
// 						stackednodes = J(ndim*gml.Nobs[i,1],gml.Nres[i],.)
						
// 						ind1 = 1; ind2 = ndim
// 						for (j=1; j<=gml.Nobs[i,1]; j++) {
// 							shift = newmu_i[j,]'
// 							vcv_new = diag(newtau_i[j,])
// 							nodes = shift :+ vcv_new * asarray(gml.baseGHnodes,i)
// 							logdetChol[i] = log(sqrt(det(vcv_new*vcv_new)))
// 							asarray(gml.aghip,(i,j),nodes)
// 							stackednodes[|ind1,.\ind2,.|] = nodes'
// 							ind1 = ind1 + ndim
// 							ind2 = ind2 + ndim
// 						}

// 						//update logl extra contribution
// 						asarray(gml.aghlogl,i,rowshape(merlin_lndmvnorm(stackednodes,I(gml.Nres[i])),gml.Nobs[i,1]) :+ logdetChol)
						
// 						//update stacked nodes, and RE specific stacked nodes
// 						asarray(gml.stackednodes,i, stackednodes')

// 						res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)
// 						for (r=1;r<=gml.Nres[i];r++) {
// 							asarray(gml.aghip2,(i,r),rowshape(res[r,],gml.Nobs[i,1]))
// 						}
						
					}
				}
				else {
					asarray(gml.b,i,asarray(gml.baseGHnodes,i))				
				}
			}
		}
		else {
			gml.ndim[i] = gml.ip[i]
			if (gml.adapt[i]) {
				/*
				asarray(gml.stackednodes,i,J(1,gml.Npanels[i],asarray(gml.baseGHnodes,i)))
				for (r=1;r<=gml.Nres[i];r++) {
					asarray(gml.aghip2,(i,r),rowshape(asarray(gml.stackednodes,i)[r,],gml.Npanels[i]))
				}	
				for (j=1; j<=gml.Npanels[i]; j++) {
					asarray(gml.aghip,(i,j),asarray(gml.baseGHnodes,i))
				}
				asarray(gml.baseGHweights,i,((2:*pi()):^(gml.Nres[i]:/2):*exp(quadcolsum(asarray(gml.baseGHnodes,i):^2):/2) :* asarray(gml.baseGHweights,i)')')
				asarray(gml.aghlogl,i,J(gml.Npanels[i,1],1,merlin_lndmvnorm(asarray(gml.baseGHnodes,i)',I(gml.Nres[i]))'))
				*/
			
				rseed(gml.seed)			//reset seed
				baseip = invnormal(halton(gml.ndim[i],gml.Nres[i])')
				asarray(gml.bdraws,i,baseip)
				
				
				asarray(gml.stackednodes,i,J(1,gml.Nobs[i,1],baseip))
				for (r=1;r<=gml.Nres[i];r++) {
					asarray(gml.aghip2,(i,r),rowshape(asarray(gml.stackednodes,i)[r,],gml.Nobs[i,1]))
				}
				for (j=1; j<=gml.Nobs[i,]; j++) {
					asarray(gml.aghip,(i,j),baseip)
				}
				asarray(gml.aghlogl,i,J(gml.Nobs[i,1],gml.ndim[i],0))
				
				/*
				gml.ndim = J(gml.Nlevels,1,gml.ip)
				for (i=1; i<=gml.Nlevels; i++) {
					if (gml.Nres[i]) {
						for (j=1; j<=gml.Npanels[i,1]; j++) {
							if (gml.renormal[i]) asarray(gml.aghip,(i,j),merlin_drawnorm(J(gml.Nres[i],1,0),asarray(gml.vcvs,i),gml.ndim[i]))
							else  asarray(gml.aghip,(i,j),merlin_drawt(J(gml.Nres[i],1,0),asarray(gml.vcvs,i),gml.df[i],gml.ndim[i]))
							asarray(gml.aghip2,(i,j),asarray(gml.aghip,(i,j)))
						}
						asarray(gml.aghlogl,i,J(gml.Npanels[i,1],gml.ndim[i],0))
					}
				}*/
			}
			else {
				rseed(gml.seed)			//reset seed
				asarray(gml.bdraws,i,invnormal(halton(gml.ndim[i],gml.Nres[i])'))
				//asarray(gml.bdraws,i,invnormal(runiform(gml.Nres[i],gml.ndim[i])))
				/*draw = runiform(gml.Nres[i],gml.ndim[i]:/2)
				draw = draw,(1:-draw)
				asarray(gml.bdraws,i,invnormal(draw))*/
				
				//if (gml.renormal[i]) 	bnodes = merlin_drawnorm(J(gml.Nres[i],1,0),asarray(gml.vcvs,i),gml.ndim[i])
				//else 					bnodes = merlin_drawt(J(gml.Nres[i],1,0),asarray(gml.vcvs,i),gml.df[i],gml.ndim[i])
				asarray(gml.b,i,asarray(gml.bdraws,i))
			}
		
		}
	}

	//note; Npanels is indexed by level and model -> they are the same across all models (until ob level which isn't used here)
}

end
