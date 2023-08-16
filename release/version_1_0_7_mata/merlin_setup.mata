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

struct merlin_ereturn_struct {
	`SR' levelvars				//level vars 
	`RC' Nreparams				//# of VCV parameters at each level
	`SC' reeqns					//eqn names for re parameters
	`SC' reivscale				//inverse scale of re parameters
	`SC' relabel				//labels for re parameters
	`SM' userfunctions			
	`SC' Nvars
	`RS' merlin
}

struct merlin_struct {

// [op] - denotes an option which is not always completely filled in

	`pgml' Pgml					//self-pointer
	
	//model definition
	`TR' y						//response variables
	`SC' familys				//family of each model
	`SC' responses				//response var names 							-> includes event and ltruncated if there
	`SC' links					//link function of each model
	`PC' invlinks				//inverse link function of each model
	`SC' failures		//[op]	//failure variable name of each model 			-> for eret list
	`SC' ltruncated		//[op]	//ltruncated variable name of each model		-> for eret list
	`RC' hasltrunc				//has left truncation/delayed entry
	`RS' hasanyltrunc			//flag for any left truncation/delayed entry
	`PC' ltruncP
	`SS' touse					//global touse variable
	`SR' modeltouses
	`SR' levelvars				//level varables from highest to lowest
	`RM' NotNull				//flag whether the model contributes to the log likelihood
	`RC' nobhaz					//presence of bhazard for each model
	`TR' bhazards				//stored bhazard variable for each model
	`RM' hastvars				//flag - timevar() specified
	`TR' timevars				//stored timevar()
	
	`RS' Nmodels				//# of models
	`RS' Nlevels				//# of levels (includes ob. level)
	`RS' Nrelevels				//# of levels with random effects
	`RS' N						//# of observations (ob. level, main touse)
	`RC' Nobs					//# of observations for each model = level by models
	`TR' xbindex				//ob. level index for each outcome, based on main touse core index
	
	`RC' Nap			//[op]	//# of anciliary parameters in a user model
	`TR' apxb
	`TR' distancb
	`TR' OrdIndexes		//[op]	//resp index for cuts for ordinal model
	`RC' qtiles			//[op]	//quantiles, if specified
	
	`PC' Plnl, Pxb				//pointers to logL and linear predictor functions for each model
	`RS' todo					//gf0,1,2
	`PC' userhaz		//[op]	//pointer to user-specified hazard functions
	`PC' userchaz		//[op]	//pointer to user-specified cumulative hazard functions
	`PC' userloghaz		//[op]	//pointer to user-specified log hazard functions
	
	//complex syntax
	`RS' nogen					//do not generate Stata variables for elements
	`RR' myb					//copy ml's b vector so I can update it with my numeric scores
	`RC' Ncmps					//number of components in each complex syntax
	`TR' elindex				//index of each element within each component, for each model
	`TR' Nels
	`RM' eqnindex
	`RC' hascons
	`TR' elinfo
	`RC' Ndistancp
	`PC' expvalP
	`RS' Nb						//total number of parameters in model
	`RR' initbindex				//starting value indexes to post fixed effects model results into b for full model fit
	`RM' initdistapindex
	`RM' initapindex			//starting value indexes to post ap into b
	`RM' initvarindex			//index for random effect variance parameters
	`SC' cmplabels
	`TR' CmpBIndex
	
	`TR' dapmodels				//colvector of model indices to call - 0 if no model. Indexed by model.
	
	//integration and random effects
	`RM' covariances
	`TR' panelindexes
	`TR' vcvs
	`TR' baseGHnodes
	`TR' baseGHweights
	`RC' ndim
	`RS' adapt
	`RS' iter
	`RS' atol
	`RS' showadapt
	`RS' df
	`TR' b, bdraws
	`TR' Li_ip
	`TR' aghip
	`TR' aghip2
	`TR' aghlogl
	`TR' stackednodes
	`TR' adpanelindexes
	`PS' Pupdateip
	`SS' seed
	`RC' usegh
	`RC' ip
	`RC' renormal
	`RC' Nres
	`RC' NI
	`TR' haznodes
	`TR' hazweights
	`RC' hazNnodes
	`TR' latlevs				//unique REs at each level for eret list
	`RM' lev_mod_latind
	
	//delayed entry
	`RS' ltflag
	`TR' Li_ip_lt
	`TR' aghip_lt
	`TR' aghip2_lt
	`TR' aghlogl_lt
	`TR' stackednodes_lt
	`TR' b_lt
	
	//gets updated
	`RS' lnfi1
	`RS' model
	`RS' modtoind
	`RR' qind
	
	//for predictions
	`RS' predict				
	`RS' fixedonly				//predictions based only on fixed effects - so skip res in utils
	
	//ereturn list
	`Egml' E
	
	//====================================================================================================//
	//development
	
	`RC' hasRExb
	`RS' updateMC				
	
	//sample weights, modelled weights
	`RS' hasweights
	`RS' Nmodels2
	`TR' weights
	
	//imputation models
	`RM' IsImputed				//flag colvector - model is for imputation
	`TR' ImputeIndex			//stores index for observed and expected responses
	
	//for predictms
	`SS' allvars
	
	//overall penalty
	`RC' haspenalty				//first element is flag, second identifies type
	
	//penalisation for ps elements
	`RS' Npenals
	`RC' penals
	`RM' H
	
	//score
	`TR' Score
	
}

void merlin_setup(`SS' GML,`SS' touse)
{
	`pgml' 	pGML
	`gml' 	gml
	`RS' 	debug

	//get started
	debug 			= st_local("debug")!=""
	if (debug) 		printf("-> Starting merlin_setup()\n")
	pGML 			= crexternal(GML)
	gml.Nmodels		= strtoreal(st_local("neq"))
	gml.N 			= strtoreal(st_local("Nobs"))				
	gml.touse 		= touse
	gml.Pgml 		= pGML
	gml.model 		= 1															//gets updated
	gml.modtoind 	= 1															//gets updated
	gml.nogen		= st_local("nogen")!=""
	gml.Npenals		= 0
	
	
	if (st_local("arthur")!="") 		gml.E.merlin	= 3
	else if (st_local("excalibur")!="") gml.E.merlin	= 2
	else 								gml.E.merlin	= 1
	
	//=======================================================================================================================//
	// refilling struct from fitted merlin model - called by predict or predictms

	gml.predict 		= st_local("predict")!=""
	gml.fixedonly		= 0
	if (gml.predict & st_local("npredict")!="") gml.N			= strtoreal(st_local("npredict"))	
	
	//=======================================================================================================================//
	// models info
	
	gml.familys 	= tokens(st_local("familylist"))'
	gml.links 		= tokens(st_local("linklist"))'								//defunct (mainly) -> to bring back
	gml.NotNull		= J(1,2,((gml.familys:!="null") :* (gml.familys:!="re")))	//contributes to logl, has main xb
	
	//survival outcome stuff
	gml.failures	= J(1,gml.Nmodels,"")
	gml.ltruncated	= J(1,gml.Nmodels,"")
	gml.hasltrunc	= J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		gml.failures[1,i] 	= st_local("failure"+strofreal(i))
		gml.ltruncated[1,i] = st_local("ltruncated"+strofreal(i))
		gml.hasltrunc[i] 	= gml.ltruncated[1,i]!=""
	}
	gml.hasanyltrunc 	= sum(gml.hasltrunc)
	gml.ltflag			= 0														//gets updated when needed
	
	//error checks
	merlin_error_checks(gml)
	gml.hastvars 	= J(gml.Nmodels,2,0)
	merlin_check_clp(gml)														//hastvars gets filled up - models x fe,re
	
		if (debug) printf("error checks done\n")
	
	merlin_get_cluster_varnames(gml)
	gml.E.levelvars = gml.levelvars[1,]											//to post
	gml.Nrelevels	= cols(gml.levelvars)
	gml.Nlevels 	= gml.Nrelevels + 1											//includes ob level
	gml.Nobs		= J(gml.Nlevels,gml.Nmodels,.)
			
		if (debug) printf("cluster variables extracted\n")	
	
	merlin_build_touses(gml)
	
	//now touses are done, can start to read in everything
	
	merlin_get_levelvars(gml)
	merlin_get_ys(gml)	
	merlin_get_timevars(gml)
	
		if (debug) printf("touses, responses, and model indices done\n")
	
	levelvars = gml.levelvars,(st_local("coreindex")\st_local("coreindex"))

	//now for each level get unique M#
	merlin_parse_unique_latents1(gml,levelvars)
	latslevmod 	= merlin_parse_unique_latents2(gml,levelvars)

		if (debug) printf("latents parsed\n")
	
	//initialise VCV matrices
	merlin_init_vcvs(gml)

	//family parameters setup
	gml.distancb 	= asarray_create("real",2)
	gml.Ndistancp 	= J(gml.Nmodels,1,0)
	merlin_rcs_setup(gml)
	merlin_ordinal_setup(gml)
	merlin_gp_setup(gml)
	merlin_get_n_dist_ancp(gml)					//declares dist anc p ml eqns in locals	
	
	//build components and post any new variables for ml equations and index equations
	gml.E.Nvars = J(gml.Nmodels,1,"")
	merlin_noconstants(gml)
	merlin_ancp_init(gml)
	gml.CmpBIndex = asarray_create("real",1)
	merlin_build_cmps(gml)													
	merlin_setup_EV_p(gml)	
	merlin_get_cmps_labels(gml)
	
	if (gml.Npenals) gml.penals = J(gml.Npenals,1,1)		//!!fixed at 1
	
		if (debug) printf("components built\n")	
	
	merlin_penalties(gml)
	
	//setup appropriate logl pointers xbs 
	merlin_get_xb_p(gml)
	merlin_get_logl_p(gml)

		if (debug) printf("pointers done\n")	
	
	//covariance structures
	merlin_parse_covstructures(gml)
	merlin_parse_vcv_eqns(gml)

		if (debug) printf("VCV structures and equations done\n")
	
	//id vars setup
	gml.panelindexes = asarray_create("real",2)

	if (gml.Nlevels>1) {
		for (j=1;j<=gml.Nmodels;j++) {
			ids = st_data(.,levelvars[1,],gml.modeltouses[1,j])
			for (i=1;i<gml.Nlevels;i++) {
				psetup 				= panelsetup(uniqrows(ids[,1..(i+1)]),i)
				gml.Nobs[i,j]	 	= rows(psetup)
				asarray(gml.panelindexes,(i,j),psetup)
			}
		}
	}
	
		if (debug) printf("panel arrays done\n")
	
	//initialise hazard integation
	gml.updateMC = st_local("updatemc")!=""
	merlin_init_numint(gml)
	merlin_hazNI_init(gml)
	
	//initialise random effects integation
	gml.iter = 0
	if (gml.Nlevels>1) {
		merlin_init_integration(gml, levelvars, gml.modeltouses)
		if (gml.hasanyltrunc) merlin_init_ip_ltrunc(gml, levelvars, gml.modeltouses)
	}

		if (debug) printf("numerical integration initialised\n")
	
	merlin_get_weights(gml,gml.modeltouses,levelvars)
	
// 		if (debug) printf("weights done\n")

	gml.todo = strtoreal(substr(st_local("evaltype"),3,1))
	
	//extra stuff needed for particular models
	
	merlin_get_quantiles(gml)
	merlin_get_bhazards(gml)
	
	//ml equations 
		
	space = " "
	for (i=1; i<=gml.Nmodels; i++) {
		mlspeci	= st_local("xb"+strofreal(i))
		mlspeci	= mlspeci+space+st_local("dap_eqns"+strofreal(i))
		mlspeci	= mlspeci+space+st_local("ap_eqns"+strofreal(i))
		st_local("mlspec",st_local("mlspec")+space+mlspeci)
	}
	st_local("mlspec", st_local("mlspec")+space+st_local("vcveqns"))
	
		if (debug) printf("ml equations done\n")
	
	
	//=======================================================================================================================//

		//devcode check
		devcode1 = strtoreal(st_local("devcode1"))	//weights
		devcode2 = strtoreal(st_local("devcode2"))	//called from predictms
		devcode3 = strtoreal(st_local("devcode3"))	//family(prp)
		
		if (sum(gml.hasweights) & devcode1!=242566) merlin_error("Nope")
		
		if (st_local("predictms")!="" & devcode2!=678990) merlin_error("Nope")
		
		if (sum(gml.familys:=="prp") & devcode3!=144930) merlin_error("Nope")
	
	//=======================================================================================================================//	
	
	//starting values
	if (st_local("from")!="") {
		stata("confirm matrix "+st_local("from"))
		gml.myb = st_matrix(st_local("from"))
		merlin_xb(gml,st_matrix(st_local("from")))
	}
	
	if (!gml.predict) {
		merlin_starting_values(gml,latslevmod)
	}
	
		if (debug) printf("starting values done\n")	
	
	//=======================================================================================================================//
	
	//score
	if (gml.todo) gml.Score = asarray_create("real",2) 							//model, complex predictor component
	
	//=======================================================================================================================//

	swap(*pGML,gml)	

		if (debug) printf("-> Finished merlin_setup()\n")	
	
	//Done
	st_local("mlprolog","derivprolog(merlin_prolog())")
	
	//gml.hasRExb 	= J(Nidlevels,1,0)											//dev. work
	//merlin_get_impute_flag(gml)												//dev. work
	
}

`SR' merlin_get_indepvars(`RS' i)
{
	depvars = J(1,0,"")
	st_local("rest",st_local("indepvars"+strofreal(i)))
	stata("gettoken lhs rest : rest, bind")						
	depvars = depvars,st_local("lhs")
	while (st_local("rest")!="") {
		stata("gettoken lhs rest : rest, bind")
		depvars = depvars,st_local("lhs")
	}
	return(depvars)
}

void merlin_get_ys(`gml' gml)
{
	gml.y = asarray_create("real",1)
	gml.responses = J(1,0,"")
	for (i=1;i<=gml.Nmodels;i++) {
		resp = st_local("response"+strofreal(i))
		if (resp!="") {
			gml.responses = gml.responses,resp
			asarray(gml.y,i,st_data(.,tokens(resp),gml.modeltouses[1,i]))
		}
	}
}

void merlin_get_impute_flag(`gml' gml)
{
	gml.IsImputed = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) gml.IsImputed[i] = st_local("impute"+strofreal(i))!=""
	if (sum(gml.IsImputed)) {
		gml.ImputeIndex = asarray_create("real",2)
	}
}

void merlin_get_quantiles(`gml' gml)
{
	gml.qtiles = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.familys=="lquantile") {
			if (st_local("qtile"+strofreal(i))!="") {
				gml.qtiles[i] = strtoreal(st_local("qtile"+strofreal(i)))
				if (gml.qtiles[i]<0 | gml.qtiles[i]>1) {
					merlin_error("Invalid quantile")
				}
			}
			else gml.qtiles[i] = 0.5
		}
	}
}

void merlin_get_timevars(`gml' gml)
{
	gml.timevars = asarray_create("real",1)								//must be declared - used in predict

	hastvars = rowsum(gml.hastvars)
	for (i=1;i<=gml.Nmodels;i++) {
		if (hastvars[i] & st_local("timevar"+strofreal(i))=="") {
			errprintf("timevar() required\n")
			exit(198)
		}
	}
	for (i=1;i<=gml.Nmodels;i++) {
		if (st_local("timevar"+strofreal(i))!="") {
			//indexed by main touse in case one model calls another - > gets indexed
			asarray(gml.timevars,i,st_data(.,st_local("timevar"+strofreal(i)),gml.touse))
		}
	}
}

void merlin_parse_unique_latents1(`gml' gml, `SM' levelvars)
{
	//notes; not effected by @'s
	//	     not affected by EV[]

	//unique latents at each level
	gml.latlevs = asarray_create("real",1)		
	gml.Nres = J(gml.Nlevels,1,0)
	
	for (i=1;i<=gml.Nlevels;i++) {
		lats = J(0,1,"")
		//go through everything for each level
		for (j=1;j<=gml.Nmodels;j++) {		
			depvars = merlin_get_indepvars(j)
			Ndv 	= cols(depvars)
			for (k=1;k<=Ndv;k++) {
				dv = strtrim(depvars[1,k])
				pos = 1	
				while (pos) {
					pos = strpos(dv,"#")
					if (pos) {
						dv2 = substr(dv,1,pos-1)
						dv = substr(dv,pos+1,.)
					}
					else dv2 = dv
					posid = strpos(dv2,levelvars[2,i])
					if (posid) lats = lats\substr(dv2,1,posid-1)	
				}	
			}
		}
		lats 		= uniqrows(lats)			//they could appear in multiple models and multiple times in same model
		gml.Nres[i] = rows(uniqrows(lats))
		asarray(gml.latlevs,i,lats)
	}
}

`TR' merlin_parse_unique_latents2(`gml' gml, `SM' levelvars)
{
	//not effected by @'s
	latslevmod = asarray_create("real",2)					//all latents at each level in each model
	gml.lev_mod_latind = J(gml.Nlevels,gml.Nmodels,.) 		//indicator; any REs at each level in each model
	for (i=1;i<=gml.Nlevels;i++) {
		//go through everything for each level
		for (j=1;j<=gml.Nmodels;j++) {
			asarray(latslevmod,(i,j),J(0,1,""))
			depvars = merlin_get_indepvars(j)
			Ndv 	= cols(depvars)
			for (k=1;k<=Ndv;k++) {
				dv = strtrim(depvars[1,k])
				pos = strpos(dv,levelvars[2,i])
				if (pos) {
					asarray(latslevmod,(i,j),(asarray(latslevmod,(i,j))\substr(dv,1,pos-1)))
				}
			}
			gml.lev_mod_latind[i,j] = (asarray(latslevmod,(i,j))!=J(0,1,""))
		}
	}
	return(latslevmod)
}

void merlin_init_vcvs(`gml' gml)
{
	gml.vcvs = asarray_create("real",1)
	for (i=1;i<=gml.Nlevels;i++) asarray(gml.vcvs,i,I(gml.Nres[i]))
}



void merlin_get_bhazards(`gml' gml)
{
	gml.nobhaz = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) gml.nobhaz[i] = st_local("bhaz"+strofreal(i))==""
	if (!sum(gml.nobhaz)) {
		gml.bhazards = asarray_create("real",1)
		for (i=1;i<=gml.Nmodels;i++) {
			if (!gml.nobhaz[i]) {
				asarray(gml.bhazards,i,st_data(.,st_local("bhaz"+strofreal(i)),gml.modeltouses[1,i]))
			}
		}
	}
}

void merlin_parse_covstructures(`gml' gml)
{
	`SR' covs
	
	covs = tokens(st_local("covariance2"))
	if (covs==J(1,0,"")) {
		//default independent
		gml.covariances = J(1,gml.Nlevels,(1\0\0))
	}
	else {
		//check not more than Nlevels
		ncovs = cols(covs)
		if (ncovs>1 & ncovs!=gml.Nlevels) {
			errprintf("Error in covariance()\n")		
			exit(198)
		}
		if (ncovs==1) 	gml.covariances = J(1,gml.Nlevels,(covs:=="diagonal"\covs:=="exchangeable"\covs:=="unstructured"))
		else 			gml.covariances = covs:=="diagonal"\covs:=="exchangeable"\covs:=="unstructured"
	}
	stata("tempname vcvmat")
	st_matrix(st_local("vcvmat"),gml.covariances)
}

void merlin_parse_vcv_eqns(`gml' gml)
{
	`SR' vcveqns
	`RR' vcvstartvals
	`SS' vars
	
	vcveqns = J(1,0,"")
	gml.E.Nreparams = J(gml.Nlevels,1,0)
	gml.E.reeqns = J(gml.Nlevels,1,"")
	gml.E.reivscale = J(gml.Nlevels,1,"")
	gml.E.relabel = J(gml.Nlevels,1,"")
	hasresvs = st_local("restartvalues")!=""
	
	gml.initvarindex = J(2,0,0)
	eqn = gml.Nb + 1
	
	if (hasresvs) {
		sv = tokens(st_local("restartvalues"))
		nsv = cols(sv)
		if (mod(nsv,2)) {
			merlin_error("Error in restartvalues()")
		}
		nsv = nsv/2
		ind1 = 1
		ind2 = 2
		re = J(1,0,"")
		val = J(1,0,.)
		for (i=1;i<=nsv;i++) {
			re = re,sv[1,ind1]
			val = val,strtoreal(sv[1,ind2])
			ind1 = ind1+2
			ind2 = ind2+2
		}
	}

	for (lev=1;lev<=gml.Nlevels;lev++) {
		lats = asarray(gml.latlevs,lev)	
		strlev = strofreal(lev)
		//ind or unstr
		if (gml.covariances[1,lev] | gml.covariances[3,lev]) {
			for (r=1;r<=gml.Nres[lev];r++) {
				gml.initvarindex = gml.initvarindex,(eqn++\0)
				todo = 1
				vars = merlin_get_vcv_vars(gml, lev)
				vcveqns = vcveqns,"(lns"+strlev+"_"+strofreal(r)+": "+vars+")"
				gml.E.Nreparams[lev] 	= gml.E.Nreparams[lev] + 1
				gml.E.reeqns[lev] 	= gml.E.reeqns[lev] + "lns"+strlev+"_"+strofreal(r)+" "
				gml.E.reivscale[lev] = gml.E.reivscale[lev] + "exp "
				gml.E.relabel[lev] 	= gml.E.relabel[lev] + "sd("+lats[r]+") " 
				if (hasresvs) {
					for (isv=1;isv<=nsv;isv++) {
						if (lats[r]==re[1,isv] & todo==1) {
							gml.initvarindex[2,cols(gml.initvarindex)] = log(sqrt(val[1,isv]))
							todo = 0
						}
					}
				}
			}
		}
		else {
			gml.initvarindex = gml.initvarindex,eqn++
			vars = merlin_get_vcv_vars(gml, lev)
			vcveqns = vcveqns,"(lns"+strlev+"_1: "+vars+")"
			gml.E.reeqns[lev] = gml.E.reeqns[lev] + "lns"+strlev+"_1 "
			gml.E.reivscale[lev] = gml.E.reivscale[lev] + "exp "
			gml.E.relabel[lev] 	= gml.E.relabel[lev] + "sd " 
			gml.E.Nreparams[lev] = gml.E.Nreparams[lev] + 1
			if (hasresvs) {
				for (sv=1;sv<=nsv;sv++) {
					if (lats[r]==re[1,sv]) {
						gml.initvarindex[2,cols(gml.initvarindex)] = log(sqrt(val[1,isv]))
					}
				}
			}
		}
		
		//exch and nres>1
		if (gml.covariances[2,lev] & gml.Nres[lev]>1) {
			eqn++
			vcveqns = vcveqns,"(art"+strlev+"_1_1: )"
			gml.E.reeqns[lev] = gml.E.reeqns[lev] + "art"+strlev+"_1_1 "
			gml.E.reivscale[lev] = gml.E.reivscale[lev] + "tanh "
			gml.E.relabel[lev] 	= gml.E.relabel[lev] + "corr " 
			gml.E.Nreparams[lev] = gml.E.Nreparams[lev] + 1
		}
		else if (gml.covariances[3,lev]) {
			ind = 1
			while (ind<gml.Nres[lev]) {
				for (r=ind+1;r<=gml.Nres[lev];r++) {
					eqn++
					vcveqns = vcveqns,"(art"+strlev+"_"+strofreal(ind)+"_"+strofreal(r)+": )"
					gml.E.reeqns[lev] = gml.E.reeqns[lev] + "art"+strlev+"_"+strofreal(ind)+"_"+strofreal(r)+" "
					gml.E.reivscale[lev] = gml.E.reivscale[lev] + "tanh "
					gml.E.relabel[lev] 	= gml.E.relabel[lev] + "corr("+lats[r]+","+lats[ind]+") " 
					gml.E.Nreparams[lev] = gml.E.Nreparams[lev] + 1
				}
				ind++
			}			
		}		
	}
	//post to Stata
	st_local("vcveqns",invtokens(vcveqns))
	
	gml.Nb = eqn-1
	
}

`SS' merlin_get_vcv_vars(`gml' gml, `RS' lev)
{	
	vars = J(1,1,"")
	//check models for any re predictors
	lats = asarray(gml.latlevs,lev)		
	for (mod=1;mod<=gml.Nmodels;mod++) {
		if (gml.familys[mod]=="re") {
			re = st_local("re"+strofreal(mod))
			for (r=1;r<=cols(lats);r++) {
				if (re==lats[1,r]) {
					gml.hasRExb[lev] = 1
					vars = st_local("xb"+strofreal(mod))
				}
			}
		}
	} 
	return(vars)
}

void merlin_get_weights(`gml' gml, `SR' modeltouses, levelvars)
{
	
	nlevs = gml.Nlevels + 1
	gml.hasweights = J(nlevs,1,0)	//inc ob level
	
	if (st_local("weights")!="") {
	
		wts = tokens(st_local("weights"))
		if (cols(wts)>nlevs) {
			errprintf("Incorrect number of weight variables\n")
			exit(198)
		}
		gml.weights = asarray_create("real",2)
		
		//observation level done separately to allow separate model touses
		//-> note levels are indexed from highest 1, to second lowest gml.Nlevels in logL calcs,
		//-> but lowest to highest should be specified in weights() option
		gml.hasweights[nlevs] = 1
		for (i=1;i<=gml.Nmodels;i++) {
			asarray(gml.weights,(nlevs,i),st_data(.,wts[1,1],modeltouses[1,i]))
		}
		if (cols(wts)>1) {
			for (i=2;i<=cols(wts);i++) {
				nlevs--
				gml.hasweights[nlevs] = 1
				//make touse to get one row per cluster
				stata("tempvar wtsindex"+strofreal(i)+" wtstouse"+strofreal(i))
				stata("qui egen "+st_local("wtsindex"+strofreal(i))+"= group("+invtokens(levelvars[1,1..nlevs])+")")
				stata("qui bys "+st_local("wtsindex"+strofreal(i))+": gen byte "+st_local("wtstouse"+strofreal(i))+"=_n==1 if "+st_local("touse"))		
				asarray(gml.weights,(nlevs,1),st_data(.,wts[1,i],st_local("wtstouse"+strofreal(i))))				
			}
		}
	}
}

void merlin_get_n_dist_ancp(`gml' gml)
{
	familys = gml.familys
	gml.initdistapindex = J(2,0,.)
	gml.dapmodels = asarray_create("real",1)
	
	for (i=1;i<=gml.Nmodels;i++) {
		stri = strofreal(i)
		dapmodel = strtoreal(tokens(st_local("dapmodel"+stri)))
		if (dapmodel==J(1,0,.)) {
			dapmodel = 0
			
			if (familys[i]=="weibull" | familys[i]=="gaussian" | familys[i]=="beta" | familys[i]=="negbinomial" | familys[i]=="lquantile" | familys[i]=="gamma" | familys[i]=="gompertz") {
				gml.Ndistancp[i] = 1
				st_local("dap_eqns"+stri,"/dap"+stri+"_1")
			}
			else if (familys[i]=="rp" | familys[i]=="rcs" | familys[i]=="prp") {
				gml.Ndistancp[i] = cols(tokens(st_local("rcsvars"+stri)))
				st_local("dap_eqns"+stri,"(_rcs"+stri+":" +st_local("rcsvars"+stri)+", nocons)")
			}
			else if (familys[i]=="ordinal") {
				for (k=1; k<=gml.Ndistancp[i]; k++) {
					st_local("dap_eqns"+stri, st_local("dap_eqns"+stri)+" /dap"+stri+"_"+strofreal(k))
				}
			}
			else if (familys[i]=="gp") {
				gml.Ndistancp[i] = 2 + (st_local("noresidual"+strofreal(i))=="")
				for (k=1; k<=gml.Ndistancp[i]; k++) st_local("dap_eqns"+stri, st_local("dap_eqns"+stri)+" /dap"+stri+"_"+strofreal(k))
			}
			
		}
		asarray(gml.dapmodels,i,dapmodel)
	}
}

void merlin_get_xb_p(`gml' gml)
{
	`RS' nm
	`PC' pmat
	
	for (i=1;i<=gml.Nmodels;i++) {
		if (!(st_local("xb"+strofreal(i))=="" & st_local("constant"+strofreal(i))=="noconstant") & gml.familys[i]!="re") {
			gml.NotNull[i,2] = 1
		}
	}	
	
	familys = gml.familys
	nm = rows(familys)
	pmat = J(nm,1,&nm)
	
	for (i=1;i<=nm;i++) {
		if (familys[i]=="exponential" | familys[i]=="poisson" | familys[i]=="bernoulli" | familys[i]=="user") {
			pmat[i] = &merlin_xb_1()
		}
		else if (familys[i]=="weibull" | familys[i]=="gaussian" | familys[i]=="beta" | familys[i]=="negbinomial" | familys[i]=="lquantile" | familys[i]=="gamma") {
			pmat[i] = &merlin_xb_2()
		}
		else if (familys[i]=="gompertz") {
			pmat[i] = &merlin_xb_3()
		}
		else if (familys[i]=="rcs" | familys[i]=="rp" | familys[i]=="prp") {
			pmat[i] = &merlin_xb_rcs()
		}
		else if (familys[i]=="null" & gml.NotNull[i,2]) {
			pmat[i] = &merlin_xb_1()
		}
		else if (familys[i]=="ordinal") {
			pmat[i] = &merlin_xb_ord()
		}
		else if (familys[i]=="gp") {
			pmat[i] = &merlin_xb_gp()
		}
	}
	gml.Pxb = pmat
}

void merlin_get_logl_p(`gml' gml)
{
	`RS' nm
	`PC' pmat
	
	familys = gml.familys
	links 	= gml.links
	nm 		= rows(familys)
	pmat = pmatlinks = pmathaz = pmatloghaz = pmatchaz = J(nm,1,&nm)
	
	//get left truncation pointers to survival functions
	if (gml.hasanyltrunc) ltruncpmat = pmat

	gml.E.userfunctions = J(nm,4,"")
	
	for (i=1; i<=nm; i++) {
	
		pmatlinks[i] = &merlin_identity()			
		
		if (familys[i]=="exponential") {
			pmat[i] = &merlin_logl_exp()
			if (gml.hasltrunc[i]) 			ltruncpmat[i] = &merlin_p_exp_ch()
		}
		else if (familys[i]=="weibull") 	{
			if (st_local("mweight"+strofreal(i))=="") 	pmat[i] = &merlin_logl_weibull()
			else 										pmat[i] = &merlin_p_weibull_s()
			if (gml.hasltrunc[i]) 			ltruncpmat[i] = &merlin_p_weibull_ch()
		}
		else if (familys[i]=="gompertz") {
			pmat[i] = &merlin_logl_gompertz()
			if (gml.hasltrunc[i]) 			ltruncpmat[i] = &merlin_p_gompertz_ch()
		}
		else if (familys[i]=="gaussian") 	{
			pmat[i] = &merlin_logl_gaussian()
			if (links[i]!="identity") 		pmatlinks[i] = &merlin_exp()
		}
		else if (familys[i]=="poisson") {
			pmat[i] = &merlin_logl_poisson()
			pmatlinks[i] = &merlin_exp()
		}
		else if (familys[i]=="bernoulli") {
			pmat[i] = &merlin_logl_bernoulli()
			if (links[i]=="logit") 			pmatlinks[i] = &merlin_invlogit()
			else if (links[i]=="probit") 	pmatlinks[i] = &merlin_normal()
			else 							pmatlinks[i] = &merlin_invcll()
		}
		else if (familys[i]=="beta") {
			pmat[i] = &merlin_logl_beta()
		}
		else if (familys[i]=="negbinomial")	{
			pmat[i] = &merlin_logl_negbinomial()
		}
		else if (familys[i]=="ordinal")	{
			if (links[i]=="logit") 			pmat[i] = &merlin_logl_ologit()
			else 							pmat[i] = &merlin_logl_oprobit()
		}
		else if (familys[i]=="lquantile") 	{
			pmat[i] = &merlin_logl_qtile()
		}
		else if (familys[i]=="rp") {
			pmat[i] = &merlin_logl_rp()
			if (gml.hasltrunc[i]) ltruncpmat[i] = &merlin_p_rp_ch()
		}
		else if (familys[i]=="prp") {
			gml.Npenals = gml.Npenals :+ 1
			pmat[i] = &merlin_logl_prp()
			if (gml.hasltrunc[i]) ltruncpmat[i] = &merlin_p_rp_ch()
		}
		else if (familys[i]=="rcs")	{
			pmat[i] = &merlin_logl_rcs()
			if (gml.hasltrunc[i]) ltruncpmat[i] = &merlin_p_rcs_ch()
		}
		else if (familys[i]=="gamma")	{
			pmat[i] = &merlin_logl_gamma()
			pmatlinks[i] = &merlin_exp()
		}
		else if (familys[i]=="user")		{
			if (st_local("hazfunction"+strofreal(i))=="" & st_local("loghazfunction"+strofreal(i))=="") {
				stata("mata: loglf = &"+st_local("loglfunction"+strofreal(i))+"()")
				external loglf
				pmat[i] = loglf
				gml.E.userfunctions[i,1] = st_local("loglfunction"+strofreal(i))
			}
			else {
				if (st_local("hazfunction"+strofreal(i))!="") {
					stata("mata: hazf = &"+st_local("hazfunction"+strofreal(i))+"()")
					external hazf
					pmathaz[i] = hazf
					gml.E.userfunctions[i,2] = st_local("hazfunction"+strofreal(i))
					
					if (st_local("chazfunction"+strofreal(i))=="") {
						pmat[i] = &merlin_logl_userhaz()
						if (gml.hasltrunc[i]) ltruncpmat[i] = &merlin_p_userh_ch()
					}
					else {
						stata("mata: chazf = &"+st_local("chazfunction"+strofreal(i))+"()")
						external chazf
						pmatchaz[i] = chazf
						pmat[i] = &merlin_logl_userhazchaz()
						gml.E.userfunctions[i,3] = st_local("chazfunction"+strofreal(i))
						if (gml.hasltrunc[i]) ltruncpmat[i] = &merlin_p_userhch_ch()
					}
				}
				else {
					stata("mata: hazf = &"+st_local("loghazfunction"+strofreal(i))+"()")
					external hazf
					pmatloghaz[i] = hazf
					gml.E.userfunctions[i,4] = st_local("loghazfunction"+strofreal(i))
					if (st_local("chazfunction"+strofreal(i))=="") {
						pmat[i] = &merlin_logl_userloghaz()
						if (gml.hasltrunc[i]) ltruncpmat[i] = &merlin_p_userlogh_ch()
					}
					else {
						stata("mata: chazf = &"+st_local("chazfunction"+strofreal(i))+"()")
						external chazf
						pmatchaz[i] = chazf
						pmat[i] = &merlin_logl_userloghazchaz()
						gml.E.userfunctions[i,3] = st_local("chazfunction"+strofreal(i))
						if (gml.hasltrunc[i]) ltruncpmat[i] = &merlin_p_userloghch_ch()
					}
				}
			}
		}
		else if (familys[i]=="null") {
			if 		(links[i]=="logit") pmatlinks[i] = &merlin_invlogit()
			else if (links[i]=="atanh") pmatlinks[i] = &merlin_tanh()
		}
		else if (familys[i]=="gp") {
			if (gml.Ndistancp[i]==2) 	pmat[i] = &merlin_logl_gp_noresid()
			else 						pmat[i] = &merlin_logl_gp()
		}
		
	}
	gml.Plnl 		= pmat
	gml.userhaz 	= pmathaz
	gml.userloghaz 	= pmatloghaz
	gml.userchaz 	= pmatchaz
	gml.invlinks 	= pmatlinks
	gml.ltruncP		= ltruncpmat
}

void merlin_ordinal_setup(`gml' gml)
{
	if (sum(gml.familys:=="ordinal")) {
		gml.OrdIndexes = asarray_create("real",2)
		for (i=1;i<=gml.Nmodels;i++)  {
			if (gml.familys[i]=="ordinal") {
			
				//get Ndf from response var i.e. # of cut points
				y 		= asarray(gml.y,i)
				ks 		= uniqrows(y)
				nks 	= rows(ks)
				kindex 	= 1::nks
				oindex 	= J(gml.Nobs[gml.Nlevels,i],2,.)
				for (j=1;j<=gml.Nobs[gml.Nlevels,i];j++) {
					oindex[j,1] = selectindex(ks:==y[j])
				}
				oindex[,2] 			= oindex[,1] :- 1
				Ndf 				= rows(ks) - 1
				gml.Ndistancp[i] 	= Ndf
				
				asarray(gml.OrdIndexes,(i,4),oindex)
				asarray(gml.OrdIndexes,(i,1),selectindex(oindex[,1]:==1))
				asarray(gml.OrdIndexes,(i,2),selectindex((oindex[,1]:!=1) :* (oindex[,2]:!=Ndf)))
				asarray(gml.OrdIndexes,(i,3),selectindex(oindex[,2]:==Ndf))

			}
		}
	}
}

void merlin_ancp_init(`gml' gml)
{
	gml.Nap = J(gml.Nmodels,1,0)
	gml.initapindex = J(2,0,.)
	for (i=1;i<=gml.Nmodels;i++)  {
		gml.Nap[i] 	= strtoreal(st_local("nap"+strofreal(i)))
		apeqns		= J(1,0,"")
		for (a=1;a<=gml.Nap[i];a++) {
			apeqns = apeqns," " + "/ap"+strofreal(i)+"_"+strofreal(a)
		}
		st_local("ap_eqns"+strofreal(i),invtokens(apeqns))
	}
	
	if (sum(gml.Nap)) gml.apxb = asarray_create("real",2)
}

void merlin_hazNI_init(`gml' gml)
{
	if (sum(gml.NI)) {
		gml.hazNnodes = J(gml.Nmodels,1,0)
		gml.haznodes = gml.hazweights = asarray_create("real",1)
		for (i=1;i<=gml.Nmodels;i++) {
			if (gml.NI[i]) {
				gml.hazNnodes[i] = 30		//!!
				gq = merlin_gq(gml.hazNnodes[i],"legendre")
// 				gq	= (0.991455371120813,-0.991455371120813,0.949107912342759,-0.949107912342759,0.864864423359769,-0.864864423359769,0.741531185599394,-0.741531185599394,0.586087235467691,-0.586087235467691,0.405845151377397,-0.405845151377397,0.207784955007898,-0.207784955007898,0)'
// 				gq = gq,( 0.022935322010529,0.022935322010529,0.063092092629979,0.063092092629979,0.104790010322250,0.104790010322250,0.140653259715525,0.140653259715525,0.169004726639267,0.169004726639267,0.190350578064785,0.190350578064785,0.204432940075298,0.204432940075298,0.209482141084728)'
				
				//need response var stime
				ys = asarray(gml.y,i)
				asarray(gml.haznodes,i, ys[,1] :* J(gml.Nobs[gml.Nlevels,i],1,gq[,1]'):/2 :+ ys[,1]:/2)
				asarray(gml.hazweights,i,ys[,1] :* J(gml.Nobs[gml.Nlevels,i],1,gq[,2]'):/2)
			}
		}
	}
}

void merlin_init_numint(`gml' gml)
{

	familys = gml.familys
	gml.NI = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		if (familys[i]=="exponential" | familys[i]=="weibull" | familys[i]=="gompertz" | (familys[i]=="user" & st_local("hazfunction"+strofreal(i))!="") | (familys[i]=="user" & st_local("loghazfunction"+strofreal(i))!="") | (familys[i]=="user" & st_local("failure"+strofreal(i))!="")) {
			if (st_local("timevar"+strofreal(i))!="") gml.NI[i] = 1
		}
		
		if (familys[i]=="user" & (st_local("hazfunction"+strofreal(i))!="" | st_local("loghazfunction"+strofreal(i))!="") & st_local("chazfunction"+strofreal(i))=="") {
			gml.NI[i] = 1
		}
		if (familys[i]=="rcs") {
			gml.NI[i] = 1
		}
	}

}

void merlin_noconstants(`gml' gml) 
{
	gml.hascons = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.familys[i]=="ordinal") gml.hascons[i] = 0
		else gml.hascons[i] = st_local("constant"+strofreal(i))!="noconstant"
	}
}

void merlin_rcs_setup(`gml' gml)
{
	for (i=1;i<=gml.Nmodels;i++) {
		gml.model = i
		rcsopts = st_local("rcsopts"+strofreal(i))
		if (rcsopts!="") merlin_rcs_build(gml,rcsopts)
	}
}

void merlin_rcs_build(`gml' gml, `SS' rcsopts)
{	

	y = asarray(gml.y,gml.model)

	if (gml.predict) { 													//called from predict
		strk 		= strofreal(gml.model)
		knots 		= strtoreal(tokens(st_global("e(rcs"+strk+")")))
		hasorthog 	= st_global("e(orthog"+strk+")")!=""
		if (hasorthog) {
			rmat 	= st_matrix("e(rcsrmat_"+strk+")")
		}
	}
	else {

		//parse options
		stata("local 0 , "+rcsopts)
		stata("syntax , [DF(string) Knots(string) NOORTHog]")
		
		knots 		= st_local("knots")
		df 			= st_local("df")
		hasorthog 	= st_local("noorthog")==""
		
		if (knots!="" & df!="") merlin_error("Can't specify both df() and knots()")
		
		
		
		if (df!="") {
			df 		= strtoreal(df)
			tv 		= select(y[,1],y[,2])
			tv 		= log(tv)
			tv	 	= sort(tv,1)
			nrows 	= rows(tv)
			if (df==1) 			index = 1\nrows
			else {
				if (df==2) 		index = 50
				else if (df==3) index = 33.3333333333\66.66666666
				else if (df==4) index = 25\50\75
				else if (df==5) index = 20\40\60\80
				else if (df==6) index = 17\33\50\67\83
				else if (df==7) index = 14\29\43\57\71\86
				else if (df==8) index = 12.5\25\37.5\50\62.5\75\87.5
				else if (df==9) index = 11.1\22.2\33.3\44.4\55.6\66.7\77.8\88.9
				else if (df==10) index = 10\20\30\40\50\60\70\80\90 
				index = 1\round(index :/100 :* nrows)\nrows
			}
			knots = tv[index]'
		}
		else knots = strtoreal(tokens(knots))

		if (hasorthog)	rmat 	= merlin_orthog(merlin_rcs(log(y[,1]),knots)) 
	}

	if (hasorthog) 	rcsvars = merlin_rcs(log(y[,1]),knots,0,rmat)
	else 			rcsvars = merlin_rcs(log(y[,1]),knots)
	
	//store stuff
	asarray(gml.distancb,(gml.model,2),rcsvars)
	asarray(gml.distancb,(gml.model,3),knots)
	asarray(gml.distancb,(gml.model,4),hasorthog)
	if (hasorthog) asarray(gml.distancb,(gml.model,5),rmat)

	if (gml.familys[gml.model]=="rp" | gml.familys[gml.model]=="prp") {
		if (hasorthog) 	drcsvars = merlin_rcs(log(y[,1]),knots,1,rmat)
		else 			drcsvars = merlin_rcs(log(y[,1]),knots,1)
		asarray(gml.distancb,(gml.model,6),drcsvars) 
	}
	
	//post to Stata and ml equation local
	stub = "_rcs"+strofreal(gml.model)+"_"
	stata("cap drop "+stub+"*")
	Nvars = cols(rcsvars)
	names = J(1,0,"")
	eqnames = J(1,0,"")
	for (r=1;r<=Nvars;r++) {
		names = names,(stub+strofreal(r))
	}

	id = st_addvar("double",names)
	st_store(.,id,gml.modeltouses[gml.model],rcsvars)
	printf("variables created: "+stub+"1 to "+stub+strofreal(Nvars)+"\n")

	st_local("rcsvars"+strofreal(gml.model),invtokens(names))
	
}

void merlin_setup_EV_p(`gml' gml)
{
	pmat = J(gml.Nmodels,1,&gml.Nmodels)
	for (i=1; i<=gml.Nmodels; i++) {
		f = gml.familys[i]
		if 		(f=="gaussian")		pmat[i] = &merlin_gaussian_expval()
		else if (f=="bernoulli") 	pmat[i] = &merlin_bernoulli_expval()
		else if (f=="exp") 			pmat[i] = &merlin_exp_expval()
		else if (f=="gamma")		pmat[i] = &merlin_gamma_expval()
		else if (f=="weibull")		pmat[i] = &merlin_weibull_expval()
		else if (f=="poisson")		pmat[i] = &merlin_poisson_expval()
		else if (f=="beta")			pmat[i] = &merlin_beta_expval()
		else if (f=="null")			pmat[i] = &merlin_null_expval()
	}
	gml.expvalP = pmat
}

`RM' merlin_null_expval(`gml' gml, | `RC' t)
{
	if (args()==1) 	return((*gml.invlinks[gml.model])(merlin_util_xzb(gml)))
	else 			return((*gml.invlinks[gml.model])(merlin_util_xzb(gml,t)))	
}

void merlin_penalties(`gml' gml)
{
	gml.haspenalty = J(3,1,0)
	pen = st_local("penalty")
	gml.haspenalty[1] = pen!=""
	
	if (gml.haspenalty[1]) {
	
		if (gml.Npenals) merlin_error("penalty() can't be used with a model with ps() elements")
		
		if 		(pen=="lasso") gml.haspenalty[2] = 1
		else if (pen=="ridge") gml.haspenalty[2] = 2
		else merlin_error("Unknown penalty()")
		
		if (st_local("lambda")=="") {
			gml.haspenalty[3] = 0.1
		}
		else gml.haspenalty[3] = strtoreal(st_local("lambda"))
	}
	
}

void merlin_gp_setup(`gml' gml)
{
// 	gml.gpresid = J(gml.Nmodels,1,.)
// 	for (i=1; i<=gml.Nmodels; i++) gml.gpresid[i] = st_local("noresidual"+strofreal(i))==""
	
}

end


