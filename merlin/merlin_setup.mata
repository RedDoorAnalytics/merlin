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

version 15.1

mata:

struct merlin_ereturn_struct {
	`SR' levelvars		//level vars 
	`RC' Nreparams		//# of VCV parameters at each level
	`SC' reeqns		//eqn names for re parameters
	`SC' reivscale		//inverse scale of re parameters
	`SC' relabel		//labels for re parameters
	`SM' userfunctions			
	`SC' Nvars
	`RS' merlin
}

struct merlin_struct {

// [op] - denotes an option which is not always completely filled in

	`pgml' Pgml			//self-pointer
	
	//model definition
	`TR' y				//response variables
	`SC' familys			//family of each model
	`SC' responses			//response var names 
                                        // -> includes event and ltruncated if there
	`SC' links			//link function of each model
	`PC' invlinks			//inverse link function of each model
	
	`RC' issurv
	`SC' failures		//[op]	//failure variable name of each model 			-> for eret list
	`SC' ltruncated		//[op]	//ltruncated variable name of each model		-> for eret list
	`SC' linterval		//[op]	//linterval variable name of each model			-> for eret list
	`TR' surv_index
	`RM' Nsurv
	`RS' survind
	
	`RC' hasltrunc			//has left truncation/delayed entry -> handle with counting process notation
	`RS' hasanyltrunc		//flag for any counting process left truncation/delayed entry
        `RS' hasmargltrunc              //requires marginal ltruncation calc.
	`RS' mltmodel			//which model spec'd the marginal ltruncation
	`RC' mltllindex			//index for subbing in the logl contrib. at the cluster level
	
	`RC' haslint			//has interval censoring
	`RS' hasanylint			//flag for any interval censoring
	
	`SS' touse			//global touse variable name
	`SR' modeltouses		//model specific touse variable names
	`SR' levelvars			//level varables from highest to lowest
	`RM' NotNull			//flag whether the model contributes to the log likelihood
	`RM' hasbh			//presence of bhazard or bhfile for each model
	`TR' bhazards			//stored bhazard variable for each model
	`SC' bhvarnames			//bhazard() varnames
	`RM' istimedep			//flag - complex syntax is time-dependent (either on timevar(), or survival outcome)
	`TR' timevars			//stored timevar()
	`SC' tvarnames			//for eret list
	
	`RS' Nmodels			//# of models
	`RS' Nlevels			//# of levels (includes ob. level)
	`RS' Nrelevels			//# of levels with random effects
	`RS' N				//# of observations (ob. level, main touse)
	`RC' Nobs			//# of observations for each model = level by models
	`TR' xbindex			//ob. level index for each outcome, based on main touse core index
	
	`RC' Nap			//[op]	//# of anciliary parameters in a user model
	`TR' apxb
	`TR' distancb
	`TR' OrdIndexes		//[op]	//resp index for cuts for ordinal model
	`RC' qtiles			//[op]	//quantiles, if specified
	
	`PC' Plnl, Pxb			//pointers to logL and linear predictor functions for each model
	`RS' todo			//gf0,1,2
	`PC' Ph				//[op]	//pointer to hazard functions
	`PC' Pch			//[op]	//pointer to cumulative hazard functions
	`PC' Plogh			//[op]	//pointer to log hazard functions
	`PC' Pcdf
	`PM' Puserst
	`RS' nocox			//for ltruncation, as it's handled within the logl function for a Cox model, not externally
	`TR' cox_index
	
	//complex syntax
	`RS' nogen			//do not generate Stata variables for elements
	`RR' myb			//copy ml's b vector so I can update it with my numeric scores
	`RC' Ncmps			//number of components in each complex syntax
	`TR' elindex			//index of each element within each component, for each model
	`TR' elvarname			//varname for each element (used in matching when at() from predict mediation etc.)
	`TR' Nels			//# of elements per component per model
	`RM' eqnindex			//start and end index (of overall b vector) for each model's complex predictor
	`RC' hascons			//whether each complex predictor has a constant term
	`TR' elinfo
	`RC' Ndistancp			//# of distributional ancillary parameters for each model
	`PC' expvalP
	`RS' Nb				//total number of parameters in model
	`RR' initbindex			//starting value indexes to post fixed effects model results into b for full model fit
	`RM' initdistapindex
	`RM' initapindex		//starting value indexes to post ap into b
	`RM' initvarindex		//index for random effect variance parameters
	`RS' gridsearch
	`SC' cmplabels
	`TR' CmpXBIndex
	`TR' hasconstraint		//each model clp, whether each parameter is constraned 0/1
	`TR' eltvar			//RR for each mod,comp - flag for whether input var matches timevar(), ltruncated()
	
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
	`TR' c
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
	`RS' chip
	`TR' latlevs			//unique REs at each level for eret list
	`RS' adaptit			//number of iterations to allow adapting nodes
	
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
	`RS' fixedonly			//predictions based only on fixed effects - so skip res in utils
	`TR' blups
	
	//imputation models
	`RM' IsImputed			//flag colvector - model is for imputation
	`RS' hasImputed			//any imputation models
	`SC' ImputeVarname		//name of variable to be imputed in that model
	`TR' ImputeIndex		//stores index for imputation models
	`RS' imputing			//flag to extract correct index
	`RS' ImputedValue		//holds the value being integrated over
	`RM' NsurvImp			//obs numbers of survival outcomes of imputed rows
	`TR' surv_index_imp		//indexes for survival events etc for imputed rows
	`RS' ImputeIP
	
	//ereturn list
	`Egml' E
	
	`RC' simple			//flag for simple design matrix -> no [] elements
	`TR' X				//clp design matrix for simple models
	`TR' dX				//d/dt of clp design matrix for simple models
	`TR' X_bindex
	
	//score & hessian
	`RC' NHeqns
	`TR' NHbs
	`RC' skip			//number of parameters to skip (fixed effect models only)
	
	//====================================================================================================//
	//development
	
        `SS' indicatorvar               //indicator() varname from stexcess
        `RC' indicator                  //indicator view for stexcess
        
	//competing risks with interval censoring
	`RS' hastmat, iccrfudge
	`RM' tmat
	
	//sample weights, modelled weights
	`RS' hasweights
	`RS' Nmodels2
	`TR' weights
	
	//for predictms
	`SS' allvars
	`RS' offsetflag
	
	//overall penalty
	`RS' haspenalty			//flag
	`RS' islasso			//lasso or ridge
	`RS' lambda			//pen. parameter
	
	//morgana
	`RS' morgana	
}

void merlin_setup(`SS' GML,`SS' touse)
{
	//declarations
	`pgml' 	pGML
	`gml' 	gml
	
	//initialise
	gml.Pgml = pGML = crexternal(GML)
	
	//core setup
	merlin_setup_core(gml,touse)
	merlin_setup_family_general(gml)
	merlin_setup_check_clp(gml)
	merlin_setup_levels(gml)
	merlin_setup_error_checks(gml)
	merlin_setup_simple_flag(gml)
	merlin_setup_devstuff(gml)
	merlin_setup_devcodes(gml)
	merlin_build_touses(gml)		//can now read stuff in
	merlin_get_levelvars(gml)
	merlin_setup_evaltype(gml)
	merlin_get_ys(gml)
	merlin_setup_family_specific(gml)
        
        gml.indicatorvar = st_local("indicator")
        if (gml.indicatorvar!="") {
                st_view(gml.indicator=.,.,gml.indicatorvar,gml.touse)
        }
	
	//clp
        merlin_setup_EV_p(gml)                  //needed in case called from a user element
	merlin_get_noconstants(gml)
        merlin_get_timevars(gml)
        merlin_build_clp(gml)			//updates istimedep
        merlin_get_cmp_labels(gml)

	//survival extras
	merlin_setup_survival(gml)
	merlin_setup_cox(gml)			//keep here, needs istimedep from build_clp()
	merlin_setup_rp(gml)
        merlin_setup_marg_ltruncated(gml)
	
	//latents
	merlin_setup_latents(gml)

	//pointers
	merlin_get_xb_p(gml)
	merlin_get_logl_p(gml)

	//ml stuff
	merlin_get_weights(gml)
	merlin_setup_gf12(gml)
	merlin_get_dap_mleqns(gml)
	merlin_setup_mleqns(gml)
	merlin_setup_wrappers(gml)
	merlin_starting_values(gml)

	//Done
	swap(*pGML,gml)
}

void merlin_setup_core(`gml' gml, `SS' touse)
{
	gml.Nmodels	= strtoreal(st_local("neq"))
	gml.touse 	= touse
	gml.model 	= 1	//gets updated
	gml.modtoind 	= 1	//gets updated
	gml.fixedonly	= 0
	gml.iter 	= 0
	gml.offsetflag 	= 0
	gml.survind	= 0
	gml.imputing	= 0
	gml.nogen	= st_local("nogen")!=""
	gml.predict 	= st_local("predict")!=""
	if (gml.predict & st_local("npredict")!="") gml.N = strtoreal(st_local("npredict"))	
}	

void merlin_setup_wrappers(`gml' gml)
{	
	if 	(st_local("bors")!="")	        gml.E.merlin	= 4  //stmerlin
	else if (st_local("arthur")!="") 	gml.E.merlin	= 3  //neuralnet
	else if (st_local("excalibur")!="")     gml.E.merlin	= 2  //stmixed
	else if (st_local("sagramore")!="")     gml.E.merlin	= 6  //jm
        else if (st_local("mordred")!="")       gml.E.merlin    = 7  //stexcess
	else {
		if (gml.Nrelevels)		gml.E.merlin	= 1  //merlin
		else 				gml.E.merlin	= 5
	}
	gml.morgana = st_global("c(prefix)")=="morgana"
}

void merlin_setup_family_general(`gml' gml)
{
	gml.familys 	= tokens(st_local("familylist"))'
	gml.links 	= tokens(st_local("linklist"))'		//defunct (mainly) -> to bring back
	gml.NotNull	= J(1,2,((gml.familys:!="null") 
                                :* (gml.familys:!="re")))	//contributes to logl, has main xb
	gml.issurv	= J(gml.Nmodels,1,0)
	gml.failures	= J(1,gml.Nmodels,"")
	gml.ltruncated	= J(1,gml.Nmodels,"")
	gml.hasltrunc	= J(gml.Nmodels,1,0)
	gml.linterval	= J(1,gml.Nmodels,"")
	gml.haslint	= J(gml.Nmodels,1,0)
	gml.hasbh 	= J(gml.Nmodels,2,0)
	gml.bhvarnames	= J(gml.Nmodels,1,"")
	for (i=1;i<=gml.Nmodels;i++) {
		gml.failures[1,i] 	= st_local("failure"+strofreal(i))
		gml.issurv[i]		= gml.failures[1,i]!=""
		gml.ltruncated[1,i]     = st_local("ltruncated"+strofreal(i))
		gml.hasltrunc[i] 	= gml.ltruncated[1,i]!=""
		gml.linterval[1,i]      = st_local("linterval"+strofreal(i))
		gml.haslint[i] 		= gml.linterval[1,i]!=""
		gml.bhvarnames[i]	= st_local("bhaz"+strofreal(i))
		gml.hasbh[i,1] 		= gml.bhvarnames[i]!=""	
		gml.hasbh[i,2] 		= st_local("bhfile"+strofreal(i))!=""
	}
	gml.hasanyltrunc = sum(gml.hasltrunc)
	gml.ltflag	 = 0			//gets updated when needed
	gml.hasanylint 	 = sum(gml.haslint)
	gml.hastmat 	 = st_local("transmatrix")!=""
	if (gml.hasanylint>1 & !gml.hastmat) {
                merlin_error("You can only have one model with interval-censoring")
        }
	gml.distancb 	= asarray_create("real",2)
	gml.Ndistancp 	= J(gml.Nmodels,1,0)
	merlin_setup_ap(gml)
	merlin_penalties(gml)
}

void merlin_setup_ap(`gml' gml)
{
	gml.Nap = J(gml.Nmodels,1,0)
	gml.initapindex = J(2,0,.)
	for (i=1;i<=gml.Nmodels;i++)  {
		gml.Nap[i] 	= strtoreal(st_local("nap"+strofreal(i)))
		apeqns		= J(1,0,"")
		for (a=1;a<=gml.Nap[i];a++) {
			apeqns = apeqns," " + "(ap"+strofreal(i)+"_"+strofreal(a)+":)"
		}
		st_local("ap_eqns"+strofreal(i),invtokens(apeqns))
	}
	if (sum(gml.Nap)) gml.apxb = asarray_create("real",2)
}

//development stuff
void merlin_setup_devstuff(`gml' gml)
{
	//transmatrix,  intcens
	if (gml.hastmat) merlin_setup_tmat(gml)
	gml.iccrfudge = st_local("iccrfudge")!=""
	
	merlin_setup_any_impute_models(gml)
}

//family parameters setup
void merlin_setup_family_specific(`gml' gml)
{	
	merlin_get_n_dist_ancp(gml)					//declares dist anc p 
	merlin_setup_pwexp(gml)
	merlin_setup_ordinal(gml)
	merlin_get_quantiles(gml)
	merlin_gp_setup(gml)
	merlin_setup_impute(gml)					//imodels
}

`SR' merlin_get_indepvars(`RS' i)
{
	indepvars = J(1,0,"")
	st_local("rest",st_local("indepvars"+strofreal(i)))
	stata("gettoken lhs rest : rest, bind")						
	indepvars = indepvars,st_local("lhs")
	while (st_local("rest")!="") {
		stata("gettoken lhs rest : rest, bind")
		indepvars = indepvars,st_local("lhs")
	}
	return(indepvars)
}

void merlin_get_ys(`gml' gml)
{
	gml.y = asarray_create("real",1)
	gml.responses = J(1,0,"")
	for (i=1;i<=gml.Nmodels;i++) {
		resp = st_local("response"+strofreal(i))
		if (resp!="") {
			gml.responses = gml.responses,resp
			asarray(gml.y,i,st_data(.,tokens(resp),gml.touse))
		}
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
	gml.tvarnames 	= J(gml.Nmodels,1,"")
	gml.timevars 	= asarray_create("real",1)		//must be declared - used in predict
	gml.istimedep	= J(gml.Nmodels,2,0)
	
	for (mod=1;mod<=gml.Nmodels;mod++) {
		timevar = st_local("timevar"+strofreal(mod))
		if (timevar!="") {
			gml.tvarnames[mod] 		= timevar
			gml.istimedep[mod,1] 	= 1
		}
		else if (gml.issurv[mod]) 	gml.tvarnames[mod] = tokens(gml.responses[mod])[1]
		
		if (gml.tvarnames[mod]!="") {
			asarray(gml.timevars,mod,st_data(.,gml.tvarnames[mod],gml.touse))
		}
	}
}

void merlin_init_vcvs(`gml' gml)
{
	gml.vcvs = asarray_create("real",1)
	for (i=1;i<=gml.Nlevels;i++) asarray(gml.vcvs,i,I(gml.Nres[i]))
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
		else 		gml.covariances = covs:=="diagonal"\covs:=="exchangeable"\covs:=="unstructured"
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
				vars = "" //merlin_get_vcv_vars(gml, lev)
				vcveqns = vcveqns,"(lns"+strlev+"_"+strofreal(r)+": "+vars+")"
				gml.E.Nreparams[lev] 	= gml.E.Nreparams[lev] + 1
				gml.E.reeqns[lev] 	= gml.E.reeqns[lev] + "lns"+strlev+"_"+strofreal(r)+" "
				gml.E.reivscale[lev]    = gml.E.reivscale[lev] + "exp "
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
			if (gml.Nres[lev]) {
				gml.initvarindex = gml.initvarindex,(eqn++\0)
				vars = "" //merlin_get_vcv_vars(gml, lev)
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
		}
		
		//exch and nres>1
		if (gml.covariances[2,lev] & gml.Nres[lev]>1) {
			gml.initvarindex = gml.initvarindex,(eqn++\0)
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
					gml.initvarindex = gml.initvarindex,(eqn++\0)
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

// `SS' merlin_get_vcv_vars(`gml' gml, `RS' lev)
// {	
// 	vars = J(1,1,"")
// 	//check models for any re predictors
// 	lats = asarray(gml.latlevs,lev)		
// 	for (mod=1;mod<=gml.Nmodels;mod++) {
// 		if (gml.familys[mod]=="re") {
// 			re = st_local("re"+strofreal(mod))
// 			for (r=1;r<=cols(lats);r++) {
// 				if (re==lats[1,r]) {
// 					gml.hasRExb[lev] = 1
// 					vars = st_local("xb"+strofreal(mod))
// 				}
// 			}
// 		}
// 	} 
// 	return(vars)
// }

void merlin_get_weights(`gml' gml)
{
	Nlevs 			= gml.Nlevels 
	gml.hasweights 	= J(Nlevs,1,0)	//inc ob level
	
	if (st_local("weights")!="") {
	
		wts = tokens(st_local("weights"))
		if (cols(wts)>Nlevs) {
			errprintf("Incorrect number of weight variables\n")
			exit(198)
		}
		gml.weights = asarray_create("real",2)
		
		//observation level done separately to allow separate model touses
		//-> note levels are indexed from highest 1, to second lowest gml.Nlevels in logL calcs,
		//-> but lowest to highest should be specified in weights() option
		gml.hasweights[Nlevs] = 1
		for (i=1;i<=gml.Nmodels;i++) {
			asarray(gml.weights,(Nlevs,i),st_data(.,wts[1,1],gml.modeltouses[1,i]))
		}
		if (cols(wts)>1) {
			for (i=2;i<=cols(wts);i++) {
				Nlevs--
				gml.hasweights[Nlevs] = 1
				//make touse to get one row per cluster
				stata("tempvar wtsindex"+strofreal(i)+" wtstouse"+strofreal(i))
				stata("qui egen "+st_local("wtsindex"+strofreal(i))+"= group("+invtokens(gml.levelvars[1,1..Nlevs])+")")
				stata("qui bys "+st_local("wtsindex"+strofreal(i))+": gen byte "+st_local("wtstouse"+strofreal(i))+"=_n==1 if "+st_local("touse"))		
				asarray(gml.weights,(Nlevs,1),st_data(.,wts[1,i],st_local("wtstouse"+strofreal(i))))				
			}
		}
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
		if (familys[i]=="exponential" | familys[i]=="poisson" | familys[i]=="bernoulli" | familys[i]=="ibernoulli" | familys[i]=="user") {
			pmat[i] = &merlin_xb_1()
		}
		else if (familys[i]=="weibull" | familys[i]=="gaussian" | familys[i]=="igaussian" | familys[i]=="beta" | familys[i]=="negbinomial" | familys[i]=="lquantile" | familys[i]=="gamma" | familys[i]=="lognormal" | familys[i]=="loglogistic") {
			pmat[i] = &merlin_xb_2()
		}
		else if (familys[i]=="gompertz") {
			pmat[i] = &merlin_xb_3()
		}
		else if (familys[i]=="ggamma") {
			pmat[i] = &merlin_xb_ggamma()
		}
		else if (familys[i]=="rp") {
			pmat[i] = &merlin_xb_rcs()
		}
		else if (familys[i]=="pwexponential") {
			pmat[i] = &merlin_xb_pwexp()
		}
		else if (familys[i]=="aft") {
			pmat[i] = &merlin_xb_aft()
		}
		else if (familys[i]=="cox") {
			pmat[i] = &merlin_xb_1()
		}
		else if (familys[i]=="loghazard" | familys[i]=="logchazard" | familys[i]=="plogchazard" | familys[i]=="addhazard") {
			pmat[i] = &merlin_xb_1()
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
	`RS' nm, predictms
	`PC' pmat
	
	familys 	= gml.familys
	links 		= gml.links
	nm 			= rows(familys)
	pmat 		= pmatlinks = J(nm,1,&nm)
	pmathaz 	= pmat
	pmatloghaz 	= pmat
	pmatchaz 	= pmat
	pmatcdf 	= pmat
	pmatuserst	= J(nm,3,&nm)					//h,logh,ch
	
	predictms 	= st_local("galahad")!=""
	
	gml.E.userfunctions = J(nm,4,"")
	
	for (i=1; i<=nm; i++) {
	
		pmatlinks[i] = &merlin_identity()			
		
		if (familys[i]=="exponential") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
			else			pmat[i] = &merlin_logl_exp()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_exp_logh()
			else		pmatloghaz[i] 	= &merlin_exp_logh()
			pmatchaz[i] 	= &merlin_exp_ch()
			pmatcdf[i] 	= &merlin_exp_cdf()
		}
		else if (familys[i]=="pwexponential") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
			else 			pmat[i] = &merlin_logl_pwexp()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_pwexp_logh()
			else		pmatloghaz[i] 	= &merlin_pwexp_logh()
			pmatchaz[i] 	= &merlin_pwexp_ch()
			pmatcdf[i] 	= &merlin_pwexp_cdf()
		} 
		else if (familys[i]=="weibull") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_weibull_ml()
			else			pmat[i] = &merlin_logl_weibull()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_weibull_logh()
			else		pmatloghaz[i] 	= &merlin_weibull_logh()
			pmatchaz[i] 	= &merlin_weibull_ch()
			pmatcdf[i] 	= &merlin_weibull_cdf()
		}
		else if (familys[i]=="gompertz") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
			else 			pmat[i] = &merlin_logl_gompertz()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_gompertz_logh()
			else  		pmatloghaz[i] 	= &merlin_gompertz_logh()
			pmatchaz[i] 	= &merlin_gompertz_ch()
			pmatcdf[i] 	= &merlin_gompertz_cdf()
		}
		else if (familys[i]=="lognormal") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_lognormal_ml()
			else 			pmat[i] = &merlin_logl_survival_ob()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_lognormal_logh()
			else		pmatloghaz[i] 	= &merlin_lognormal_logh()
			pmatchaz[i] 	= &merlin_lognormal_ch()
		}
		else if (familys[i]=="loglogistic") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_loglogistic_ml()
			else 			pmat[i] = &merlin_logl_survival_ob()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_loglogistic_logh()
			else		pmatloghaz[i] 	= &merlin_loglogistic_logh()
			pmatchaz[i] 	= &merlin_loglogistic_ch()
		}
		else if (familys[i]=="gaussian" | familys[i]=="igaussian") 	{
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_gaussian_ml()
			else			pmat[i] = &merlin_logl_gaussian()
			if (links[i]!="identity") 	pmatlinks[i] = &merlin_exp()
		}
		else if (familys[i]=="poisson") {
			if (gml.Nrelevels)      pmat[i] = &merlin_logl_poisson_ml()
                        else                    pmat[i] = &merlin_logl_poisson()
			pmatlinks[i] = &merlin_exp()
		}
		else if (familys[i]=="bernoulli" | familys[i]=="ibernoulli") {
			if (gml.Nrelevels)      pmat[i] = &merlin_logl_bernoulli_ml()
                        else                    pmat[i] = &merlin_logl_bernoulli()
			if (links[i]=="logit") 		pmatlinks[i] = &merlin_invlogit()
			else if (links[i]=="probit") 	pmatlinks[i] = &merlin_normal()
			else 				pmatlinks[i] = &merlin_invcll()
		}
		else if (familys[i]=="beta") {
			pmat[i] = &merlin_logl_beta()
		}
		else if (familys[i]=="negbinomial")	{
			pmat[i] = &merlin_logl_negbinomial()
		}
		else if (familys[i]=="ordinal")	{
			if (links[i]=="logit") 	pmat[i] = &merlin_logl_ologit()
			else 			pmat[i] = &merlin_logl_oprobit()
		}
		else if (familys[i]=="lquantile") 	{
			pmat[i] = &merlin_logl_qtile()
		}
		else if (familys[i]=="logchazard") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
			else 			pmat[i] = &merlin_logl_logchazard()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_logchazard_logh()
			else		pmatloghaz[i] 	= &merlin_logchazard_logh()
			pmatchaz[i] 	= &merlin_logchazard_ch()
			pmatcdf[i] 	= &merlin_logchazard_cdf()
		}
		else if (familys[i]=="rp") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
			else 			pmat[i] = &merlin_logl_rp()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_rp_logh()
			else		pmatloghaz[i] 	= &merlin_rp_logh()
			pmatchaz[i] 	= &merlin_rp_ch()
			pmatcdf[i] 	= &merlin_rp_cdf()
		}
		else if (familys[i]=="aft") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
                        else                    pmat[i] = &merlin_logl_survival_ob()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_aft_logh()
			else		pmatloghaz[i] 	= &merlin_aft_logh()
			pmatchaz[i] 	= &merlin_aft_ch()
			pmatcdf[i] 	= &merlin_aft_cdf()
		}
		else if (familys[i]=="cox") {
			pmat[i] 		= &merlin_logl_cox()
		}
		else if (familys[i]=="loghazard") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
			else			pmat[i] = &merlin_logl_loghazard()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_loghazard_logh()
			else		pmatloghaz[i] 	= &merlin_loghazard_logh()
			pmatchaz[i] 	= &merlin_loghazard_ch()
			pmatcdf[i] 	= &merlin_loghazard_cdf()
		}
		else if (familys[i]=="addhazard") {
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_survival()
			else			pmat[i] = &merlin_logl_addhazard()
			if (predictms) 	pmatloghaz[i] 	= &merlin_p_addhazard_logh()
			else		pmatloghaz[i] 	= &merlin_addhazard_logh()
			pmatchaz[i] 	= &merlin_addhazard_ch()
			pmatcdf[i] 	= &merlin_addhazard_cdf()
		}
		else if (familys[i]=="ggamma")	{
			if (gml.Nrelevels) 	pmat[i] = &merlin_logl_ggamma_ml()
			else			pmat[i] = &merlin_logl_survival_ob()
			if (predictms) 	        pmatloghaz[i] 	= &merlin_p_ggamma_logh()
			else			pmatloghaz[i] 	= &merlin_ggamma_logh()
			pmatchaz[i] 	= &merlin_ggamma_ch()
		}
		else if (familys[i]=="gamma")	{
			pmat[i] = &merlin_logl_gamma()
			pmatlinks[i] = &merlin_exp()
		}
		else if (familys[i]=="user")		{
			if (st_local("hazfunction"+strofreal(i))=="" & 
                                st_local("loghazfunction"+strofreal(i))=="") {
				stata("mata: loglf = &"+st_local("loglfunction"+strofreal(i))+"()")
				external loglf
				pmat[i] = loglf
				gml.E.userfunctions[i,1] = st_local("loglfunction"+strofreal(i))
			}
			else {
				if (gml.Nrelevels)      pmat[i] = &merlin_logl_survival()
                                else                    pmat[i] = &merlin_logl_survival_ob()
				if (st_local("hazfunction"+strofreal(i))!="") {
					stata("mata: hazf = &"+st_local("hazfunction"+strofreal(i))+"()")
					external hazf
					pmatuserst[i,1] = hazf
					if (predictms) 	pmatloghaz[i] 	= &merlin_p_userh_logh()
					else		pmatloghaz[i] 	= &merlin_userhaz_logh()
					gml.E.userfunctions[i,2] = st_local("hazfunction"+strofreal(i))
					
					if (st_local("chazfunction"+strofreal(i))!="") {
						stata("mata: chazf = &"+st_local("chazfunction"+strofreal(i))+"()")
						external chazf
						pmatuserst[i,3] = chazf
						pmatchaz[i]    	= &merlin_userhazchaz_ch()
						pmatcdf[i]  	= &merlin_userhazchaz_cdf()
						gml.E.userfunctions[i,3] = st_local("chazfunction"+strofreal(i))
					}
					else {
						pmatchaz[i] = &merlin_userhaz_ch()
						pmatcdf[i]  = &merlin_userhaz_cdf()
					}
				}
				else {
					stata("mata: hazf = &"+st_local("loghazfunction"+strofreal(i))+"()")
					external hazf
					pmatuserst[i,2] = hazf
					if (predictms) 	pmatloghaz[i] 	= &merlin_p_userlogh_logh()
					else		pmatloghaz[i]	= &merlin_userloghaz_logh()
					gml.E.userfunctions[i,4] = st_local("loghazfunction"+strofreal(i))
					if (st_local("chazfunction"+strofreal(i))!="") {
						stata("mata: chazf = &"+st_local("chazfunction"+strofreal(i))+"()")
						external chazf
						pmatuserst[i,3] = chazf
						pmatchaz[i]    	= &merlin_userloghazchaz_ch()
						pmatcdf[i]  	= &merlin_userloghazchaz_cdf()
						gml.E.userfunctions[i,3] = st_local("chazfunction"+strofreal(i))
					}
					else {
						pmatchaz[i] = &merlin_userloghaz_ch()
						pmatcdf[i] 	= &merlin_userloghaz_cdf()
					}
				}
			}
		}
		else if (familys[i]=="null") {
			if 	(links[i]=="logit") pmatlinks[i] = &merlin_invlogit()
			else if (links[i]=="atanh") pmatlinks[i] = &merlin_tanh()
		}
		else if (familys[i]=="gp") {
			if (gml.Ndistancp[i]==2)        pmat[i] = &merlin_logl_gp_noresid()
			else 				pmat[i] = &merlin_logl_gp()
		}
		
	}
	
	gml.Plnl 	 = pmat
	gml.invlinks    = pmatlinks
	gml.Ph 		= pmathaz
	gml.Plogh 	= pmatloghaz
	gml.Pch 	= pmatchaz
	gml.Pcdf	= pmatcdf
	gml.Puserst     = pmatuserst
}

void merlin_get_n_dist_ancp(`gml' gml)
{
	familys = gml.familys
	gml.initdistapindex = J(2,0,.)
	
	for (i=1;i<=gml.Nmodels;i++) {
		stri = strofreal(i)
		if (familys[i]=="weibull" | familys[i]=="gaussian" | 
                        familys[i]=="igaussian" | familys[i]=="beta" | 
                        familys[i]=="negbinomial" | familys[i]=="lquantile" | 
                        familys[i]=="gamma" | familys[i]=="gompertz" | 
                        familys[i]=="lognormal" | familys[i]=="loglogistic") {
			gml.Ndistancp[i] = 1
		}
		else if (familys[i]=="ggamma") {
			gml.Ndistancp[i] = 2
		}
		else if (familys[i]=="rp") {
			rcsopts = st_local("rcsopts"+strofreal(i))
			stata("local 0 , "+rcsopts)
			stata("syntax , [DF(string) Knots(string) NOORTHog]")
			if (st_local("df")!="") gml.Ndistancp[i] = strtoreal(st_local("df"))
			else 			gml.Ndistancp[i] = cols(tokens(st_local("knots")))-1
		}
		else if (familys[i]=="aft") {
			rcsopts = st_local("rcsopts"+strofreal(i))
			stata("local 0 , "+rcsopts)
			stata("syntax , [DF(string) Knots(string) NOORTHog]")
			if (st_local("df")!="") {
                                gml.Ndistancp[i] = strtoreal(st_local("df"))+1
                        }
			else {
                                gml.Ndistancp[i] = cols(tokens(st_local("knots")))+2
                        }
		}
		else if (familys[i]=="gp") {
			gml.Ndistancp[i] = 2 + (st_local("noresidual"+strofreal(i))=="")
		}
	}
}

void merlin_get_dap_mleqns(`gml' gml)
{
	familys = gml.familys
	
	for (i=1;i<=gml.Nmodels;i++) {
		stri = strofreal(i)
		if (familys[i]=="weibull" | familys[i]=="gaussian" | 
                        familys[i]=="igaussian" | familys[i]=="beta" | 
                        familys[i]=="negbinomial" | familys[i]=="lquantile" | 
                        familys[i]=="gamma" | familys[i]=="gompertz" | 
                        familys[i]=="lognormal" | familys[i]=="loglogistic") {
                        gml.Ndistancp[i] = 1
			st_local("dap_eqns"+stri,"(dap"+stri+"_1:)")
		}
		else if (familys[i]=="ggamma") {
			st_local("dap_eqns"+stri,
                                "(dap"+stri+"_1:) (dap"+stri+"_2:)")
		}
		else if (familys[i]=="rp") {
			st_local("dap_eqns"+stri,"(_rcs"+stri+":" + 
                                st_local("rcsvars"+stri)+", nocons)")
		}
		else if (familys[i]=="aft") {
			st_local("dap_eqns"+stri,"(_rcs"+stri+":" +
                                st_local("rcsvars"+stri)+")")
		}
		else if (familys[i]=="ordinal" | familys[i]=="pwexponential") {
			for (k=1; k<=gml.Ndistancp[i]; k++) {
				st_local("dap_eqns"+stri, 
                                        st_local("dap_eqns"+stri)+
                                        " (dap"+stri+"_"+strofreal(k)+":)")
			}
		}
		else if (familys[i]=="gp") {
			for (k=1; k<=gml.Ndistancp[i]; k++) {
                                st_local("dap_eqns"+stri, 
                                        st_local("dap_eqns"+stri)+
                                        " (dap"+stri+"_"+strofreal(k)+":)")
                        }
		}
	}
}


void merlin_get_noconstants(`gml' gml) 
{
	gml.hascons = J(gml.Nmodels,1,0)
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.familys[i]=="ordinal" | gml.familys[i]=="cox") gml.hascons[i] = 0
		else gml.hascons[i] = st_local("constant"+strofreal(i))!="noconstant"
	}
}

void merlin_setup_EV_p(`gml' gml)
{
	pmat = J(gml.Nmodels,1,&gml.Nmodels)
	for (i=1; i<=gml.Nmodels; i++) {
		f = gml.familys[i]
		if 	(f=="gaussian")		pmat[i] = &merlin_gaussian_expval()
		else if (f=="bernoulli") 	pmat[i] = &merlin_bernoulli_expval()
		else if (f=="gamma")		pmat[i] = &merlin_gamma_expval()
		else if (f=="poisson")		pmat[i] = &merlin_poisson_expval()
		else if (f=="beta")		pmat[i] = &merlin_beta_expval()
		else if (f=="null")		pmat[i] = &merlin_null_expval()
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
	gml.haspenalty = st_local("penalty")!=""
	gml.islasso = st_local("penalty")=="lasso"
	if (gml.haspenalty) {
		if (st_local("lambda")=="") gml.lambda = 0.1
		else 						gml.lambda = strtoreal(st_local("lambda"))
	}
}

void merlin_gp_setup(`gml' gml)
{
// 	gml.gpresid = J(gml.Nmodels,1,.)
// 	for (i=1; i<=gml.Nmodels; i++) gml.gpresid[i] = st_local("noresidual"+strofreal(i))==""
}

void merlin_setup_tmat(`gml' gml) 
{
	gml.tmat     = st_matrix(st_local("transmatrix"))	
}

void merlin_setup_simple_flag(`gml' gml)
{
	gml.simple = J(gml.Nmodels,1,.)
	for (i=1;i<=gml.Nmodels;i++) {
		cmps = merlin_get_indepvars(i)
		gml.simple[i] = !any(strpos(cmps,"["))
	}
}

void merlin_setup_gf12(`gml' gml)
{
	
	if (gml.todo & !gml.predict) {
	
		gml.NHeqns 	= J(gml.Nmodels,1,0)
		NHbs 		= asarray_create("real",1)
		gml.skip	= J(gml.Nmodels,1,0)
		
		for (j=1;j<=gml.Nmodels;j++) {
			
			if 		(gml.familys[j]=="gaussian") {
				gml.NHeqns[j] = 2
				asarray(NHbs,j,((gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1)\1))
			}
			else if (gml.familys[j]=="exponential") {
				gml.NHeqns[j] = 1
				asarray(NHbs,j,(gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1))
			}
			else if (gml.familys[j]=="bernoulli") {
				gml.NHeqns[j] = 1
				asarray(NHbs,j,(gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1))
			}
			else if (gml.familys[j]=="weibull" | gml.familys[j]=="gompertz") {
				gml.NHeqns[j] = 2
				asarray(NHbs,j,((gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1)\1))
			}
			else if (gml.familys[j]=="loghazard" | gml.familys[j]=="addhazard") {
				gml.NHeqns[j] = 1
				asarray(NHbs,j,(gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1))
			}
			else if (gml.familys[j]=="logchazard") {
				gml.NHeqns[j] = 1
				asarray(NHbs,j,(gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1))
			}
			else if (gml.familys[j]=="rp" | gml.familys[j]=="pwexponential") {
				gml.NHeqns[j] = 2
				asarray(NHbs,j,((gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1)\gml.Ndistancp[j]))
			}
			else if (gml.familys[j]=="poisson") {
				gml.NHeqns[j] = 1
				asarray(NHbs,j,(gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1))
			}
			else if (gml.familys[j]=="cox") {
				gml.NHeqns[j] = 1
				asarray(NHbs,j,(gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1))
			}
			else if (gml.familys[j]=="user") {
				gml.NHeqns[j] = 1 + gml.Nap[j]
				nbetas = (gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1)
				if (gml.Nap[j]) nbetas = nbetas\J(gml.Nap[j],1,1)
				asarray(NHbs,j,nbetas)
			}
			else if (gml.familys[j]=="null") {
				gml.NHeqns[j] = 1
				asarray(NHbs,j,(gml.eqnindex[j,2] - gml.eqnindex[j,1] + 1))
			}
			
			if (j>1) {
				for (pm=1;pm<j;pm++) gml.skip[j] = gml.skip[j-1] + sum(asarray(NHbs,pm))
			}	
		}
		
		gml.NHbs = NHbs
	}
}

void merlin_setup_evaltype(`gml' gml)
{
	if (st_local("evaltype")=="") {
		gf = "gf0"
		if (sum(gml.simple)==gml.Nmodels) {
			check = J(gml.Nmodels,1,0)
			for (i=1;i<=gml.Nmodels;i++) {
				f = gml.familys[i]
				if (f=="addhazard" | f=="cox" | f=="exponential" | f=="weibull" | f=="rp" | f=="loghazard" | f=="logchazard" | f=="gaussian" | f=="poisson" | f=="pwexponential") {
					check[i] = 2
				}
				if ((f=="pwexponential" | f=="exponential" | f=="gompertz" | f=="weibull" | f=="loghazard") & gml.hasanylint) {
					check[i] = 1
				}
				if ((f=="exponential" | f=="weibull" | f=="logchazard" | f=="pwexponential") & gml.hasbh[i,1]) {
					check[i] = 0
				}
				if (gml.hasbh[i,2]) {
					check[i] = 0
				}
				if (f=="rp" & gml.hasbh[i,1]) {
					check[i] = 0
				}
			}
			gf = "gf"+strofreal(min(check))
		}
		if (gml.hastmat) gf = "gf0"
		st_local("evaltype",gf)
	}
	gml.todo = strtoreal(substr(st_local("evaltype"),3,1))
}

void merlin_setup_mleqns(`gml' gml)
{
	space = " "
	for (i=1; i<=gml.Nmodels; i++) {
		mlspeci	= st_local("xb"+strofreal(i))
		mlspeci	= mlspeci+space+st_local("dap_eqns"+strofreal(i))
		mlspeci	= mlspeci+space+st_local("ap_eqns"+strofreal(i))
		st_local("mlspec",st_local("mlspec")+space+mlspeci)
	}
	st_local("mlspec", st_local("mlspec")+space+st_local("vcveqns"))
}

void merlin_cleanup(`SS' GML)
{
	rmexternal(GML)
// 	`pgml' PM
// 	11
// 	PM = findexternal(GML)
// 	PM
// 	22
// 	PM->Pgml = NULL
// 	33
// 	rmexternal(GML)
// 	44
}

end
