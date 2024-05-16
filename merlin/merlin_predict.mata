*! version 1.0.0 ?????2016

local gml 	struct merlin_struct scalar
local pgml	pointer(struct merlin_struct scalar) scalar
local RS 	real scalar
local SS 	string scalar
local PS 	pointer scalar
local RM 	real matrix
local SM	string matrix
local PC 	pointer colvector
local SC 	string colvector
local SR	string rowvector
local TR	transmorphic
local RC	real colvector
local RR	real rowvector
local PM	pointer matrix

version 15.1

mata:

void merlin_predict(`SS' object, `SS' newvar, `SS' touse, `SS' stat, `SS' predtype)
{
	`gml' gml
	swap(gml,*findexternal(object))
	
	merlin_predict_setup(gml,stat,touse)
	stand	= st_local("standardise")!=""

	if 	(stat=="reffects") 	pred = merlin_predict_blups(gml)
	else if (stat=="reses") 	pred = merlin_predict_blups(gml,1)
	else {
		pf 	= merlin_p_getpf(gml,stat)
		pred 	= merlin_predict_core(gml,pf,predtype,stand)
	}
	id = st_addvar("double",tokens(newvar))
	if (stand & !gml.issurv[gml.model])	{
		st_store(1,id,pred)
	}
	else	st_store(.,id,touse,pred)
}

void merlin_predict_setup(`gml' gml, `SS' stat, `SS' touse)
{
	devcode7 = st_local("devcode7")
	if (st_local("panel")!="" & devcode7!="bc928r72crncfye98fyqc9r398ry") {
		merlin_error("Nope")
	}

	gml.myb = st_matrix(st_local("best"))	//get parameter estimates
	merlin_xb(gml,gml.myb)			//and fill up (also updates NI points)
	
	k = strtoreal(st_local("outcome"))
	gml.model = gml.modtoind = k		//keep here merlin_xb() changes it

	if (stat!="tprob" & stat!="tlos") {
		merlin_predict_error_check(gml,stat)
	}
	
	if (st_global("e(cmd2)")=="stexcess") {
		gml.Puserst[1,3] = &merlin_p_stexcess_ch()
		gml.Puserst[1,1] = &merlin_p_stexcess_h()
	}
	
	//update N
	if (st_local("npredict")=="") {
		gml.N = gml.Nobs[gml.Nlevels,k]
	}

	//fill in timevar if specified
	if (st_local("timevar")!="") {
		asarray(gml.timevars,gml.model,st_data(.,st_local("timevar"),touse))
	}
	else {
// 		if (st_global("e(timevar"+st_local("outcome")+")")=="") {
// 			if (gml.familys[k]!="null") {
// 				//timevar in struct gets called, so needs to be outcome for survival models
// 				asarray(gml.timevars,gml.model,asarray(gml.y,gml.model)[,1])
// 			}
// 			else asarray(gml.timevars,gml.model,J(gml.N,1,0))	//fudge
// 		}
	}
	
}

void merlin_predict_error_check(`gml' gml, `SS' stat)
{
	f = gml.familys[gml.model]
	
	issurv = f=="addhazard" | f=="pwexponential" | f=="loglogistic" | 
                f=="lognormal" | f=="cox" | f=="rp" | f=="aft" | 
                f=="exponential" | f=="weibull" | f=="gompertz" | 
                f=="ggamma" | f=="loghazard" | f=="logchazard" | 
                f=="plogchazard" | 
                (f=="user" & st_global("e(failure"+strofreal(gml.model)+")")!="")
	
	if (!issurv) {
		if (stat=="hazard" | stat=="survival" | stat=="chazard" | 
                        stat=="logchazard" | stat=="cif" | stat=="rmst" | 
                        stat=="timelost") {
			merlin_error(stat+" not allowed with family("+f+")")
		}
	}
	
	if (stat=="mu") {
		if (f=="addhazard" | f=="loglogistic" | f=="lognormal" | 
                        f=="cox" | f=="rp" | f=="exponential" | f=="weibull" | 
                        f=="gompertz" | f=="ggamma" | f=="lquantile" | 
                        f=="loghazard" | f=="logchazard" | f=="user" | 
                        f=="ordinal") {
			merlin_error("mu not supported with family("+f+")")
		}
	}
	
	if (stat!="eta" & 
                st_global("e(llfunction"+strofreal(gml.model)+")")!="" &
		st_global("e(cmd2)")!="stexcess") {
		merlin_error("eta only valid prediction with family(user, llfunction())")
	}
	
}

`PS' merlin_p_getpf(`gml' gml, `SS' stat)
{

	if 	(stat=="eta")		return(&merlin_p_eta())
	else if (stat=="cif") 		return(&merlin_p_cif())
	else if (stat=="totalsurvival") return(&merlin_p_totalsurvival())
	else if (stat=="rmst")		return(&merlin_p_rmst())
	else if (stat=="timelost")	return(&merlin_p_timelost())
	else if (stat=="totaltimelost")	return(&merlin_p_totaltimelost())
	else if (stat=="tprob")		return(&merlin_p_transprob())
	else if (stat=="tlos")		return(&merlin_p_los())
	else if (stat=="user") {
		stata("mata: _merlin_userf = &"+st_local("userfunction")+"()")
		external _merlin_userf
		return(_merlin_userf)
	}
	else {
		f = gml.familys[gml.model]
		if 	(f=="exponential") {
			if 	(stat=="hazard")	return(&merlin_p_exp_h())
			else if (stat=="survival")	return(&merlin_p_exp_s())
			else if (stat=="chazard")	return(&merlin_p_exp_ch())
			else if (stat=="logchazard")	return(&merlin_p_exp_logch())
			else if (stat=="mu") 		return(&merlin_p_exp_mu())
                        else if (stat=="density") 	return(&merlin_p_exp_dens())
		}
		else if (f=="weibull") {
			if 	(stat=="hazard")	return(&merlin_p_weibull_h())
			else if (stat=="survival")	return(&merlin_p_weibull_s())
			else if (stat=="chazard")	return(&merlin_p_weibull_ch())
			else if (stat=="logchazard")	return(&merlin_p_weibull_logch())
			else if (stat=="mu") 		return(&merlin_p_weibull_mu())
                        else if (stat=="density") 	return(&merlin_p_weibull_dens())
		}
		else if (f=="pwexponential") {
			if 	(stat=="hazard")	return(&merlin_p_pwexp_h())
			else if (stat=="survival")	return(&merlin_p_pwexp_s())
			else if (stat=="chazard")	return(&merlin_p_pwexp_ch())
			else if (stat=="logchazard")	return(&merlin_p_pwexp_logch())
			else if (stat=="mu") 		return(&merlin_p_pwexp_mu())
                        else if (stat=="density") 	return(&merlin_p_pwexp_dens())
		}
		else if (f=="gompertz") {
			if 	(stat=="hazard")	return(&merlin_p_gompertz_h())
			else if (stat=="survival")	return(&merlin_p_gompertz_s())
			else if (stat=="chazard")	return(&merlin_p_gompertz_ch())
			else if (stat=="logchazard")	return(&merlin_p_gompertz_logch())
                        else if (stat=="density")	return(&merlin_p_gompertz_dens())
		}
		else if (f=="ggamma") {
			if 	(stat=="hazard")	return(&merlin_p_ggamma_h())
			else if (stat=="survival")	return(&merlin_p_ggamma_s())
			else if (stat=="chazard")	return(&merlin_p_ggamma_ch())
			else if (stat=="logchazard")	return(&merlin_p_ggamma_logch())
                        else if (stat=="density")	return(&merlin_p_ggamma_pdf())
		}
		else if (f=="lognormal") {
			if 	(stat=="hazard")	return(&merlin_p_lognormal_h())
			else if (stat=="survival")	return(&merlin_p_lognormal_s())
			else if (stat=="chazard")	return(&merlin_p_lognormal_ch())
			else if (stat=="logchazard")	return(&merlin_p_lognormal_logch())
                        else if (stat=="density")	return(&merlin_p_lognormal_pdf())
		}
		else if (f=="loglogistic") {
			if 	(stat=="hazard")	return(&merlin_p_loglogistic_h())
			else if (stat=="survival")	return(&merlin_p_loglogistic_s())
			else if (stat=="chazard")	return(&merlin_p_loglogistic_ch())
			else if (stat=="logchazard")	return(&merlin_p_loglogistic_logch())
                        else if (stat=="density")	return(&merlin_p_loglogistic_pdf())
		}
		else if (f=="loghazard") {
			if 	(stat=="hazard")	return(&merlin_p_loghazard_h())
			else if (stat=="survival")	return(&merlin_p_loghazard_s())
			else if (stat=="chazard")	return(&merlin_p_loghazard_ch())
			else if (stat=="logchazard")	return(&merlin_p_loghazard_logch())
                        else if (stat=="density")	return(&merlin_p_loghazard_dens())
		}
		else if (f=="addhazard") {
			if 	(stat=="hazard")	return(&merlin_p_addhazard_h())
			else if (stat=="survival")	return(&merlin_p_addhazard_s())
			else if (stat=="chazard")	return(&merlin_p_addhazard_ch())
			else if (stat=="logchazard")	return(&merlin_p_addhazard_logch())
                        else if (stat=="density")	return(&merlin_p_addhazard_dens())
		}
		else if (f=="logchazard" | f=="plogchazard") {
			if 	(stat=="hazard")	return(&merlin_p_logchazard_h())
			else if (stat=="survival")	return(&merlin_p_logchazard_s())
			else if (stat=="chazard")	return(&merlin_p_logchazard_ch())
			else if (stat=="logchazard")	return(&merlin_p_logchazard_logch())
                        else if (stat=="density")	return(&merlin_p_logchazard_dens())
		}
		else if (f=="rp") {
			if 	(stat=="hazard")	return(&merlin_p_rp_h())
			else if (stat=="survival")	return(&merlin_p_rp_s())
			else if (stat=="chazard")	return(&merlin_p_rp_ch())
			else if (stat=="logchazard")	return(&merlin_p_rp_logch())
                        else if (stat=="density")	return(&merlin_p_rp_dens())
		}
		else if (f=="aft") {
			if 	(stat=="hazard")	return(&merlin_p_aft_h())
			else if (stat=="survival")	return(&merlin_p_aft_s())
			else if (stat=="chazard")	return(&merlin_p_aft_ch())
			else if (stat=="logchazard")	return(&merlin_p_aft_logch())
                        else if (stat=="density")	return(&merlin_p_aft_dens())
		}
		else if (f=="cox") {
			if 	(stat=="basehazard")	return(&merlin_p_cox_h0())
			else if (stat=="hazard")        return(&merlin_p_cox_h())
			else if (stat=="survival")	return(&merlin_p_cox_s())
			else if (stat=="chazard")	return(&merlin_p_cox_ch())
			else if (stat=="logchazard")	return(&merlin_p_cox_logch())
                        else if (stat=="density")	return(&merlin_p_cox_dens())
		}
		else if (f=="bernoulli") {
			if (stat=="mu") return(&merlin_p_bernoulli_mu())
		}
		else if (f=="beta") {
			if (stat=="mu") return(&merlin_p_beta_mu())
		}
		else if (f=="gamma") {
			if (stat=="mu") return(&merlin_p_gamma_mu())
		}
		else if (f=="gaussian") {
			if (stat=="mu") return(&merlin_p_gaussian_mu())
		}
		else if (f=="negbinomial") {
			if (stat=="mu") return(&merlin_p_negbinomial_mu())
		}
		else if (f=="user") {
			strk 	= strofreal(gml.model)
			hash	= st_global("e(hfunction"+strk+")")!=""
			hasch 	= st_global("e(chfunction"+strk+")")!=""
			
			if (hasch) {
				if (hash) {
					if (stat=="hazard")		return(&merlin_p_userhch_h())
					else if (stat=="survival")	return(&merlin_p_userhch_s())
					else if (stat=="chazard")	return(&merlin_p_userhch_ch())
					else if (stat=="logchazard")	return(&merlin_p_userhch_logch())
                                        else if (stat=="density")	return(&merlin_p_userhch_dens())
				}
				else {
					if 	(stat=="hazard")	return(&merlin_p_userloghch_h())
					else if (stat=="survival")	return(&merlin_p_userloghch_s())
					else if (stat=="chazard")	return(&merlin_p_userloghch_ch())
					else if (stat=="logchazard")	return(&merlin_p_userloghch_logch())
                                        else if (stat=="density")	return(&merlin_p_userloghch_dens())
				}
			}
			else {
				if (hash) {
					if 	(stat=="hazard")	return(&merlin_p_userh_h())
					else if (stat=="survival")	return(&merlin_p_userh_s())
					else if (stat=="chazard")	return(&merlin_p_userh_ch())	
					else if (stat=="logchazard")	return(&merlin_p_userh_logch())
                                        else if (stat=="density")	return(&merlin_p_userh_dens())
				}
				else {
					if      (stat=="hazard")	return(&merlin_p_userlogh_h())
					else if (stat=="survival")	return(&merlin_p_userlogh_s())
					else if (stat=="chazard")	return(&merlin_p_userlogh_ch())
					else if (stat=="logchazard")	return(&merlin_p_userlogh_logch())
                                        else if (stat=="density")	return(&merlin_p_userlogh_dens())
				}
			}
			
		}
		
		
	}
}

`RC' merlin_predict_core(`gml' gml, `PS' pf, `SS' predtype, `RS' stand)
{
	gml.fixedonly = 1
	if (predtype=="fitted") {
		if (st_local("panel")!="") {
			merlin_predict_getspecblups(gml)
			gml.fixedonly = 3
		}
		else {
			merlin_predict_getblups(gml)	//fills up all blups 
			gml.fixedonly = 2
		}
		//reset
		gml.model = gml.modtoind = strtoreal(st_local("outcome"))
	}

	gml.survind = 0
	if (st_local("overoutcome")!="") {
		
		pfo = merlin_predict_overoutcome_setup(gml)
		
		if (gml.issurv[gml.model] & stand) {
		
			t 	= asarray(gml.timevars,gml.model)
			Nt 	= rows(t)
			Nobs	= gml.Nobs[gml.Nlevels,gml.model]
			pred 	= J(Nt,1,.)

			for (i=1;i<=Nt;i++) {
				asarray(gml.timevars,gml.model,J(Nobs,1,t[i]))
				pred[i] = mean(merlin_predict_overoutcome(gml,pfo,pf))
			}

		}
		else {
			pred = merlin_predict_overoutcome(gml,pfo,pf)
			if (stand) pred = mean(pred)
		}
		
	}
	else {
		if (gml.issurv[gml.model] & stand) {
		
			t 		= asarray(gml.timevars,gml.model)
			Nt 		= rows(t)
			Nobs	= gml.Nobs[gml.Nlevels,gml.model]
			
			
			if (gml.familys[1,1]=="cox") {
				pred = merlin_p_cox_s_stand(gml)
			}
			else {
				pred 	= J(Nt,1,.)
				for (i=1;i<=Nt;i++) {
					asarray(gml.timevars,gml.model,J(Nobs,1,t[i]))
					pred[i] = mean((*pf)(gml))
				}
			}
			
		}
		else {
			if (predtype=="marginal") {
                                pred = merlin_predict_marginal(1,gml,pf)
                        }
                        else    pred = (*pf)(gml)
			if (stand) pred = mean(pred)
		}
		
	}
	return(pred)
}

//fixed portion only
`RC' merlin_predict_fixedonly(	`gml' gml,	///
				`PS' pf)	//	
{
	gml.fixedonly = 1
	return((*pf)(gml))
}

//include EB means
`RC' merlin_predict_fitted(	`gml' gml,	///
				`PS' pf)	//	
{
	gml.fixedonly = 2
	return((*pf)(gml))
}

//single subject specific
`RC' merlin_predict_ebspec(	`gml' gml,	///
                                `PS' pf)	//
{
	gml.fixedonly = 3
	return((*pf)(gml))
}

//integrate over random effects
`RC' merlin_predict_marginal(	`RS' index,	///	-level-
                                `gml' gml,	///
                                `PS' pf)	//	
{
	gml.fixedonly = 0
        index2 = index + 1
	if (index<gml.Nrelevels) {
		result = J(gml.Nobs[gml.Nlevels,1],gml.ndim[index],0)
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			result[,q] = merlin_predict_marginal(index2,gml,pf)
		}
	}
	else result = (*pf)(gml)

        //GHQ or MCI
	if (gml.usegh[index]) {
                return(result * asarray(gml.baseGHweights,index))
        }
	else 	return(quadrowsum(result):/gml.ndim[index])
}

`RM' merlin_p_eta(`gml' gml)
{
	if (gml.istimedep[gml.model,1]) {
		t = merlin_util_timevar(gml)
		return(merlin_util_xzb(gml,t))
	}
	else return(merlin_util_xzb(gml))
}

`RM' merlin_p_cif(`gml' gml, | `RC' t)
{
	
	if (gml.familys[gml.model]=="cox") {
		return(merlin_p_cox_cif(gml))
	}
	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	refmod 	= gml.model

	//get pointers to cause-specific h and ch functions -> gml.model indexes this
	hf 		= merlin_p_getpf(gml, "hazard")
	chfs 	= J(1,0,NULL)
		
	if (st_local("causes")=="") {
		modind = J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="addhazard" | f=="pwexponential" | f=="loglogistic" | f=="lognormal" | f=="rp" | f=="aft" | f=="exponential" | f=="weibull" | f=="gompertz" | f=="ggamma" | f=="loghazard" | f=="logchazard" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (issurv) {
				chfs 	= chfs,merlin_p_getpf(gml, "chazard")
				modind 	= modind,gml.model
			}
		}
		Nsurvmodels = cols(modind)
	}
	else {
		modind = strtoreal(tokens(st_local("causes")))
		Nsurvmodels = cols(modind)
		for (mod=1; mod<=Nsurvmodels; mod++) {
			modm = modind[mod]
			gml.model = modm
			f = gml.familys[modm]
			issurv = f=="addhazard" | f=="pwexponential" | f=="loglogistic" | f=="lognormal" | f=="rp" | f=="aft" | f=="exponential" | f=="weibull" | f=="gompertz" | f=="ggamma" | f=="loghazard" | f=="logchazard" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (!issurv) merlin_error("model "+strofreal(modm)+" in causes() is not a survival model")
			chfs 	= chfs,merlin_p_getpf(gml, "chazard")
		}
	}
	
	if (Nsurvmodels==1) {
		
		gml.model = modind
		pf = merlin_p_getpf(gml, "survival")
		return(1:-(*pf)(gml,t))
		
	}
	else {
	
		gml.model 	= refmod
		Ngq 		= 30
		gq 			= merlin_gq(Ngq,"legendre")
		result 		= J(gml.N,1,0)
		qp			= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2

		for (q=1; q<=Ngq; q++) {						//cif integral
			gml.model 	= refmod
			haz			= (*hf)(gml,qp[,q])
			ochres 		= J(gml.N,1,0)

			for (k=1; k<=Nsurvmodels; k++) {			//overall cumulative hazard integral
				gml.model = modind[k]
				ochres = ochres :+ (*chfs[k])(gml,qp[,q])
			}
			result = result :+ haz :* exp(-ochres) :* gq[q,2] :* t :/ 2
		}

		gml.model 	= refmod
		return(result)	
	
	}	
	
}

`RM' merlin_p_totalsurvival(`gml' gml, | `RC' t)
{
	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	refmod 	= gml.model

	//get pointers to cause-specific h and ch functions -> gml.model indexes this
	chfs 	= J(1,0,NULL)
		
	if (st_local("causes")=="") {
		modind = J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="addhazard" | f=="pwexponential" | f=="loglogistic" | f=="lognormal" | f=="rp" | f=="aft" | f=="exponential" | f=="weibull" | f=="gompertz" | f=="ggamma" | f=="loghazard" | f=="logchazard" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (issurv) {
				chfs 	= chfs,merlin_p_getpf(gml, "chazard")
				modind 	= modind,gml.model
			}
		}
		Nsurvmodels = cols(modind)
	}
	else {
		modind = strtoreal(tokens(st_local("causes")))
		Nsurvmodels = cols(modind)
		for (mod=1; mod<=Nsurvmodels; mod++) {
			modm = modind[mod]
			gml.model = modm
			f = gml.familys[modm]
			issurv = f=="addhazard" | f=="pwexponential" | f=="loglogistic" | f=="lognormal" | f=="rp" | f=="aft" | f=="exponential" | f=="weibull" | f=="gompertz" | f=="ggamma" | f=="loghazard" | f=="logchazard" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (!issurv) merlin_error("model "+strofreal(modm)+" in causes() is not a survival model")
			chfs 	= chfs,merlin_p_getpf(gml, "chazard")
		}
	}
	
	if (Nsurvmodels==1) {
		
		gml.model = modind
		pf = merlin_p_getpf(gml, "survival")
		return(1:-(*pf)(gml,t))
		
	}
	else {
	
		gml.model 	= refmod
		result 		= J(gml.N,1,0)

		ochres 		= J(gml.N,1,0)
		for (k=1; k<=Nsurvmodels; k++) {			//overall cumulative hazard integral
			gml.model = modind[k]
			result = result :+ (*chfs[k])(gml,t)
		}
		gml.model = refmod
		return(exp(-result))	
	}	
	
}

`RM' merlin_p_timelost(`gml' gml, | `RC' t)
{
	if (gml.familys[gml.model]=="cox") {
		return(merlin_p_cox_timelost(gml))
	}
	
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	Ngq 	= 30
	gq 		= merlin_gq(Ngq,"legendre")
	qp		= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	
	//integrate cif
	result 	= J(gml.N,1,0)
	for (q=1; q<=Ngq; q++) {						
		result = result :+ merlin_p_cif(gml,qp[,q]) :* gq[q,2] :* t :/ 2
	}
	return(result)	
}

`RM' merlin_p_totaltimelost(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	refmod 	= gml.model
	
	if (st_local("causes")=="") {
		modind 	= J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="addhazard" | f=="pwexponential" | f=="loglogistic" | f=="lognormal" | f=="rp" | f=="aft" | f=="exponential" | f=="weibull" | f=="gompertz" | f=="ggamma" | f=="loghazard" | f=="logchazard" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (issurv) {
				modind 	= modind,gml.model
			}
		}
		Nsurvmodels = cols(modind)
	}
	else {
		modind = strtoreal(tokens(st_local("causes")))
		Nsurvmodels = cols(modind)
		for (mod=1; mod<=Nsurvmodels; mod++) {
			modm = modind[mod]
			gml.model = modm
			f = gml.familys[modm]
			issurv = f=="addhazard" | f=="pwexponential" | f=="loglogistic" | f=="lognormal" | f=="rp" | f=="aft" | f=="exponential" | f=="weibull" | f=="gompertz" | f=="ggamma" | f=="loghazard" | f=="logchazard" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (!issurv) merlin_error("model "+strofreal(modm)+" in causes() is not a survival model")
			modind 	= modind,modm
		}
		
	}
	
	//total time lost
	result 	= J(gml.N,1,0)
	for (q=1; q<=Nsurvmodels; q++) {		
		gml.model 	= modind[q]
		result 		= result :+ merlin_p_timelost(gml,t)
	}
	return(result)	
}

`RM' merlin_p_rmst(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	return(t:-merlin_p_totaltimelost(gml,t))	
}

`PS' merlin_predict_overoutcome_setup(`gml' gml)
{
	overmod = strtoreal(st_local("overmodel"))
	f	 	= gml.familys[overmod]
		
	if 	(f=="bernoulli") return(&merlin_p_over_bernoulli_mu())
	else if (f=="gaussian")  return(&merlin_p_over_gaussian_mu())
}

`RC' merlin_predict_overoutcome(`gml' gml, `PS' pfo, `PS' pf)
{
	return((*pfo)(gml,pf))
}

/*
- calculate blups for all panels and store within gml.blups array
*/

`RM' merlin_predict_blups(`gml' gml, | `RS' getses)
{
	if (args()==2) 	merlin_predict_getblups(gml,getses)
	else 			merlin_predict_getblups(gml)
	gml.survind = 0
	gml.model = gml.modtoind = strtoreal(st_local("outcome"))	//getblups changes it
	return(asarray(gml.blups,1)[merlin_get_adpanelindex(gml,1),])
}

void merlin_predict_getblups(`gml' gml, | `RS' getses)
{
	getses 		= args()==2
	gml.blups 	= asarray_create("real",1)
	
	//adaptive updates to calculate blups
	gml.fixedonly = 0
	oldlnl 	= quadsum(merlin_logl_panels(1,gml),1)

	(*gml.Pupdateip)(gml)
	newlnl = quadsum(merlin_logl_panels(1,gml),1)

	while (reldif(oldlnl,newlnl)>gml.atol) {
		swap(oldlnl,newlnl)
		(*gml.Pupdateip)(gml)
		newlnl = quadsum(merlin_logl_panels(1,gml),1)
	}	
	
	//now store them -> unnecessary extra step, but o/w would have to change updateip function
	if (getses) {
		for (i=1; i<=1; i++) {
			seblups 	= J(gml.Nobs[i,1],gml.Nres[i],.)
			baseweights 	= asarray(gml.baseGHweights,i)
			L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
			numer 		= asarray(gml.Li_ip,gml.qind) :/ L_i
			baseweights2 	= baseweights'	
			cholV		= cholesky(asarray(gml.vcvs,i))
			for (j=1; j<=gml.Nobs[i,1]; j++) {
				ipij 		= cholV * asarray(gml.aghip,(i,j))
				blups 		= numer[j,] * (ipij :* baseweights2)'
				vcv_new 	= (numer[j,] * (merlin_outerprod_by_col(ipij') :* baseweights)) :- blups # blups
				seblups[j,] = sqrt(diagonal(rowshape(vcv_new,gml.Nres[i])))'
			}
			asarray(gml.blups,i,seblups)
		}
	}
	else {
		for (i=1; i<=1; i++) {
			blups 		= J(gml.Nobs[i,1],gml.Nres[i],.)
			baseweights 	= asarray(gml.baseGHweights,i)
			L_i 		= asarray(gml.Li_ip,gml.qind) * baseweights
			numer 		= asarray(gml.Li_ip,gml.qind) :/ L_i
			baseweights2 	= baseweights'	
			cholV		= cholesky(asarray(gml.vcvs,i))
			for (j=1; j<=gml.Nobs[i,1]; j++) {
				ipij 		= cholV * asarray(gml.aghip,(i,j))
				blups[j,] 	= numer[j,] * (ipij :* baseweights2)'
			}
			asarray(gml.blups,i,blups)
		}
	}
	
}

`RM' merlin_predict_getspecblups(`gml' gml)
{
	gml.blups = asarray_create("real",1)
	
	stata("tempvar specb spectouse")
	if (st_local("blupif")!="") {
		stata("qui gen byte "+st_local("spectouse")+ "= "+st_global("e(levelvars)")+"=="+st_local("panel")+" & "+st_local("touse")+" & "+st_local("blupif"))
	}
	else {
		stata("qui gen byte "+st_local("spectouse")+ "= "+st_global("e(levelvars)")+"=="+st_local("panel")+" & "+st_local("touse"))
	}
	stata("qui predict "+st_local("specb")+"* if "+st_local("spectouse")+", reffects")
	specblups = st_data(.,st_local("specb")+"*",st_local("spectouse"))
	asarray(gml.blups,1,specblups[1,])
}


end
