*! version 1.0.0 ?????2016

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

void merlin_predict(`SS' object, `SS' newvar, `SS' touse, `SS' stat, `SS' predtype)
{
	`gml' gml
	gml 		= *findexternal(object)

	merlin_predict_setup(gml,stat,touse)
	pf 		= merlin_p_getpf(gml,stat)
	pred 	= merlin_predict_core(gml,pf,predtype)

	id = st_addvar("double",newvar)
	st_store(.,id,touse,pred)
	
}

void merlin_predict_setup(`gml' gml, `SS' stat, `SS' touse)
{
	gml.myb = st_matrix(st_local("best"))							//get fitted values
	merlin_xb(gml,gml.myb)											//and fill up (also updates NI points)
	
	k			= strtoreal(st_local("outcome"))
	gml.model 	= k													//keep here merlin_xb() changes it

	merlin_predict_error_check(gml,stat)
	
	//update N
	if (st_local("npredict")=="") {
		gml.N = gml.Nobs[gml.Nlevels,k]
	}
	
	//fill in timevar if specified
	if (st_local("timevar")!="") {
		asarray(gml.timevars,gml.model,st_data(.,st_local("timevar"),touse))
	}
	else {
		if (gml.familys[k]!="null") {
			//timevar in struct gets called, so needs to be outcome for survival models
			asarray(gml.timevars,gml.model,asarray(gml.y,gml.model)[,1])
		}
		else asarray(gml.timevars,gml.model,J(gml.N,1,0))	//fudge
	}
	
}

void merlin_predict_error_check(`gml' gml, `SS' stat)
{
	f = gml.familys[gml.model]
	
	issurv = f=="exponential" | f=="weibull" | f=="gompertz" | f=="rcs" | f=="rp" | f=="prp" | (f=="user" & st_global("e(failure"+strofreal(gml.model)+")")!="")
	
	if (!issurv) {
		if (stat=="hazard" | stat=="survival" | stat=="chazard" | stat=="logchazard" | stat=="cif" | stat=="rmst" | stat=="timelost") {
			merlin_error(stat+" not allowed with family("+f+")")
		}
	}
	
	if (stat=="mu") {
		if (f=="gompertz" | f=="lquantile" | f=="rcs" | f=="rp" | f=="user" | f=="ordinal") {
			merlin_error("mu not supported with family("+f+")")
		}
	}
	
	if (stat!="eta" & st_global("e(llfunction"+strofreal(gml.model)+")")!="") {
		merlin_error("eta only valid prediction with family(user, llfunction())")
	}
	
}

`PS' merlin_p_getpf(`gml' gml, `SS' stat)
{
	if 		(stat=="eta")				return(&merlin_p_eta())
	else if (stat=="cif") 				return(&merlin_p_cif())
	else if (stat=="rmst")				return(&merlin_p_rmst())
	else if (stat=="timelost")			return(&merlin_p_timelost())
	else if (stat=="totaltimelost")		return(&merlin_p_totaltimelost())
	else {
		f = gml.familys[gml.model]
		if 		(f=="exponential") {
			if 		(stat=="hazard")		return(&merlin_p_exp_h())
			else if (stat=="survival")		return(&merlin_p_exp_s())
			else if (stat=="chazard")		return(&merlin_p_exp_ch())
			else if (stat=="logchazard")	return(&merlin_p_exp_logch())
			else if (stat=="mu") 			return(&merlin_p_exp_mu())
		}
		else if (f=="weibull") {
			if 		(stat=="hazard")		return(&merlin_p_weibull_h())
			else if (stat=="survival")		return(&merlin_p_weibull_s())
			else if (stat=="chazard")		return(&merlin_p_weibull_ch())
			else if (stat=="logchazard")	return(&merlin_p_weibull_logch())
			else if (stat=="mu") 			return(&merlin_p_weibull_mu())
		}
		else if (f=="gompertz") {
			if 		(stat=="hazard")		return(&merlin_p_gompertz_h())
			else if (stat=="survival")		return(&merlin_p_gompertz_s())
			else if (stat=="chazard")		return(&merlin_p_gompertz_ch())
			else if (stat=="logchazard")	return(&merlin_p_gompertz_logch())
		}
		else if (f=="rcs") {
			if 		(stat=="hazard")		return(&merlin_p_rcs_h())
			else if (stat=="survival")		return(&merlin_p_rcs_s())
			else if (stat=="chazard")		return(&merlin_p_rcs_ch())
			else if (stat=="logchazard")	return(&merlin_p_rcs_logch())
		}
		else if (f=="rp" | f=="prp") {
			if 		(stat=="hazard")		return(&merlin_p_rp_h())
			else if (stat=="survival")		return(&merlin_p_rp_s())
			else if (stat=="chazard")		return(&merlin_p_rp_ch())
			else if (stat=="logchazard")	return(&merlin_p_rp_logch())
		}
		else if (f=="bernoulli") {
			if 		(stat=="mu") 		return(&merlin_p_bernoulli_mu())
		}
		else if (f=="beta") {
			if 		(stat=="mu") 		return(&merlin_p_beta_mu())
		}
		else if (f=="gamma") {
			if 		(stat=="mu") 		return(&merlin_p_gamma_mu())
		}
		else if (f=="gaussian") {
			if 		(stat=="mu") 		return(&merlin_p_gaussian_mu())
		}
		else if (f=="negbinomial") {
			if 		(stat=="mu") 		return(&merlin_p_negbinomial_mu())
		}
		else if (f=="user") {
			strk 	= strofreal(gml.model)
			hash	= st_global("e(hfunction"+strk+")")!=""
			hasch 	= st_global("e(chfunction"+strk+")")!=""
			
			if (hasch) {
				if (hash) {
					if 		(stat=="hazard")		return(&merlin_p_userhch_h())
					else if (stat=="survival")		return(&merlin_p_userhch_s())
					else if (stat=="chazard")		return(&merlin_p_userhch_ch())
					else if (stat=="logchazard")	return(&merlin_p_userhch_logch())
				}
				else {
					if 		(stat=="hazard")		return(&merlin_p_userloghch_h())
					else if (stat=="survival")		return(&merlin_p_userloghch_s())
					else if (stat=="chazard")		return(&merlin_p_userloghch_ch())
					else if (stat=="logchazard")	return(&merlin_p_userloghch_logch())
				}
			}
			else {
				if (hash) {
					if 		(stat=="hazard")		return(&merlin_p_userh_h())
					else if (stat=="survival")		return(&merlin_p_userh_s())
					else if (stat=="chazard")		return(&merlin_p_userh_ch())	
					else if (stat=="logchazard")	return(&merlin_p_userh_logch())
				}
				else {
					if 		(stat=="hazard")		return(&merlin_p_userlogh_h())
					else if (stat=="survival")		return(&merlin_p_userlogh_s())
					else if (stat=="chazard")		return(&merlin_p_userlogh_ch())
					else if (stat=="logchazard")	return(&merlin_p_userlogh_logch())
				}
			}
			
		}
		
		
	}
}

`RC' merlin_predict_core(`gml' gml, `PS' pf, `SS' predtype)
{
	if (predtype=="fixedonly") 	return(merlin_predict_fixedonly(gml,pf))
	else 						return(merlin_predict_marginal(1,gml,pf))
}

`RC' merlin_predict_fixedonly(	`gml' gml,		///
								`PS' pf)		//	
{
	gml.fixedonly = 1
	return((*pf)(gml))
}


`RC' merlin_predict_marginal(	`RS' index,		///	-level-
								`gml' gml,		///
								`PS' pf)		//	
{
	index2 = index+1

	if (index<gml.Nrelevels) {
		result = J(gml.Nobs[gml.Nlevels,1],gml.ndim[index],0)
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			result[,q] = merlin_predict_marginal(index2,gml,pf)
		}
	}
	else result = (*pf)(gml)

	if (gml.usegh[index]) 	return(result * asarray(gml.baseGHweights,index))	//GHQ
	else 					return(quadrowsum(result):/gml.ndim[index])			//MCI
}

`RM' merlin_p_eta(`gml' gml)
{
	t = merlin_util_timevar(gml)
	return(merlin_util_xzb(gml,t))
}

`RM' merlin_p_cif(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)
	
	refmod 	= gml.model
	
	//get pointers to cause-specific h and ch functions -> gml.model indexes this
	hf 		= merlin_p_getpf(gml, "hazard")
	chfs 	= merlin_p_getpf(gml, "chazard")
	
	modind 	= gml.model
		
	if (st_local("causes")=="") {
		
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="exponential" | f=="weibull" | f=="gompertz" | f=="rcs" | f=="rp" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (issurv & mod!=refmod) {
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
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="exponential" | f=="weibull" | f=="gompertz" | f=="rcs" | f=="rp" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (!issurv) merlin_error("model "+strofreal(mod)+" in cause() is not a survival model")
			if (mod!=refmod) {
				chfs 	= chfs,merlin_p_getpf(gml, "chazard")
				modind 	= modind,gml.model
			}
		}
		
	}
	
	gml.model 	= refmod

	Ngq 	= 30
	gq 		= merlin_gq(Ngq,"legendre")
	result 	= J(gml.N,1,0)
	
	qp		= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	
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

`RM' merlin_p_timelost(`gml' gml, | `RC' t)
{
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
	modind 	= gml.model
		
	if (st_local("causes")=="") {
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="exponential" | f=="weibull" | f=="gompertz" | f=="rcs" | f=="rp" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (issurv & mod!=refmod) {
				modind 	= modind,gml.model
			}
		}
		Nsurvmodels = cols(modind)
	}
	else {
		modind = strtoreal(tokens(st_local("causes")))
		Nsurvmodels = cols(modind)
		for (mod=1; mod<=Nsurvmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="exponential" | f=="weibull" | f=="gompertz" | f=="rcs" | f=="rp" | (f=="user" & st_global("e(llfunction"+strofreal(mod)+")")=="")
			if (!issurv) merlin_error("model "+strofreal(mod)+" in cause() is not a survival model")
			if (mod!=refmod) {
				modind 	= modind,gml.model
			}
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

end
