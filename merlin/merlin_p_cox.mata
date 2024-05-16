*! version 1.0.0 ?????2016

local GMLS 		struct merlin_struct scalar
local pGMLS		pointer(struct merlin_struct scalar) scalar
local TR 		transmorphic
local RS 		real scalar
local RC 		real colvector
local SS 		string scalar
local PS 		pointer scalar
local RR 		real rowvector
local RM 		real matrix
local PC 		pointer colvector
local PM 		pointer matrix
local SC 		string colvector

version 14.1

mata:

`RM' merlin_p_cox_h0(`GMLS' gml)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	haslt	= gml.hasltrunc[model]
	haz 	= J(merlin_get_nobs(gml,model),1,0)

	//exactly observed events
	gml.survind = 1
	
	//get index for sorted failure times
	Nfails = gml.Nsurv[model,1]
	index1 = asarray(gml.surv_index,(model,1))
	
	//xb at all survival/censoring times
	gml.survind = 0

	if (!gml.istimedep[model,1]) {	//time-independent
		expxb2 	= exp(merlin_util_xzb(gml,y[,1]))
		if (haslt) {
			for (i=1;i<=Nfails;i++) {
				atriskid = selectindex((y[index1[i],1] :<= y[,1]) :& (y[index1[i],1] :> y[,3]))
				D = sum(y[index1[i],1] :== select(y[,1],y[,2]:==1))
				haz[index1[i],] = D:/quadcolsum(expxb2[atriskid,])
			}
		}
		else {
			for (i=1;i<=Nfails;i++) {
				atriskid = selectindex(y[index1[i],1] :<= y[,1])
				D = sum(y[index1[i],1] :== select(y[,1],y[,2]:==1))
				haz[index1[i],] = D:/quadcolsum(expxb2[atriskid,])
			}
		}
	}
	else {				//time-dependent
		if (haslt) {
			for (i=1;i<=Nfails;i++) {
				atriskid = selectindex((y[index1[i],1] :<= y[,1]) :& (y[index1[i],1] :> y[,3]))
				expxb2   = exp(merlin_util_xzb(gml,J(merlin_get_nobs(gml,model),1,y[index1[i],1])))
				D	 = sum(y[index1[i],1] :== select(y[,1],y[,2]:==1))
				haz[index1[i],] = D:/quadcolsum(expxb2[atriskid,])
			}
		}
		else {
			for (i=1;i<=Nfails;i++) {
				atriskid = selectindex(y[index1[i],1] :<= y[,1])
				expxb2 	 = exp(merlin_util_xzb(gml,J(merlin_get_nobs(gml,model),1,y[index1[i],1])))
				D	 = sum(y[index1[i],1] :== select(y[,1],y[,2]:==1))
				haz[index1[i],] = D:/quadcolsum(expxb2[atriskid,])
			}
		}
	}
	
	return(haz)
}

`RM' merlin_p_cox_h(`GMLS' gml)
{		
	h0 = st_data(.,st_local("baseh"+strofreal(gml.model)),st_local("touse"))
	gml.survind = 0
	return(h0 :* exp(merlin_util_xzb(gml)))
}

`RM' merlin_p_cox_ch0(`GMLS' gml)
{	
	h0 = st_data(.,st_local("baseh"+strofreal(gml.model)),st_local("touse"))	
	N  = gml.Nobs[gml.Nlevels,gml.model]
	y  = merlin_util_depvar(gml)
	
	huniq 	= h0,y[,1]
	huniq 	= uniqrows(select(huniq,y[,2]))
	_sort(huniq,2)
	chuniq 	= runningsum(huniq[,1])

	ch = J(N,1,.)
	for (i=1;i<=N;i++) {
		ind = max(selectindex(y[i,1]:>=huniq[,2]))
		if (ind==.) ch[i] = 0	//censoring prior to first event
		else ch[i] = chuniq[ind]
	}

	return(ch)
}

`RM' merlin_p_cox_ch(`GMLS' gml)
{
	gml.survind = 0
	if (!gml.istimedep[gml.model,1]) {
		return(merlin_p_cox_ch0(gml) :* exp(merlin_util_xzb(gml)))
	}
	else {
				
		h0 = st_data(.,st_local("baseh"+strofreal(gml.model)),st_local("touse"))
		y  = merlin_util_depvar(gml)
		N  = gml.Nobs[gml.Nlevels,gml.model]
		
		huniq 	= h0,y[,1]
		huniq 	= uniqrows(select(huniq,y[,2]))
		_sort(huniq,2)
		
		Nuniq	= rows(huniq)
		base	= J(N,Nuniq,.)

		for (i=1;i<=Nuniq;i++) {
			asarray(gml.timevars,gml.model,J(N,1,huniq[i,2]))
			base[,i] = merlin_util_xzb(gml)
		}
		
		haz = (huniq[,1]' :* exp(base))'
		for (i=1;i<=N;i++) {
			haz[,i] = runningsum(haz[,i])
		}

		ch = J(N,1,.)
		for (i=1;i<=N;i++) {
			index = max(selectindex(y[i,1]:>=huniq[,2]))
			ch[i] = haz[index,i]
		}

		return(ch)
		
	}
	
}

`RM' merlin_p_cox_logch(`GMLS' gml)
{	
	return(log(merlin_p_cox_ch(gml)))
}

`RM' merlin_p_cox_s(`GMLS' gml)
{
	return(exp(-merlin_p_cox_ch(gml)))
}

`RM' merlin_p_cox_f(`GMLS' gml)
{
	return(1:-merlin_p_cox_s(gml))
}

`RM' merlin_p_cox_dens(`GMLS' gml)
{
	return(merlin_p_cox_h(gml):*merlin_p_cox_s(gml))
}

`RM' merlin_p_cox_s_stand(`GMLS' gml)
{	
	gml.survind = 0
	if (!gml.istimedep[gml.model,1]) {
		return(mean(exp(-exp(merlin_util_xzb(gml)) * merlin_p_cox_ch0(gml)'))')
	}
	else {
		
		h0 = st_data(.,st_local("baseh"+strofreal(gml.model)),st_local("touse"))
		y  = merlin_util_depvar(gml)
		N  = gml.Nobs[gml.Nlevels,gml.model]
		
		huniq 	= h0,y[,1]
		huniq 	= uniqrows(select(huniq,y[,2]))
		_sort(huniq,2)
		
		Nuniq	= rows(huniq)
		base	= J(N,Nuniq,.)

		for (i=1;i<=Nuniq;i++) {
			asarray(gml.timevars,gml.model,J(N,1,huniq[i,2]))
			base[,i] = merlin_util_xzb(gml)
		}
		
		pred = (huniq[,1]' :* exp(base))'
		for (i=1;i<=N;i++) {
			pred[,i] = exp(-runningsum(pred[,i]))
		}
		S = mean(pred')'
		index = J(N,1,.)
		for (i=1;i<=N;i++) {
			index[i] = max(selectindex(y[i,1]:>=huniq[,2]))
		}
		return(S[index])
		
	}
	
}

`RM' merlin_p_cox_cif(`GMLS' gml)
{
	Nobs	= gml.N
	refmod 	= gml.model
	
	if (st_local("causes")=="") {
		modind = J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="cox"
			if (issurv) modind = modind,gml.model
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
			issurv = f=="cox" 
			if (!issurv) merlin_error("model "+strofreal(modm)+" in causes() is not a Cox model")
		}

	}
	
	if (Nsurvmodels==1) {
		
		gml.model = modind
		return(merlin_p_cox_f(gml))
		
	}
	else {
	
		gml.model 	= refmod
		result 		= J(Nobs,1,0)
		
		bhvars = J(1,0,"")
		for (i=1;i<=Nsurvmodels;i++) bhvars = bhvars,st_local("baseh"+strofreal(modind[i]))
		
		h0 = st_data(.,bhvars,st_local("touse"))

		//need outcome() hazard at unique event times, for all obs
		ymain		= merlin_util_depvar(gml)
		N			= gml.Nobs[gml.Nlevels,gml.model]

		huniqmain	= h0[,refmod],ymain[,1]
		huniqmain 	= uniqrows(select(huniqmain,ymain[,2]))
		_sort(huniqmain,2)
		Nuniqmain	= rows(huniqmain)

		if (!gml.istimedep[gml.model,1]) {
			pred = exp(merlin_util_xzb(gml)) * huniqmain[,1]'
		}
		else {

			pred 	= J(N,Nuniqmain,0)
			t 		= asarray(gml.timevars,gml.model)
			
			for (i=1;i<=Nuniqmain;i++) {				
				asarray(gml.timevars,gml.model,J(N,1,t[i]))
				pred[,i] = merlin_util_xzb(gml)
			}
			pred = exp(pred) :* huniqmain[,1]'
		}

		//now need survival for all outcomes at all outcome times, for all obs

		for (k=1; k<=Nsurvmodels; k++) {
			
			gml.model 	= modind[k]
			y 			= merlin_util_depvar(gml)
			N			= gml.Nobs[gml.Nlevels,gml.model]
			huniq 		= h0[,gml.model],y[,1]
			huniq 		= uniqrows(select(huniq,y[,2]))
			_sort(huniq,2)
			
			if (!gml.istimedep[gml.model,1]) {
				haz = exp(merlin_util_xzb(gml)) * huniq[,1]'
			}
			else {
				Nuniq	= rows(huniq)
				haz 	= J(N,Nuniq,0)
				t 		= asarray(gml.timevars,gml.model)
				
				for (i=1;i<=Nuniq;i++) {				
					asarray(gml.timevars,gml.model,J(N,1,t[i]))
					haz[,i] = merlin_util_xzb(gml)
				}
				haz = exp(haz) :* huniq[,1]'
			}

			_transpose(haz)
			for (i=1;i<=N;i++) {
				haz[,i] = runningsum(haz[,i])
			}
			_transpose(haz)
			
			if (gml.model==refmod) {
				pred = pred :* exp(-haz)
			}
			else {

				for (i=1;i<=Nuniqmain;i++) {
					ind = max(selectindex(huniqmain[i,2]:>=huniq[,2]))
					if (ind!=.) pred[,i] = pred[,i] :* exp(-haz[,ind])
				}

			}
			
		}	

		_transpose(pred)
		for (i=1;i<=N;i++) {
			pred[,i] = runningsum(pred[,i])
		}

		gml.model 	= refmod

		cif = J(N,1,.)
		for (i=1;i<=N;i++) {
			index = max(selectindex(ymain[i,1]:>=huniqmain[,2]))
			if (index==.) 	cif[i] = 0
			else 			cif[i] = pred[index,i]
		}

		return(cif)
		

	}	
	
}

`RM' merlin_p_cox_timelost(`GMLS' gml)
{
	Nobs	= gml.N
	refmod 	= gml.model
	
	if (st_local("causes")=="") {
		modind = J(1,0,.)
		//assume all survival models contribute
		for (mod=1; mod<=gml.Nmodels; mod++) {
			gml.model = mod
			f = gml.familys[mod]
			issurv = f=="cox"
			if (issurv) modind 	= modind,gml.model
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
			issurv = f=="cox" 
			if (!issurv) merlin_error("model "+strofreal(modm)+" in causes() is not a Cox model")
		}
	}
	
	//calc.
	
	gml.model 	= refmod
	result 		= J(Nobs,1,0)
	
	bhvars = J(1,0,"")
	for (i=1;i<=Nsurvmodels;i++) bhvars = bhvars,st_local("baseh"+strofreal(modind[i]))
	
	h0 = st_data(.,bhvars,st_local("touse"))

	//need outcome() hazard at unique event times, for all obs
	ymain		= merlin_util_depvar(gml)
	N			= gml.Nobs[gml.Nlevels,gml.model]

	huniqmain	= h0[,refmod],ymain[,1]
	huniqmain 	= uniqrows(select(huniqmain,ymain[,2]))
	_sort(huniqmain,2)
	Nuniqmain	= rows(huniqmain)

	if (!gml.istimedep[gml.model,1]) {
		pred = exp(merlin_util_xzb(gml)) * huniqmain[,1]'
	}
	else {

		pred 	= J(N,Nuniqmain,0)
		t 		= asarray(gml.timevars,gml.model)
		
		for (i=1;i<=Nuniqmain;i++) {				
			asarray(gml.timevars,gml.model,J(N,1,t[i]))
			pred[,i] = merlin_util_xzb(gml)
		}
		pred = exp(pred) :* huniqmain[,1]'
	}

	//now need survival for all outcomes at all outcome times, for all obs

	for (k=1; k<=Nsurvmodels; k++) {
		
		gml.model 	= modind[k]
		y 			= merlin_util_depvar(gml)
		N			= gml.Nobs[gml.Nlevels,gml.model]
		huniq 		= h0[,gml.model],y[,1]
		huniq 		= uniqrows(select(huniq,y[,2]))
		_sort(huniq,2)
		
		if (!gml.istimedep[gml.model,1]) {
			haz = exp(merlin_util_xzb(gml)) * huniq[,1]'
		}
		else {
			Nuniq	= rows(huniq)
			haz 	= J(N,Nuniq,0)
			t 		= asarray(gml.timevars,gml.model)
			
			for (i=1;i<=Nuniq;i++) {				
				asarray(gml.timevars,gml.model,J(N,1,t[i]))
				haz[,i] = merlin_util_xzb(gml)
			}
			haz = exp(haz) :* huniq[,1]'
		}

		_transpose(haz)
		for (i=1;i<=N;i++) {
			haz[,i] = runningsum(haz[,i])
		}
		_transpose(haz)
		
		if (gml.model==refmod) {
			pred = pred :* exp(-haz)
		}
		else {
			for (i=1;i<=Nuniqmain;i++) {
				ind = max(selectindex(huniqmain[i,2]:>=huniq[,2]))
				if (ind!=.) pred[,i] = pred[,i] :* exp(-haz[,ind])
			}
		}
		
	}	

	//cif
	_transpose(pred)
	for (i=1;i<=N;i++) {
		pred[,i] = runningsum(pred[,i])
	}

	//area under cif
	
	gml.model = refmod

	cif = J(N,1,.)
	for (i=1;i<=N;i++) {
		index = max(selectindex(ymain[i,1]:>=huniqmain[,2]))
		if (index==.) 	cif[i] = 0
		else 	cif[i] = pred[index,i]
	}

	return(cif)
}

end
