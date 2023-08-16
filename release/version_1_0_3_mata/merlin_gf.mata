*! version 1.0.0 MJC

local TS 	transmorphic scalar
local RS 	real scalar
local RC 	real colvector
local RM 	real matrix
local Pgml 	pointer(struct merlin_struct scalar) scalar
local Sgml 	struct merlin_struct scalar

version 12.1
mata:

//gf core function
void merlin_gf(	`TS' M,
				`RS' todo,
				`RR' b,
				`RC' lnfi,
				`RM' G,
				`RM' H)
{
	`Sgml' gml

	gml  		= moptimize_util_userinfo(M,1)
	gml.todo 	= todo

	gml.myb = b
	merlin_xb(gml,b)															//linear predictors

	if (gml.Nrelevels) 	{
		
		lnfi = merlin_logl_panels(1,gml)										//logl
		gml.Pgml->lnfi1 = lnfi

		if (gml.hasanyltrunc) {													//ltruncation
			gml.ltflag = 1														//gets correct random effect nodes
			lnfi = lnfi :- merlin_ltrunc(gml)
		}
		
	}
	else {
		
		lnfi = merlin_logl_ob(gml)												//logl ob level
		
		if (gml.hasanyltrunc) lnfi = lnfi :- merlin_ltrunc_ob(gml)				//ltruncation

	}

	if (todo==0) return
	
}

`RM' merlin_logl(`Sgml' gml)
{	
	return((*gml.Plnl[gml.model])(gml))
}

`RC' merlin_logl_panels(`RS' index,		///	-level-
						`Sgml' gml)		//	
{
	`RS' index2
	`RM' res, panelindex

	index2 = index+1
	
	res = J(gml.Nobs[index,1],gml.ndim[index],0)
	
	if (index<gml.Nrelevels) {
		panelindex = asarray(gml.panelindexes,(index,1))
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			res[,q] = panelsum(merlin_logl_panels(index2,gml),panelindex)
		}
	}
	else {
		for (j=1;j<=gml.Nmodels;j++) {
			gml.model = gml.modtoind = j
			if (gml.NotNull[j,1]) {
				res2 = merlin_logl(gml)
				if (gml.hasweights[index2])	res2 = res2 :* asarray(gml.weights,(index2,j))
				res = res :+ panelsum(res2,asarray(gml.panelindexes,(index,j)))
			}
		}
	}

	if (gml.adapt[index]) res = res :+ asarray(gml.aghlogl,index)
	expres = exp(res)
	if (gml.adapt[index] | gml.todo) asarray(gml.Li_ip,gml.qind,expres) 

	if (gml.usegh[index]) {		//GHQ
		if (gml.hasweights[index]) return(asarray(gml.weights,(index,1)) :* log(expres * asarray(gml.baseGHweights,index)))
		else return(log(expres * asarray(gml.baseGHweights,index)))
	}
	else {						//MCI
		if (gml.hasweights[index]) return(asarray(gml.weights,(index,1)) :* log(quadrowsum(expres):/gml.ndim[index]))
		else return(log(quadrowsum(expres):/gml.ndim[index]))	
	}
}

`RC' merlin_ltrunc(`Sgml' gml)
{
	result = J(gml.Nobs[1,1],1,0)
	
	for (i=1; i<=gml.Nmodels; i++) {
		
		if (gml.hasltrunc[i]) {
			gml.model = gml.modtoind = i
			result = result :+ merlin_ltrunc_panels(1,gml)
		}
	
	}

	return(result)
}

`RC' merlin_ltrunc_panels(	`RS' index,		///	-level-
							`Sgml' gml)		//	
{
	`RS' index2
	`RM' res, panelindex

	index2 = index+1
	
	res = J(gml.Nobs[index,1],gml.ndim[index],0)
	
	if (index<gml.Nrelevels) {
		panelindex = asarray(gml.panelindexes,(index,1))
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			res[,q] = panelsum(merlin_ltrunc_panels(index2,gml),panelindex)
		}
	}
	else {		
		t0 		= merlin_util_depvar(gml)[,3]
		res2 	= -(*gml.ltruncP[gml.model])(gml,t0)
		//hanlde t0 = 0
		t0index = selectindex(t0:==0)
		nt0 = rows(t0index)
		if (nt0) res2[t0index,] = J(nt0,gml.ndim[index],0)
		if (gml.hasweights[index2])	res2 = res2 :* asarray(gml.weights,(index2,gml.model))
		res = panelsum(res2,asarray(gml.panelindexes,(index,gml.model)))	
	}

	if (gml.adapt[index]) res = res :+ asarray(gml.aghlogl_lt,index)
	expres = exp(res)
	if (gml.adapt[index] | gml.todo) asarray(gml.Li_ip_lt,gml.qind,expres) 

	if (gml.usegh[index]) {		//GHQ
		if (gml.hasweights[index]) return(asarray(gml.weights,(index,1)) :* log(expres * asarray(gml.baseGHweights,index)))
		else return(log(expres * asarray(gml.baseGHweights,index)))
	}
	else {						//MCI
		if (gml.hasweights[index]) return(asarray(gml.weights,(index,1)) :* log(quadrowsum(expres):/gml.ndim[index]))
		else return(log(quadrowsum(expres):/gml.ndim[index]))	
	}
}

`RC' merlin_ltrunc_ob(`Sgml' gml)
{
	result = 0
	
	for (i=1; i<=gml.Nmodels; i++) {
		if (gml.hasltrunc[i]) {
			gml.model 	= gml.modtoind = i
			t0 			= merlin_util_depvar(gml)[,3]
			if (gml.hasweights[1]) 	ltres = (*gml.ltruncP[gml.model])(gml,t0) :* asarray(gml.weights,(1,i))
			else 					ltres = (*gml.ltruncP[gml.model])(gml,t0)
			result 	= result :- quadsum(ltres)
		}
	}

	return(result)
}

end
