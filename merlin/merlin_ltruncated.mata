*! version 1.0.0 MJC

local TS 	transmorphic scalar
local RS 	real scalar
local RC 	real colvector
local RM 	real matrix
local Pgml 	pointer(struct merlin_struct scalar) scalar
local Sgml 	struct merlin_struct scalar

version 14.2
mata:

`RC' merlin_ltrunc(`Sgml' gml)
{
	gml.ltflag  = 1		//gets correct random effect nodes
	gml.survind = 7
	gml.model = gml.modtoind = gml.mltmodel
	result = merlin_ltrunc_panels(1,gml)
	gml.survind = 0		//reset
	return(result)
}

`RC' merlin_ltrunc_panels(`RS' index,		///	-level-
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
		t0 	= merlin_util_depvar(gml)[,3]
		res 	= -(*gml.Pch[gml.model])(gml,t0)
//                 rows(res2), cols(res2)
//                 exit(1986)
//		
//                 //create null vector/matrix to index into 
//                
// 		t0index = selectindex(t0:==0)
//                 t0index
// 		nt0 = rows(t0index)
// 		if (nt0) res2[t0index,] = J(nt0,gml.ndim[index],0)
// 		if (gml.hasweights[index2]) {
//                         res2 = res2 :* 
//                                 asarray(gml.weights,(index2,gml.model))
//                 }
// 		res = panelsum(res2,asarray(gml.panelindexes,(index,gml.model)))	
	}

	if (gml.adapt[index]) res = res :+ asarray(gml.aghlogl_lt,index)
	
	expres = exp(res)
	if (gml.adapt[index] | gml.todo) asarray(gml.Li_ip_lt,gml.qind,expres) 

	if (gml.usegh[index]) {		                //GHQ
		if (gml.hasweights[index]) {
                        return(asarray(gml.weights,(index,1)) :* 
                                log(expres * asarray(gml.baseGHweights,index)))
		}
                else return(log(expres * asarray(gml.baseGHweights,index)))
	}
	else {				                //MCI
		if (gml.hasweights[index]) {
                        return(asarray(gml.weights,(index,1)) :* 
                                log(quadrowsum(expres):/gml.ndim[index]))
                }
		else return(log(quadrowsum(expres,1):/gml.ndim[index]))	
	}
}	
	
`RC' merlin_ltrunc_ob(`Sgml' gml)
{
	result = 0
	gml.survind = 4
	for (i=1; i<=gml.Nmodels; i++) {
		if (gml.hasltrunc[i] & gml.Nsurv[i,4]) {
			gml.model 	= gml.modtoind = i
			t0 			= merlin_util_depvar(gml)[,3]
			if (gml.hasweights[1]) 	{
				index4 		= asarray(gml.surv_index,(i,4))
				ltres = (*gml.Pch[i])(gml,t0) :* asarray(gml.weights,(1,i))[index4]
			}
			else ltres = (*gml.Pch[i])(gml,t0)
			result 	= result :- quadsum(ltres,1)
		}
	}
	gml.survind = 0
	return(result)
}

end
