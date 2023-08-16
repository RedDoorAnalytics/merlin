*! version 1.0.0 MJC

local TS 	transmorphic scalar
local RS 	real scalar
local RC 	real colvector
local RM 	real matrix
local Pgml 	pointer(struct merlin_struct scalar) scalar
local Sgml 	struct merlin_struct scalar

version 14.2
mata:

`RM' merlin_logl_survival_ob(`Sgml' gml, | `RM' G, `RM' H)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	hasbh	= gml.hasbh[model,]
	logl 	= gml.Nrelevels ? J(gml.Nobs[gml.Nlevels,model],
                        gml.ndim[gml.Nrelevels],0) : 
                        J(merlin_get_nobs(gml,model),1,0)

	//exactly observed events -> hazard function
	gml.survind = 1
	if (merlin_get_nobs(gml,model)) {
		index1 = merlin_get_surv_index(gml)
		if (gml.hastmat) {
                        logl[index1,] 	= merlin_logl_ms(gml,y[index1,1])
                }
		else {
			logl[index1,] = (*gml.Plogh[model])(gml,y[index1,1])
			if (hasbh[1]) {
				logl[index1,] = log(exp(logl[index1,]) :+ 
					merlin_util_bhazard(gml))
			}
			else if (hasbh[2]) {		
				logl[index1,] = log(exp(logl[index1,]) :* 
					merlin_util_bhazard(gml))
			} 
                }
	}

	//exactly observed events and/or right censoring
	//-> survival function
	gml.survind = 2
	if (merlin_get_nobs(gml,model)) {
		index2 		= merlin_get_surv_index(gml)
		logl[index2,] 	= logl[index2,] :- 
                                        (*gml.Pch[model])(gml,y[index2,1])
	}
        
        //delayed entry
        gml.survind = 4
        if (merlin_get_nobs(gml,model)) {
		index4 		= merlin_get_surv_index(gml)
		logl[index4,] 	= logl[index4,] :+ 
                                        (*gml.Pch[model])(gml,y[index4,3])
	}

	//interval censoring
	//-> cdf function
	gml.survind = 3
	if (merlin_get_nobs(gml,model)) {
		if (gml.hasltrunc[model]) i2 = 4
		else 			  i2 = 3
		index2 = merlin_get_surv_index(gml)	//exit times
		gml.survind = 5
		index5 = merlin_get_surv_index(gml)	//entry times
		if (gml.hastmat) {
                        logl[index3,] = merlin_logl_intcens_ms(gml,
                                                y[,(1,i2)],index3, index5)
                }
		else 			 {
			gml.survind 	= 3
			logl[index3,] 	= (*gml.Pcdf[model])(gml,y[index3,1])
			gml.survind 	= 5
			logl[index5,] 	= logl[index5,] :- 
                                                (*gml.Pcdf[model])(gml,
                                                        y[index5,i2])
			logl[index3,] 	= log(logl[index3,])
		}
	}

	//left truncation handled externally

	gml.survind = 0
	return(logl)
}


`RM' merlin_logl_survival(`Sgml' gml, | `RM' G, `RM' H)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
        haslt	= gml.hasltrunc[model]
	hasbh	= gml.hasbh[model,]
	logl 	= gml.Nrelevels ? 
                        J(gml.Nobs[gml.Nlevels,model],
                                gml.ndim[gml.Nrelevels],0) : 
                        J(merlin_get_nobs(gml,model),1,0)

	//exactly observed events
	//-> hazard function
	gml.survind = 1
	if (merlin_get_nobs(gml,model)) {
		index1 = merlin_get_surv_index(gml)
		if (gml.hastmat) {
                        logl[index1,] 	= merlin_logl_ms(gml,y[index1,1])
                }
		else {
                        logl[index1,] 	= (*gml.Plogh[model])(gml,y[index1,1])

			if (hasbh[1]) {
				logl[index1,] = log(exp(logl[index1,]) :+ 
					merlin_util_bhazard(gml))
			}
			else if (hasbh[2]) {		
				logl[index1,] = log(exp(logl[index1,]) :* 
					merlin_util_bhazard(gml))
			} 
                }
	}

	//exactly observed events and/or right censoring
	//-> survival function
	gml.survind = 2
	if (merlin_get_nobs(gml,model)) {
		index2 		= merlin_get_surv_index(gml)
                if (haslt) {
                        logl[index2,] = logl[index2,] :- 
                                                (*gml.Pch[model])(gml,
                                                        y[index2,1],
                                                        y[index2,3])
                }
                else {
                        logl[index2,] = logl[index2,] :- 
                                                (*gml.Pch[model])(gml,
                                                        y[index2,1])
                }
	}

	//interval censoring
	//-> cdf function
	gml.survind = 3
	if (merlin_get_nobs(gml,model)) {
		if (gml.hasltrunc[model]) 	i2 = 4
		else 				i2 = 3
		index2 = merlin_get_surv_index(gml)	//exit times
		gml.survind = 5
		index5 = merlin_get_surv_index(gml)	//entry times
		if (gml.hastmat) {
                        logl[index3,] = merlin_logl_intcens_ms(gml,
                                                y[,(1,i2)],index3, index5)
                }
		else 			 {
			gml.survind 	= 3
			logl[index3,] 	= (*gml.Pcdf[model])(gml,y[index3,1])
			gml.survind 	= 5
			logl[index5,] 	= logl[index5,] :- 
                                (*gml.Pcdf[model])(gml,y[index5,i2])
			logl[index3,] 	= log(logl[index3,])
		}
	}

	//left truncation (non-multiple record) handled externally
        
	gml.survind = 0
	return(logl)
}

`RM' merlin_logl_ms(`Sgml' gml, `RM' t)
{
	nobs 	= gml.Nsurv[gml.model,3]
	refmod 	= gml.model
	tmat	= gml.tmat

	//find where to go from starting state of gml.model transition
	start	= sum((1::rows(gml.tmat)) :* rowsum(refmod:==gml.tmat))
	Ns 		= sum(tmat[start,]:!=.)

	//logl contribution
	if (Ns==1) 	{
		//if there's only one place to go then it's standard log hazard
		result = (*gml.Plogh[refmod])(gml,t)
	}
	else {
		
		index 		= select(1::rows(tmat),tmat[start,]':!=.)'
		posstates 	= tmat[start,index]
		result 		= (*gml.Plogh[refmod])(gml,t)

		for (i=1;i<=Ns;i++) {
			mod = posstates[i]
			if (gml.issurv[mod]) {
				if (mod!=refmod) {
					gml.model = mod
					result = result :- 
                                                (*gml.Pch[mod])(gml,t)
					gml.model = refmod
				}	
			}
		}
	}
	
	return(result)
}


`RM' merlin_logl_intcens_ms(`Sgml' gml, `RM' t1, `RC' index3, `RC' index5)
{
	nobs 	= gml.Nsurv[gml.model,3]
	refmod 	= gml.model
	tmat	= gml.tmat
	
	//find where to go from starting state of gml.model transition
	start	= sum((1::rows(gml.tmat)) :* rowsum(refmod:==gml.tmat))
	Ns 		= sum(tmat[start,]:!=.)

	//logl contribution
	if (Ns==1) 	{
		//if there's only one place to go then it's standard CDF
		gml.survind = 3
		result = (*gml.Pcdf[model])(gml,t1[index3,1])
		gml.survind = 5
		result[index5,] = result[index5,] :- 
                                        (*gml.Pcdf[model])(gml,t1[index5,2])
		result = log(result)
	}
	else {
		
		index 		= select(1::rows(tmat),tmat[start,]':!=.)'
		posstates 	= tmat[start,index]
		
		if (gml.iccrfudge) {	//comparison for paper

			gml.survind = 3
			result = J(gml.Nobs[refmod],1,0)
			result[index3,] = (*gml.Pcdf[refmod])(gml,t1[index3,1])
			
			gml.survind = 5
			result[index5,] = result[index5,] :- 
                                                (*gml.Pcdf[refmod])(gml,
                                                        t1[index5,2])

			result[index3,] = log(result[index3,])

			for (i=1;i<=Ns;i++) {
				mod = posstates[i]
				if (gml.issurv[mod] & mod!=refmod) {
                                        gml.model = mod
                                        result[index5,] = result[index5,] :- 
                                                        (*gml.Pch[mod])(gml,
                                                                t1[index5,2])
                                        gml.model = refmod
				}
			}
			result = result[index3,]
		}
		else {			// correct logl

			gml.survind = 3	// left 0s handled fine with 
                                        // numerical integration
			result 	= J(nobs,1,0)
			Ngq 	= gml.chip
			gq 	= merlin_gq(Ngq,"legendre")
			qp	= (t1[index3,1]:-t1[index3,2]) :/ 2 :* 
                                        J(nobs,1,gq[,1]') :+ 
                                        (t1[index3,2]:+t1[index3,1]):/2
			logqw	= log((t1[index3,1]:-t1[index3,2]) :/ 2 :* 
                                        J(nobs,1,gq[,2]'))

			for (q=1;q<=Ngq;q++) {
				cr = J(nobs,1,0)
				for (i=1;i<=Ns;i++) {
					mod = posstates[i]
					if (gml.issurv[mod]) {
						if (mod==refmod) {
                                cr = cr :+ (*gml.Plogh[mod])(gml,qp[,q]) :- 
                                        (*gml.Pch[mod])(gml,qp[,q])
						}
						else {
                                gml.model = mod
                                cr = cr :- (*gml.Pch[mod])(gml,qp[,q])
                                gml.model = refmod
						}	
					}
				}
				result = result :+ exp(cr :+ logqw[,q])
			}
			result = log(result)
		}
	}
	
	return(result)
}

end
