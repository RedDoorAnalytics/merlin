*! version 1.0.0 ?????2016

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
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

// expval = (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
// return(lnnormalden(merlin_util_depvar(gml), expval, asarray(gml.distancb,(gml.model,1))))

`RM' merlin_logl_gaussian(	`gml' gml, 	///
                                `RM' G, 	///
                                `RM' H)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml,model)
	logl 	= J(Nobs,1,0)

	//core
	xb	= merlin_util_xzb(gml)
	sd	= merlin_util_dap(gml,1)

	//==================================================================//
	// log likelihood
	
                logl = lnnormalden(y,xb,sd)

                if (gml.haspenalty) logl = logl :- merlin_get_penalty(gml):/gml.N
                        
        if (gml.todo==0) return(logl)

	//===========================================================================================================================//
	// score
		
		//indexes for covariates and baseline spline equations
		NHbs 	= asarray(gml.NHbs,model)
		sindex1 = (1..NHbs[1]) :+ gml.skip[model]
		sindex2	= ((NHbs[1]+1)..(NHbs[1]+NHbs[2]))  :+ gml.skip[model]
		
		//core
                index   = merlin_get_index(gml)
		x	= merlin_util_xz(gml)
		
		G[index,sindex1] = x :* (y :- xb)  :/ (sd:^2)
		G[index,sindex2] = -1 :+ (y :- xb):^ 2 :/ (sd :^ 2)
		
		if (gml.haspenalty) {
			cmpbix	= asarray(gml.CmpXBIndex,(gml.model,2))
			bindex 	= 1..(max(cmpbix)-1)
			G[index,bindex] = G[index,bindex] :- merlin_get_deriv_penalty(gml):/gml.N
		}
	
	if (gml.todo==1) return(logl)

	//===========================================================================================================================//
	// hessian

		//x indices
		xindex1 = xindex2 = J(1,0,.)
		for (i=1; i<=NHbs[1]; i++) {
			refind = 1
			while (refind<=i) {
				xindex1 = xindex1,i
				xindex2 = xindex2,refind
				refind++
			}
		}
		Hx = J(Nobs,cols(xindex1),0)
		
		//x,gamma indices
		xindex3 = xindex4 = J(1,0,.)
		for (i=1; i<=NHbs[2]; i++) {
			for (j=1;j<=NHbs[1];j++) {
				xindex3 = xindex3,i
				xindex4 = xindex4,j
			}
		}
		Hxs = J(Nobs,cols(xindex3),0)
		
		//gamma indices
		xindex5 = xindex6 = J(1,0,.)
		for (i=1; i<=NHbs[2]; i++) {
			refind = 1
			while (refind<=i) {
				xindex5 = xindex5,i
				xindex6 = xindex6,refind
				refind++
			}
		}
		Hs = J(Nobs,1,0)

		Hx = - x[,xindex1] :* x[,xindex2] :/ (sd:^2)
		Hxs = -2 :* x[,xindex4] :* (y :- xb)  :/ (sd:^2)
		Hs = -2 :* (y :- xb):^ 2 :/ (sd :^ 2)
		
		//build H
		Hxsum 	= quadcolsum(Hx,1)
		Hxssum	= quadcolsum(Hxs,1)
		Hssum 	= quadcolsum(Hs,1)
		
		Hm	= J(NHbs[1],NHbs[1],.)	
		el 	= 1
		for (e1=1;e1<=NHbs[1];e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) Hm[e1,e1] = Hxsum[el++]
				else 		Hm[e2,e1] = Hm[e1,e2] = Hxsum[el++]
				e2++
			}
		}
		
		Hoff = J(NHbs[2],NHbs[1],.)
		el = 1				
		for (e1=1;e1<=NHbs[2];e1++) {
			for (e2=1;e2<=NHbs[1];e2++) {
				Hoff[e1,e2] = Hxssum[el++]
			}					
		} 
		
		H2	= J(NHbs[2],NHbs[2],.)	
		el 	= 1
		for (e1=1;e1<=NHbs[2];e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) H2[e1,e1] = Hssum[el++]
				else 		H2[e2,e1] = H2[e1,e2] = Hssum[el++]
				e2++
			}
		}
		
		Hm = Hm,Hoff'\Hoff,H2
		
		Hindex1 = 1 + gml.skip[model] 
		Hindex2 = NHbs[1] + NHbs[2] + gml.skip[model] 
		H[|Hindex1,Hindex1\Hindex2,Hindex2|] = Hm
		
		return(logl)
}

`RM' merlin_logl_gaussian_ml(`gml' gml)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs    = merlin_get_nobs(gml,model)
	logl 	= J(Nobs,gml.ndim[gml.Nrelevels],0)

	//core
	xb	= merlin_util_xzb(gml)
	sd	= merlin_util_dap(gml,1)

	//===========================================================================================================================//
	// log likelihood
	
		logl = lnnormalden(y,xb,sd)

		if (gml.haspenalty) logl = logl :- merlin_get_penalty(gml):/gml.N
		
	if (gml.todo==0) return(logl)
}

`RM' merlin_logl_gauss(`gml' gml)
{
	expval = (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
	return(lnnormalden(merlin_util_depvar(gml), expval, asarray(gml.distancb,(gml.model,1))))
}

`RM' merlin_gaussian_expval(`gml' gml, | `RC' t)
{
	if (args()==1) 	return((*gml.invlinks[gml.model])(merlin_util_xzb(gml)))
	else 		return((*gml.invlinks[gml.model])(merlin_util_xzb(gml,t)))
}

`RM' merlin_gaussian_resid(`gml' gml, `RC' t)
{
	expval = (*gml.invlinks[gml.model])(merlin_util_xzb(gml,t))
	return(merlin_util_depvar(gml) :- expval)
}

`RM' merlin_gaussian_rstand(`gml' gml, `RC' t)
{
	expval = (*gml.invlinks[gml.model])(merlin_util_xzb(gml,t))
	return((merlin_util_depvar(gml) :- expval):/asarray(gml.distancb,(gml.model,1)))
}

`RM' merlin_logl_igaussian(`gml' gml)
{
	expval = (*gml.invlinks[gml.model])(merlin_util_xzb(gml))
	return(lnnormalden(gml.ImputedValue, expval, asarray(gml.distancb,(gml.model,1))))
}

end
