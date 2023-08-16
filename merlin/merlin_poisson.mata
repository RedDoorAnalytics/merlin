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

`RM' merlin_logl_poisson_ml(	`gml' gml)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs    = merlin_get_nobs(gml,model)
	logl 	= J(Nobs,gml.ndim[gml.Nrelevels],0)
        
	//core
	xb	= merlin_util_xzb(gml)

	//===========================================================================================================================//
	// log likelihood
	
        logl = y :* xb :- exp(xb) :- lngamma(y:+1)

	if (gml.todo==0) return(logl)
}

`RM' merlin_logl_poisson(`gml' gml, 	///
			`RM' G, 	///
			`RM' H)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml,model)
	logl 	= J(Nobs,1,0)

	//core
	xb	= merlin_util_xzb(gml)

	//===========================================================================================================================//
	// log likelihood
	
		logl = y :* xb :- exp(xb) :- lngamma(y:+1)

	if (gml.todo==0) return(logl)

	//===========================================================================================================================//
	// score
		
		//indexes for covariates and baseline spline equations
		NHbs 	= asarray(gml.NHbs,model)
		sindex	= (1..NHbs[1]) :+ gml.skip[model]
		
		//core
                index   = merlin_get_index(gml)
		x	= merlin_util_xz(gml)
		
		G[index,sindex] = x :* (y :- exp(xb))
		
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
	
		Hx = - x[,xindex1] :* exp(xb) :* x[,xindex2]
		
		//build H
		Hxsum 	= quadcolsum(Hx,1)
		
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
		
		Hindex1 = 1 + gml.skip[model] 
		Hindex2 = NHbs[1] + gml.skip[model] 
		H[|Hindex1,Hindex1\Hindex2,Hindex2|] = Hm
		
		return(logl)
}


`RM' merlin_logl_pois(`gml' gml)
{
	y 		= merlin_util_depvar(gml)
	xzb 	= merlin_util_xzb(gml)
	return(y :* xzb :- exp(xzb) :- lngamma(y:+1))
}

`RM' merlin_poisson_expval(`gml' gml, | `RC' t)
{
	if (args()==1) 	return(exp(merlin_util_xzb(gml)))
	else 			return(exp(merlin_util_xzb(gml,t)))
}

end
