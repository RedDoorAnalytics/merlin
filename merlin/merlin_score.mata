*! version 1.0.0 MJC

local TS 	transmorphic scalar
local RS 	real scalar
local RC 	real colvector
local RM 	real matrix
local Pgml 	pointer(struct merlin_struct scalar) scalar
local gml 	struct merlin_struct scalar

version 14.2
mata:

/*
scores for multilevel model (without any ?EV[], ?XB[]

- one function for complex predictor
- one function for each distributional ancillary parameter
- random effects done with deriv()
*/

void merlin_score_panels(`gml' gml, `RC' lnfi, `RM' G)
{
	analytic = 1
	G 	 = J(gml.Nobs[1,1],gml.Nb,0)
	bindex 	 = 1
	
	if (analytic) {
		
		Li = exp(lnfi)
		for (mod=1;mod<=gml.Nmodels;mod++) {
			
			gml.model = gml.modtoind = mod
			NHbs 	  = asarray(gml.NHbs,mod)	
			Nbs	  = sum(NHbs)
			
			//clp
			hasconstr = asarray(gml.hasconstraint,mod)
			for (el=1;el<=NHbs[1];el++) {
				if (!hasconstr[el]) {
					G[,bindex] = merlin_panels(1,gml,	
						&merlin_weibull_score_loglambda(),bindex) 
				}
				bindex++
			}
			
			//dap
			if (gml.Ndistancp[mod]) {
				for (eqn=1;eqn<=gml.Ndistancp[mod];eqn++) {
					G[,bindex] = merlin_panels(1,gml,&merlin_weibull_score_loggamma())
					bindex++
				}
			}
			//ap
			
		}

		G[,1..Nbs] = G[,1..Nbs] :/ Li
	}
	else {
	
		for (mod=1;mod<=gml.Nmodels;mod++) {
			
			gml.model 	= gml.modtoind = mod
			NHbs 		= asarray(gml.NHbs,mod)	
			
			for (eqn=1;eqn<=gml.NHeqns[mod];eqn++) {
				if (eqn==1) {
					hasconstr = asarray(gml.hasconstraint,mod)
					for (el=1;el<=NHbs[eqn];el++) {
						if (!hasconstr[el]) G[,bindex] = merlin_dscore_panels(gml,bindex)
						bindex++
					}
				}
				else {
					for (el=1;el<=NHbs[eqn];el++) {
						G[,bindex] = merlin_dscore_panels(gml,bindex)
						bindex++
					}
				}
			}
		}
		
	}
	
	//vcv - derivatives found numerically
	eqnind 	= bindex
	Nres 	= gml.Nres
	covariances = gml.covariances
	
	for (i=1 ; i < gml.Nlevels ; i++) {

		if (Nres[i]) {

			if (covariances[1,i]) {
				for (j=1;j<=Nres[i];j++) {
					G[,eqnind] = merlin_deriv(gml,eqnind)
					eqnind++
				}
			}
			else if (covariances[2,i]) {
				G[,eqnind] = merlin_deriv(gml,eqnind)
					eqnind++
			}
			else if (covariances[3,i]) {
				for (j=1;j<=Nres[i];j++) {
					G[,eqnind] = merlin_deriv(gml,eqnind)
					eqnind++
				}
				for (j=1;j<=Nres[i];j++) {
					ind = 1
					while (ind<j) {
						G[,eqnind] = merlin_deriv(gml,eqnind)
						eqnind++
						ind++
					}
				}
			}
			else {
				G[,eqnind] = merlin_deriv(gml,eqnind)
				eqnind++
			}			
		}
	}
	
}

`RC' merlin_panels(`RS' index,		/// -level-
		   `gml' gml,		/// -merlin object-
		   `PS' func,		/// -function to call-
		   | `RS' Xindex)	// -design matrix element-
{
	`RS' index2
	`RM' res, panelindex
	
	if (args()==3) Xindex = 1
	
	index2 = index+1
	
	resq = J(gml.Nobs[index,1],gml.ndim[index],0)
	
	if (index<gml.Nrelevels) {
		panelindex = asarray(gml.panelindexes,(index,1))
		for (q=1;q<=gml.ndim[index];q++) {
			gml.qind[1,index2] = q
			resq[,q] = panelsum(merlin_panels(index2,gml,func,Xindex),panelindex)
		}
	}
	else {
		for (j=1;j<=gml.Nmodels;j++) {
			gml.model = gml.modtoind = j
			if (gml.NotNull[j,1]) {
				resq2 = (*func)(gml,Xindex)
				if (gml.hasweights[index2])	resq2 = resq2 :* asarray(gml.weights,(index2,j))
				resq = resq :+ panelsum(resq2,asarray(gml.panelindexes,(index,j)))
			}
		}
	}
 
	resq = resq :* asarray(gml.Li_ip,gml.qind) 

	if (gml.usegh[index]) {			//GHQ
		if (gml.hasweights[index]) 	return(asarray(gml.weights,(index,1)) :* resq * asarray(gml.baseGHweights,index))
		else 				return(resq * asarray(gml.baseGHweights,index))
	}
	else {					//MCI
		if (gml.hasweights[index]) 	return(asarray(gml.weights,(index,1)) :* quadrowsum(resq):/gml.ndim[index])
		else 				return(quadrowsum(resq):/gml.ndim[index])
	}
}
void merlin_logl_panels_d(`RS' myb,`gml' gml, `RS' b, `RC' lnf)
{
	gml.myb[b] = myb
	merlin_xb(gml,gml.myb)
	lnf = merlin_logl_panels(1,gml)
}

`RM' merlin_deriv(`gml' gml, `RS' eqnind)
{
	
	myb = gml.myb
	D   = deriv_init()
	deriv_init_evaluator(D, &merlin_logl_panels_d())
	deriv_init_evaluatortype(D, "v")
// 	deriv_init_search(D,"off")
	deriv_init_argument(D, 1, gml)
	gml.myb = myb		
	deriv_init_params(D, myb[eqnind])
	deriv_init_argument(D, 2, eqnind)
	ddb = deriv(D, 1)
	gml.myb = myb
	return(deriv_result_scores(D))
}


`RC' merlin_dscore_panels(`gml' gml, `RS' eqnind)
{
	b = gml.myb[eqnind]				
	hstep = c("epsdouble")^(1/3)
	if (abs(b)<=1 & abs(b)>0) {
		hstep = abs(b)*hstep
	}
	copymyb = gml.myb
	gml.myb[eqnind] = copymyb[eqnind] :+ hstep
	merlin_xb(gml,gml.myb)
	lnf1 = merlin_logl_panels(1,gml)
	gml.myb[eqnind] = copymyb[eqnind] :- hstep
	merlin_xb(gml,gml.myb)
	lnf2 = merlin_logl_panels(1,gml)
	gml.myb = copymyb
	return((lnf1:-lnf2):/(2:*hstep))
}

end
