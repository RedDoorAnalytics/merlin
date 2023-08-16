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

//Evaluate linear predictors
void merlin_xb(`GMLS' gml, `RR' b)
{
	eqnind = 1
	for (i=1;i<=gml.Nmodels;i++) {
		gml.model = i
		if (gml.NotNull[i,2]) {
			(*gml.Pxb[i])(gml,b,eqnind)
			if (gml.Nap[i]) {
				for (k=1;k<=gml.Nap[i];k++) {
					asarray(gml.apxb,(i,k),b[eqnind++])	
				}
			}
		}
	}
	merlin_fillvcv(gml,b,eqnind)
}

void merlin_xb_1(`GMLS' gml, `RR' b, `RS' eqnind)
{
	eqnind = gml.eqnindex[gml.model,2] + 1
}

void merlin_xb_2(`GMLS' gml, `RR' b, `RS' eqnind)
{
	eqnind = gml.eqnindex[gml.model,2] + 1
	asarray(gml.distancb,(gml.model,1),exp(b[eqnind++]))
}

void merlin_xb_3(`GMLS' gml, `RR' b, `RS' eqnind)
{
	eqnind = gml.eqnindex[gml.model,2] + 1
	asarray(gml.distancb,(gml.model,1),b[eqnind++])
}

void merlin_xb_rcs(`GMLS' gml, `RR' b, `RS' eqnind)
{
	mod 	= gml.model
	eqnind 	= gml.eqnindex[mod,2] + 1
	asarray(gml.distancb,(mod,1),b[|eqnind\(eqnind+gml.Ndistancp[mod]-1)|]')		//spline coefficients
	eqnind = eqnind + gml.Ndistancp[mod]
}

void merlin_xb_ord(`GMLS' gml, `RR' b, `RS' eqnind)
{
	mod 	= gml.model
	eqnind 	= gml.eqnindex[mod,2] + 1
	Ncuts	= gml.Ndistancp[mod]
	cuts 	= J(Ncuts,1,.)
	for (k=1; k<=Ncuts; k++) cuts[k] = b[eqnind++]
	asarray(gml.distancb,(mod,1),cuts)
}

// fill variance-covariance parameters
void merlin_fillvcv(	`GMLS' gml, 
					`RR' b,
					`RS' eqnind)
{
	`RM' covariances, vcv, sdcor
	`RC' Nres
	`RS' adapt

	covariances = gml.covariances	//indep \ exch \ unstr
	Nres 		= gml.Nres
	adapt 		= gml.adapt

	for (i=1 ; i < gml.Nlevels ; i++) {

		if (Nres[i]) {

			/*if (gml.hasRExb[i]) {

				vcv = sdcor = J(gml.Nobs[i,1],rows(vech(I(Nres[i]))),0) 
				if (covariances[1,i]) {
					index = 1
					for (j=1;j<=Nres[i];j++) {
						xb = exp(moptimize_util_xb(M,b,eqnind++))
						if (rows(xb)>1) sdcor[,index] = select(xb,gml.oneidvars[,i])
						else 			sdcor[,index] = J(gml.Nobs[i,1],1,xb)
						vcv[,index] = sdcor[,index] :^2
						index = index + Nres[i] - (j-1)
					}
				}
			}
			else {*/

				sdcor = vcv = asarray(gml.vcvs,i)
				if (covariances[1,i]) {
					for (j=1;j<=Nres[i];j++) {
						sdcor[j,j] = exp(b[eqnind++])
						vcv[j,j] = sdcor[j,j] :^2
					}
				}
				else if (covariances[2,i]) {
					var_xb = exp(b[eqnind++]):^2
					cov_xb = var_xb :* tanh(b[eqnind++])
					for (j=1;j<=Nres[i];j++) {
						ind = 1
						while (ind<j) {
							vcv[ind,j] = vcv[j,ind] = cov_xb
							ind++  
						}	
						vcv[j,j] = var_xb
					}
				}
				else if (covariances[3,i]) {
					for (j=1;j<=Nres[i];j++) {
						sdcor[j,j] = exp(b[eqnind++])
						vcv[j,j] = sdcor[j,j] :^2
					}
					for (j=1;j<=Nres[i];j++) {
						ind = 1
						while (ind<j) {
							sdcor[ind,j] = sdcor[j,ind] = tanh(b[eqnind++])
							ind++
						}
						ind = 1
						while (ind<j) {
							vcv[ind,j] = vcv[j,ind] = sdcor[ind,ind]:*sdcor[j,j]:*sdcor[ind,j]
							ind++
						}
					}
				
				}
				else {
					var_xb = exp(b[eqnind++]):^2
					for (j=1;j<=Nres[i];j++) vcv[j,j] = var_xb
				}			
			//}
			asarray(gml.vcvs,i,vcv)

			merlin_update_ip(gml,i)
		}
	}
	
}

void merlin_update_ip(`GMLS' gml, `RS' i)
{
	//update nodes
	if (gml.adapt[i]) {
		res = cholesky(asarray(gml.vcvs,i)) * asarray(gml.stackednodes,i)	
		for (r=1;r<=gml.Nres[i];r++) {
			asarray(gml.aghip2,(i,r),rowshape(res[r,],gml.Nobs[i,1]))
		}	
	}
	else {
		if (gml.usegh[i]) {
// 			if (gml.hasRExb[i]) {
// 				vcv = asarray(gml.vcvs,i)
// 				ndim = gml.ndim[i]
// 				ind1 = 1; ind2 = ndim
// 				stackednodes = J(gml.ndim[i]*gml.Nobs[i,1],gml.Nres[i],.)
// 				for (j=1;j<=gml.Nobs[i,1];j++) {
// 					nodes = cholesky(invvech(vcv[j,])) * asarray(gml.baseGHnodes,i)
// 					stackednodes[|ind1,.\ind2,.|] = nodes'
// 					ind1 = ind1 + ndim
// 					ind2 = ind2 + ndim
// 				}
// 				stackednodes = stackednodes'			//!!fix
// 				for (r=1;r<=gml.Nres[i];r++) {
// 					asarray(gml.aghip2,(i,r),rowshape(stackednodes[r,],gml.Nobs[i,1]))
// 				}				
// 			}
// 			else {
				bnodes = cholesky(asarray(gml.vcvs,i)) * asarray(gml.baseGHnodes,i)
				asarray(gml.b,i,bnodes)	
// 			}
		}
		else {
			if (gml.renormal[i]) 	{
				asarray(gml.b,i,cholesky(asarray(gml.vcvs,i)) * asarray(gml.bdraws,i))
			}
			else {	
				rseed(gml.seed)			//reset seed
				bnodes = merlin_drawt(J(gml.Nres[i],1,0),asarray(gml.vcvs,i),gml.df[i],gml.ndim[i])
				asarray(gml.b,i,bnodes)	
			}
		}
			
	}	
}

end
