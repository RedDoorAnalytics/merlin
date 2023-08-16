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

`RM' merlin_cox_score_clp(`gml' gml)
{	
	model 	= gml.model
	y 		= merlin_util_depvar(gml)
	S 		= J(merlin_get_nobs(gml,model),asarray(gml.NHbs,model),0)

	//exactly observed events
	gml.survind = 1
	
	//get index for sorted failure times
	if (gml.imputing) 	{
		Nfails = gml.NsurvImp[model,1]
		index1 = asarray(gml.surv_index_imp,(model,1))
	}
	else	{
		Nfails = gml.Nsurv[model,1]
		index1 = asarray(gml.surv_index,(model,1))
	}
	
	//x at failure times
	S[index1,] = merlin_util_xz(gml)

	//xb at all survival/censoring times
	gml.survind = 0

	if (!(gml.hastvars[model,1] | gml.hastvars[model,2])) {

		expxb2 		= exp(merlin_util_xzb(gml,y[,1]))
		x2			= merlin_util_xz(gml)
		for (i=1;i<=Nfails;i++) {
			atriskid 		= selectindex(y[index1[i],1] :<= y[,1])
			S[index1[i],] 	= S[index1[i],] :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,]) :/ quadcolsum(expxb2[atriskid,])
		}
		
	}
	else {
		
		//time-dependent
		for (i=1;i<=Nfails;i++) {
			//so inefficient...........
			atriskid 			= selectindex(y[index1[i],1] :<= y[,1])
			expxb2 				= exp(merlin_util_xzb(gml,J(merlin_get_nobs(gml,model),1,y[index1[i],1])))
			x2					= merlin_util_xz(gml,J(merlin_get_nobs(gml,model),1,y[index1[i],1]))
			S[index1[i],] 		= S[index1[i],] :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,]) :/ (quadcolsum(expxb2[atriskid,]))
		}
	
	}

	return(S)
}

end
