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

`RM' merlin_logl_cox_firth(`gml' gml, `RM' S, `RM' H)
{	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml)
	haslt	= gml.hasltrunc[model]

	logl 	= J(Nobs,1,0)
	
	nb = asarray(gml.NHbs,gml.model)
	ind1 = ind2 = J(1,0,.)
	for (i=1; i<=nb; i++) {
		refind = 1
		while (refind<=i) {
			ind1 = ind1,i
			ind2 = ind2,refind
			refind++
		}
	}
	H = J(Nobs,cols(ind1),0)
	
	//exactly observed events
	gml.survind = 1
	Nfails = merlin_get_nobs(gml)
	index1 = merlin_get_surv_index(gml)

	//xb at failure times
	logl[index1,] = merlin_util_xzb_simple(gml)

	//x at failure times
	S 		= J(Nobs,asarray(gml.NHbs,model),0)	
	S[index1,] 	= merlin_util_xz_simple(gml)

	//xb at all survival/censoring times
	gml.survind = 0

	if (!gml.istimedep[model,1]) {									//time-independent
		expxb2 		= exp(merlin_util_xzb_simple(gml))
		x2			= merlin_util_xz_simple(gml)	
		if (haslt) {
			for (i=1;i<=Nfails;i++) {
				atriskid 	 = selectindex((y[index1[i],1] :<= y[,1]) :& (y[index1[i],1] :> y[,3]))
				expxb2sum 	 = quadcolsum(expxb2[atriskid,])
				expxb2sumx2	 = quadcolsum(expxb2[atriskid,] :* x2[atriskid,])
				logl[index1[i],] = logl[index1[i],] :- log(expxb2sum)
				S[index1[i],] 	= S[index1[i],] :- expxb2sumx2 :/ expxb2sum
				H[index1[i],] = expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ 	///
					(expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* 	///
					x2[atriskid,ind1] :* x2[atriskid,ind2]) :/ expxb2sum
			}
		}
		else {
			for (i=1;i<=Nfails;i++) {
				atriskid  = selectindex(y[index1[i],1] :<= y[,1])
				expxb2sum = quadcolsum(expxb2[atriskid,])
				logl[index1[i],] = logl[index1[i],] :- log(expxb2sum)
				expxb2sumx2	= quadcolsum(expxb2[atriskid,] :* x2[atriskid,])
				S[index1[i],] 	= S[index1[i],] :- expxb2sumx2 :/ expxb2sum
				H[index1[i],]	= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ 	///
					(expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* 	///
					x2[atriskid,ind1] :* x2[atriskid,ind2],1) :/ expxb2sum
			}
		}
	}
	else {															//time-dependent
		if (haslt) {
			for (i=1;i<=Nfails;i++) {
				atriskid 	= selectindex((y[index1[i],1] :<= y[,1]) :& (y[index1[i],1] :> y[,3]))
				expxb2 		= exp(merlin_util_xzb_simple(gml,J(Nobs,1,y[index1[i],1])))
				expxb2sum 	= quadcolsum(expxb2[atriskid,],1)
				logl[index1[i],] = logl[index1[i],] :- log(expxb2sum)
				x2		= merlin_util_xz_simple(gml,J(Nobs,1,y[index1[i],1]))
				S[index1[i],] 	= S[index1[i],] :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,]) :/ expxb2sum
				expxb2sumx2	= quadcolsum(expxb2[atriskid,] :* x2[atriskid,],1)
				H[index1[i],]	= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ 	///
					(expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* 	///
					x2[atriskid,ind1] :* x2[atriskid,ind2],1) :/ expxb2sum
			}
		}
		else {
			for (i=1;i<=Nfails;i++) {
				atriskid 	= selectindex(y[index1[i],1] :<= y[,1])
				expxb2 		= exp(merlin_util_xzb_simple(gml,J(Nobs,1,y[index1[i],1])))
				expxb2sum 	= quadcolsum(expxb2[atriskid,],1)
				logl[index1[i],] = logl[index1[i],] :- log(expxb2sum)
				x2		= merlin_util_xz_simple(gml,J(Nobs,1,y[index1[i],1]))
				S[index1[i],] 	= S[index1[i],] :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,],1) :/ expxb2sum
				expxb2sumx2	= quadcolsum(expxb2[atriskid,] :* x2[atriskid,],1)
				H[index1[i],]	= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ ///
					(expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* 	///
					x2[atriskid,ind1] :* x2[atriskid,ind2],1) :/ expxb2sum
			}
		}
	}
	
	Hxsum 	= quadcolsum(H,1)
	H	= J(nb,nb,.)	
	el 	= 1
	for (e1=1;e1<=nb;e1++) {
		e2 = 1
		while (e2<=e1) {
			if (e1==e2) H[e1,e1] = Hxsum[el++]
			else 	H[e2,e1] = H[e1,e2] = Hxsum[el++]
			e2++
		}
	}
	st_matrix("firthV",(1:/(-H)))
	return(logl :+ 0.5 :* log(det(-H)) :/ gml.N)
}

end
