*! version 1.0.0 ?????2016

local gml 	struct merlin_struct scalar
local pgml	pointer(struct merlin_struct scalar) scalar
local TR 	transmorphic
local RS 	real scalar
local RC 	real colvector
local SS 	string scalar
local PS 	pointer scalar
local RR 	real rowvector
local RM 	real matrix
local PC 	pointer colvector
local PM 	pointer matrix
local SC 	string colvector

version 14.1

mata:

`RM' merlin_logl_cox(`gml' gml, `RM' S, `RM' H)
{	
 	if (gml.firth) return(merlin_logl_cox_firth(gml,S,H))
	
	model 	= gml.model
	y 	= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml)
	haslt	= gml.hasltrunc[model]

	logl 	= J(Nobs,1,0)
	
	if (gml.todo>1) {
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
	}
	
	//exactly observed events
	gml.survind = 1
	Nfails = merlin_get_nobs(gml)
	index1 = merlin_get_surv_index(gml)

	//xb at failure times
	logl[index1,] = merlin_util_xzb_simple(gml)

	//x at failure times
	if (gml.todo) {
		S 		= J(Nobs,asarray(gml.NHbs,model),0)	
		S[index1,] 	= merlin_util_xz_simple(gml)
	}

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
				if (gml.todo) {
					S[index1[i],] 	= S[index1[i],] :- expxb2sumx2 :/ expxb2sum
					if (gml.todo>1) H[index1[i],] = expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ (expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,ind1] :* x2[atriskid,ind2]) :/ expxb2sum
				}
			}
		}
		else {
			for (i=1;i<=Nfails;i++) {
				atriskid  = selectindex(y[index1[i],1] :<= y[,1])
				expxb2sum = quadcolsum(expxb2[atriskid,])
				logl[index1[i],] = logl[index1[i],] :- log(expxb2sum)
				if (gml.todo) {
					expxb2sumx2	= quadcolsum(expxb2[atriskid,] :* x2[atriskid,])
					S[index1[i],] 	= S[index1[i],] :- expxb2sumx2 :/ expxb2sum
					if (gml.todo>1) H[index1[i],]	= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ (expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,ind1] :* x2[atriskid,ind2],1) :/ expxb2sum
				}
			}
		}
	}
	else {															//time-dependent
		if (haslt) {
			for (i=1;i<=Nfails;i++) {
				atriskid 			= selectindex((y[index1[i],1] :<= y[,1]) :& (y[index1[i],1] :> y[,3]))
				expxb2 				= exp(merlin_util_xzb_simple(gml,J(Nobs,1,y[index1[i],1])))
				expxb2sum 			= quadcolsum(expxb2[atriskid,],1)
				logl[index1[i],] 	= logl[index1[i],] :- log(expxb2sum)
				if (gml.todo) {
					x2					= merlin_util_xz_simple(gml,J(Nobs,1,y[index1[i],1]))
					S[index1[i],] 		= S[index1[i],] :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,]) :/ expxb2sum
					if (gml.todo>1) {
						expxb2sumx2			= quadcolsum(expxb2[atriskid,] :* x2[atriskid,],1)
						H[index1[i],]		= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ (expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,ind1] :* x2[atriskid,ind2],1) :/ expxb2sum
					}
				}
			}
		}
		else {
			for (i=1;i<=Nfails;i++) {
				atriskid 			= selectindex(y[index1[i],1] :<= y[,1])
				expxb2 				= exp(merlin_util_xzb_simple(gml,J(Nobs,1,y[index1[i],1])))
				expxb2sum 			= quadcolsum(expxb2[atriskid,],1)
				logl[index1[i],] 	= logl[index1[i],] :- log(expxb2sum)
				if (gml.todo) {
					x2					= merlin_util_xz_simple(gml,J(Nobs,1,y[index1[i],1]))
					S[index1[i],] 		= S[index1[i],] :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,],1) :/ expxb2sum
					if (gml.todo>1) {
						expxb2sumx2			= quadcolsum(expxb2[atriskid,] :* x2[atriskid,],1)
						H[index1[i],]		= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ (expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,ind1] :* x2[atriskid,ind2],1) :/ expxb2sum
					}
				}
			}
		}
	}
	
	if (gml.todo==2) {
		Hxsum 	= quadcolsum(H,1)
		H		= J(nb,nb,.)	
		el 		= 1
		for (e1=1;e1<=nb;e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) H[e1,e1] = Hxsum[el++]
				else 		H[e2,e1] = H[e1,e2] = Hxsum[el++]
				e2++
			}
		}
	}
	
	return(logl)
}

`RM' merlin_logl_cox2(`gml' gml, `RM' S, `RM' H)
{	
	model 	= gml.model
	y 		= merlin_util_depvar(gml)
	Nobs	= merlin_get_nobs(gml)
	haslt	= gml.hasltrunc[model]

	logl 	= J(Nobs,1,0)
	
	if (gml.todo) {
		nb = asarray(gml.NHbs,gml.model)
		if (gml.todo==2) {
			ind1 = ind2 = J(1,0,.)
			for (i=1; i<=nb; i++) {
				refind = 1
				while (refind<=i) {
					ind1 = ind1,i
					ind2 = ind2,refind
					refind++
				}
			}
			NHb = cols(ind1)
			H 	= J(Nobs,NHb,0)
		}
	}
	
	//exactly observed events
	gml.survind = 1
	Nfails = merlin_get_nobs(gml)
	index1 = merlin_get_surv_index(gml)

	//xb at failure times
	logl[index1,] = merlin_util_xzb_simple(gml)

	//x at failure times
	if (gml.todo) {
		S 			= J(Nobs,asarray(gml.NHbs,model),0)	
		S[index1,] 	= merlin_util_xz_simple(gml)
	}

	//xb at all survival/censoring times
	gml.survind = 0
	if (!gml.istimedep[model,1]) {	//time-independent
		
		expxb2 = exp(merlin_util_xzb_simple(gml))
		if (gml.todo) x2 = merlin_util_xz_simple(gml)
		
		if (haslt) {
			for (i=1;i<=Nfails;i++) {
				atriskid 		= asarray(gml.cox_index,(model,2,i+6))
				expxb2sum 		= quadcolsum(expxb2[atriskid,])
				logl[index1[i],] 	= logl[index1[i],] :- log(expxb2sum)
				if (gml.todo) {
					expxb2sumx2	= quadcolsum(expxb2[atriskid,] :* x2[atriskid,])
					S[index1[i],] 	= S[index1[i],] :- expxb2sumx2 :/ expxb2sum
					if (gml.todo>1) H[index1[i],]	= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ (expxb2sum):^2 :- quadcolsum(expxb2[atriskid,] :* x2[atriskid,ind1] :* x2[atriskid,ind2]) :/ expxb2sum
				}
			}
		}
		else {
			//done
			gml.survind 	= 2
			ord 		= merlin_get_surv_index(gml) 	//index to put times in descending order
			gml.survind	= 3
			invord 		= merlin_get_surv_index(gml)	//index to put back to original order
			sumexpxb2 	= quadrunningsum(expxb2[ord],1)[invord]
			logl[index1] 	= logl[index1] :- log(sumexpxb2[index1])
			if (gml.todo) {
				expxb2x2	= expxb2[ord] :* x2[ord,]
				sumexpxb2x2 = J(Nfails,nb,.)
				for (b=1;b<=nb;b++) {
					sumexpxb2x2[,b] = quadrunningsum(expxb2x2[,b],1)[invord,][index1,]
				}
				S[index1,] 	= S[index1,] :- sumexpxb2x2 :/ sumexpxb2[index1]
				if (gml.todo==2) {
					Hexpxb2x2		= expxb2[ord] :* x2[ord,ind1] :* x2[ord,ind2]
					Hsumexpxb2x2 	= J(Nfails,NHb,.)
					for (b=1;b<=NHb;b++) {
						Hsumexpxb2x2[,b] = quadrunningsum(Hexpxb2x2[,b],1)[invord,][index1,]
					}
					H[index1,] = sumexpxb2x2[,ind2] :* sumexpxb2x2[,ind1] :/ (sumexpxb2[index1]):^2 :- Hsumexpxb2x2 :/ sumexpxb2[index1]
				}
			}
		}
	}
        else {	//time-dependent
		for (i=1;i<=Nfails;i++) {
			Natrisk 			= asarray(gml.cox_index,(model,1,i+6))								//ltrunc handled in setup
			gml.survind			= i + 6
			expxb2 				= exp(merlin_util_xzb_simple(gml,J(Natrisk,1,y[index1[i],1])))
			expxb2sum 			= quadcolsum(expxb2,1)
			logl[index1[i],] 	= logl[index1[i],] :- log(expxb2sum)
			if (gml.todo) {
				x2					= merlin_util_xz_simple(gml,J(Natrisk,1,y[index1[i],1]))
				S[index1[i],] 		= S[index1[i],] :- quadcolsum(expxb2 :* x2) :/ expxb2sum
				if (gml.todo>1) {
					expxb2sumx2			= quadcolsum(expxb2 :* x2,1)
					H[index1[i],]		= expxb2sumx2[,ind2] :* expxb2sumx2[,ind1] :/ (expxb2sum):^2 :- quadcolsum(expxb2 :* x2[,ind1] :* x2[,ind2],1) :/ expxb2sum
				}
			}
		}	
	}

	if (gml.todo==2) {
		Hxsum 	= quadcolsum(H,1)
		H		= J(nb,nb,.)	
		el 		= 1
		for (e1=1;e1<=nb;e1++) {
			e2 = 1
			while (e2<=e1) {
				if (e1==e2) H[e1,e1] = Hxsum[el++]
				else 		H[e2,e1] = H[e1,e2] = Hxsum[el++]
				e2++
			}
		}
	}
	
	return(logl)
}

end
