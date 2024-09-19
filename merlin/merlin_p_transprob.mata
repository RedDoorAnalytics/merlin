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

`RM' merlin_p_transprob(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	//hard coded to illness-death or extended illness-death
	tmattype = merlin_p_check_tmat(gml)
	tp = strtoreal(st_local("transprob"))

	if (tmattype==1) {		//standard illness-death
		if 	(tp==1) 	pred = merlin_p_tp_p11(gml,t)
		else if (tp==2) 	pred = merlin_p_tp_p12(gml,t)
		else 			pred = merlin_p_tp_p13(gml,t)
	}
	else {					//extended illness-death
		if 	(tp==1) 	pred = merlin_p_tp_p11(gml,t)
		else if (tp==2) 	pred = merlin_p_tp_p12(gml,t)
		else if (tp==3)		pred = merlin_p_tp_p113(gml,t)
		else 			pred = merlin_p_tp_p123(gml,t)
	}
	
	return(pred)	
}

`RM' merlin_p_los(`gml' gml, | `RC' t)
{
	not = args()==1
	if (not) t = merlin_util_timevar(gml)

	//hard coded to illness-death or extended illness-death
	tmattype = merlin_p_check_tmat(gml)
	los = strtoreal(st_local("los"))

	if (tmattype==1) {		//standard illness-death
		if 	(los==1) 	pred = &merlin_p_tp_p11()
		else if (los==2) 	pred = &merlin_p_tp_p12()
		else 			pred = &merlin_p_tp_p13()
	}
	else {					//extended illness-death
		if 	(los==1) 	pred = &merlin_p_tp_p11()
		else if (los==2) 	pred = &merlin_p_tp_p12()
		else if (los==3)	pred = &merlin_p_tp_p113()
		else 			pred = &merlin_p_tp_p123()
	}
	
	return(merlin_p_get_los(gml,t,pred))	
}

`RS' merlin_p_check_tmat(`gml' gml)
{
	if (gml.tmat!=(.,1,2\.,.,3\.,.,.) & gml.tmat!=(.,1,2,.\.,.,.,3\.,.,.,.\.,.,.,.)) {
		errprintf("\n")
		exit(198)
	}
	
	if (gml.tmat==(.,1,2\.,.,3\.,.,.)) 	return(1)
	else 					return(2)
}

`RC' merlin_p_tp_p11(`gml' gml, `RC' t)
{
	gml.model = gml.modtoind = 1
	pred = (*gml.Pch[1])(gml,t)
	gml.model = gml.modtoind = 2
	pred = pred :+ (*gml.Pch[2])(gml,t)
	return(exp(-pred))
}

`RC' merlin_p_tp_p12(`gml' gml, `RC' t)
{
	Ngq 	= 30
	gq 	= merlin_gq(Ngq,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	pred 	= J(gml.N,1,0)
	for (q=1; q<=Ngq; q++) {						
		gml.model = gml.modtoind = 1
		res = (*gml.Plogh[1])(gml,qp[,q]) :- (*gml.Pch[1])(gml,qp[,q])
		gml.model = gml.modtoind = 2
		res = res :- (*gml.Pch[2])(gml,qp[,q])
		res = res :+ log(merlin_p_tp_p22(gml,qp[,q],t))
		pred = pred :+ exp(res :+ log(gq[q,2])) 
	}
	return(pred :* t :/ 2)
}

`RC' merlin_p_tp_p13(`gml' gml, `RC' t) {
	return(1 :- merlin_p_tp_p11(gml,t) :- merlin_p_tp_p12(gml,t))
}

`RM' merlin_p_tp_p22(`gml' gml,  `RC' r, `RC' t)
{
	gml.model = gml.modtoind = 3
	return(exp(-(*gml.Pch[3])(gml,t,r) :+ (*gml.Pch[3])(gml,r,r)))
}

`RM' merlin_p_tp_p23(`gml' gml,  `RC' r, `RC' t)
{
	return(1:-merlin_p_tp_p22(gml,r,t))
}

`RC' merlin_p_tp_p113(`gml' gml, `RC' t)
{
	Ngq 	= 30
	gq 	= merlin_gq(Ngq,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	pred 	= J(gml.N,1,0)
	for (q=1; q<=Ngq; q++) {						
		gml.model = gml.modtoind = 2
		res = (*gml.Plogh[2])(gml,qp[,q]) :- (*gml.Pch[2])(gml,qp[,q])
		gml.model = gml.modtoind = 1
		res = res :- (*gml.Pch[1])(gml,qp[,q])
		pred = pred :+ exp(res :+ log(gq[q,2]))
	}
	return(pred :* t :/ 2)
}

`RC' merlin_p_tp_p123(`gml' gml, `RC' t)
{
	Ngq 	= 30
	gq 	= merlin_gq(Ngq,"legendre")
	qp	= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	pred 	= J(gml.N,1,0)
	for (q=1; q<=Ngq; q++) {						
		gml.model = gml.modtoind = 1
		res = (*gml.Plogh[1])(gml,qp[,q]) :- (*gml.Pch[1])(gml,qp[,q])
		gml.model = gml.modtoind = 2
		res = res :- (*gml.Pch[2])(gml,qp[,q])
		res = res :+ log(merlin_p_tp_p23(gml,qp[,q],t))
		pred = pred :+ exp(res :+ log(gq[q,2])) 
	}
	return(pred :* t :/ 2)
}

`RC' merlin_p_get_los(`gml' gml, `RC' t, `PS' func)
{
	Ngq 	= 30
	gq 		= merlin_gq(Ngq,"legendre")
	qp		= t :/ 2 :* J(gml.N,1,gq[,1]') :+ t:/2
	pred 	= J(gml.N,1,0)
	for (q=1; q<=Ngq; q++) {						
		pred 	= pred :+ (*func)(gml,qp[,q]) :* gq[q,2] :* t :/ 2
	}
	return(pred)
}
end
