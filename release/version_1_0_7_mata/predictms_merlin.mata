
/*
stand not done yet

*/

version 14.2

local gml 	struct merlin_struct scalar
local SS 	string scalar
local RS	real scalar
local RM	numeric matrix
local RC	real colvector
local TR	transmorphic
local RR	real rowvector
local PS	pointer(struct merlin_struct scalar) scalar

mata:

`PS' predictms_merlin_setup(`RC' b, `RS' at, `RS' N, `RS' trans, | `RC' t0)
{
	//start by recalling merlin 
	
	stata("qui preserve")
	stata("cap set obs "+strofreal(N))
	stata("cap estimates restore "+st_local("modelests"+strofreal(trans)))

	//start by replacing allvars with zeroes
	allvars = tokens(st_global("e(allvars)"))
	for (i=1;i<=cols(allvars);i++) {
		stata("qui replace "+allvars[i]+" = 0 if _n<="+strofreal(N))
	}
	
	//need to replace any variables with their at#()
	ats = tokens(st_local("at"+strofreal(at)))'
	Natstodo = rows(ats)/2
	i = 1
	j = 2
	for (a=1;a<=Natstodo;a++) {
		stata("qui replace "+ats[i]+" = "+ats[j]+" if _n<="+strofreal(N))
		i = i+2
		j = j+2
	}
	
	//if I replace ltruncated() timevar with newly simulated entry time
	//--> everthing will work for multiple timescales
	//--> this must be after zeros and at, so it's replaced appropriately
	if (args()==5 & st_global("e(ltruncated1)")!="") {
		st_view(lt=.,.,st_global("e(ltruncated1)"))
		lt[|1,1\N,1|] = t0									//!! need to check any sorting doesn't screw this up
	}
	
	//post coefficients
	stata("tempname bmat")
	st_matrix(st_local("bmat"),b')
	
	//remove any options
	stata("local cmd "+st_global("e(cmdline)"))
	stata("gettoken merlin cmd : cmd")
	stata(`"gettoken cmd rhs : cmd, parse(",") bind"')
	
	//recall merlin
	stata("tempname tousem GML")
	
	cmd = "quietly"
// 	cmd = ""
	cmd = cmd + " merlin_parse"
	cmd = cmd + " " + st_local("GML")
	cmd = cmd + " " + ", touse("+st_local("tousem")+") : " + st_local("cmd")
	cmd = cmd + " " + ", predict predictms nogen"
	cmd = cmd + " " + "from("+st_local("bmat")+")"
	cmd = cmd + " " + "npredict("+strofreal(N)+")"
	cmd = cmd + " " + "devcode1("+st_local("devcode1")+")"		//weights
	cmd = cmd + " " + "devcode2("+st_local("devcode2")+")"		//merlin with predictms
	cmd = cmd + " " + "devcode3("+st_local("devcode3")+")"		//prp
	
	//integration method
	cmd = cmd + " " + "intmethod(ghermite)"
	
	
	stata(cmd)
	
// 	quietly `noisily' merlin_parse `GML'	, touse(`tousem') : `cmd'		///
// 										, 								///
// 											predict 					///
// 											nogen 						///
// 											from(`bmat') 				///
// 											`intmethods' 				///
// 											`intpoints' 				///
// 											`ptvar'
	
	//get struct
	struct merlin_struct gml
	gml = *findexternal(st_local("GML"))
	gml.fixedonly = st_local("marginal")==""
	
	stata("qui restore")
	
	return(&gml)
}

/*
for simulation
*/

function merlin_sim(	`RC' x,			/// -simulated survival time-
						`RC' U,			///	-random draw from U(0,1)- 
						`RS' x1,		///
						`RS' x2,		/// - fudged
						`RS' x3,		///
						`RS' x4,		///
						`PS' p, |		///
						`RC' t0			///	-entry time-
						)
{

	//get struct
	`gml' gml
	gml = *p
	mod = gml.model
	nodelent = args()==7
	
	//update N, Nobs and index -> used in merlin utils, and main xb()
	nobs 			= rows(x)
	gml.N = gml.Nobs[mod] = nobs												//!!model number
	asarray(gml.xbindex,mod,asarray(gml.xbindex,mod)[|1,1\nobs,1|])

	//call CH
	if (gml.NI[mod]) {															//!!model number

		Nq 	= gml.hazNnodes[mod]
		qd 	= predictms_gq(Nq)

		if (nodelent) 	nodes = x:/2 :* J(nobs,1,qd[,1]') :+ x:/2
		else  			nodes = (x:-t0):/2 :* J(nobs,1,qd[,1]') :+ (x:+t0):/2

		cumhaz = J(nobs,Nq,0)
		if (gml.familys[mod]=="exp") {
			for (q=1;q<=Nq;q++) {
				cumhaz[,q] = merlin_p_exp_h(gml,nodes[,q])
			}
		}
		else if (gml.familys[mod]=="weibull") {
			for (q=1;q<=Nq;q++) {
				cumhaz[,q] = merlin_p_weibull_h(gml,nodes[,q])
			}
		}
		else if (gml.familys[mod]=="gompertz") {
			for (q=1;q<=Nq;q++) {
				cumhaz[,q] = merlin_p_gompertz_h(gml,nodes[,q])
			}
		}
		else if (gml.familys[mod]=="rcs") {
			for (q=1;q<=Nq;q++) {
				cumhaz[,q] = merlin_p_rcs_h(gml,nodes[,q])
			}
		}
		else if (gml.familys[mod]=="user") {
			for (q=1;q<=Nq;q++) {
				cumhaz[,q] = merlin_p_userh_h(gml,nodes[,q])
			}
		}
		
		if (nodelent) 	ret1 = x:/2 :* cumhaz * qd[,2]		
		else 			ret1 = (x:-t0):/2 :* cumhaz * qd[,2]		
		
	}
	else {

		if (gml.familys[mod]=="rp") 			ret1 = merlin_p_rp_ch(gml,x)
		else if (gml.familys[mod]=="exp") 		ret1 = merlin_p_exp_ch(gml,x)
		else if (gml.familys[mod]=="weibull") 	ret1 = merlin_p_weibull_ch(gml,x)
		else if (gml.familys[mod]=="gompertz") 	ret1 = merlin_p_gompertz_ch(gml,x)
		else if (gml.familys[mod]=="user") 		ret1 = merlin_p_userh_ch(gml,x)

		if (!nodelent) {
			index = 1::nobs
			index = select(index,t0:>0)
			nobs2 = rows(index)
			if (nobs2) {
				
				gml.Nobs[mod]	= nobs2															
				asarray(gml.xbindex,mod,asarray(gml.xbindex,mod)[|1,1\nobs2,1|])
				
				if (gml.familys[mod]=="rp") {
					ret1[index,] = ret1[index,] :- merlin_p_rp_ch(gml,t0[index])
				}
				else if (gml.familys[mod]=="exp") {
					ret1[index,] = ret1[index,] :- merlin_p_exp_ch(gml,t0[index])
				}
				else if (gml.familys[mod]=="weibull") {
					ret1[index,] = ret1[index,] :- merlin_p_weibull_ch(gml,t0[index])
				}
				else if (gml.familys[mod]=="gompertz") {
					ret1[index,] = ret1[index,] :- merlin_p_gompertz_ch(gml,t0[index])
				}
				else if (gml.familys[mod]=="user") {
					ret1[index,] = ret1[index,] :- merlin_p_userh_ch(gml,t0[index])
				}

			}
			
		}
		
	}
	
	return(ret1 :+ log(U))
}

/*
	for AJ
*/

`RM' merlin_ch(	`RC' x,			///	-time-
				`PS' p			///
				)
{

	//get struct
	`gml' gml
	gml = *p
	mod = gml.model

	//update N, Nobs and index -> used in merlin utils, and main xb()
	nobs = rows(x)
	gml.N = gml.Nobs[gml.Nlevels,mod] = nobs												//!!model number
	asarray(gml.xbindex,mod,asarray(gml.xbindex,mod)[|1,1\nobs,1|])

	//call CH
	if 		(gml.familys[mod]=="rp") 		ret1 = merlin_p_rp_ch(gml,x)
	else if (gml.familys[mod]=="rcs") 		ret1 = merlin_p_rcs_ch(gml,x)
	else if (gml.familys[mod]=="exp") 		ret1 = merlin_p_exp_ch(gml,x)
	else if (gml.familys[mod]=="weibull") 	ret1 = merlin_p_weibull_ch(gml,x)
	else if (gml.familys[mod]=="gompertz") 	ret1 = merlin_p_gompertz_ch(gml,x)
	else if (gml.familys[mod]=="user") 		ret1 = merlin_p_userh_ch(gml,x)

	//correct for time = 0
	time0 = selectindex(x:==0)
	if (rows(time0) & cols(time0)) ret1[time0,] = J(rows(time0),cols(ret1),0)
	
	return(ret1)
}


end
