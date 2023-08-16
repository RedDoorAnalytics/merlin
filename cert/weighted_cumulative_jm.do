
use https://www.mjcrowther.co.uk/data/jm_example.dta, clear

merlin 	(logb rcs(time,df(3)) fp(time,pow(1))#M2[id]@1 M1[id]@1	, family(gaussian) timevar(time)) 	///
		(stime trt EV[logb]							, family(weibull, failure(died)) timevar(stime))	///
		, covariance(unstructured) restartvalues(M2 0.01)
do ./cert/predictions.do

merlin 	(logb rcs(time,df(3)) fp(time,pow(1))#M2[id]@1 M1[id]@1	, family(gaussian) timevar(time)) 	///
		(stime trt iEV[logb]						, family(weibull, failure(died)) timevar(stime))	///
		, covariance(unstructured) restartvalues(M2 0.01)
do ./cert/predictions.do

//weighted cumulative 
mata mata clear
mata:
//define my function, which needs the second input of time and we will be integrating
//it out to get the cumulative hazard required in the likelihood			
real matrix userhaz(gml, t)
{
	//declare everything we need
	real matrix	linpred, logh, gq, nodes, weights, intres
	real scalar 	gam1, alpha, sd, maxtime, nq, Nobs

	//evaluate the main complex linear predictor log(lambda) of the Weibull model
	linpred = merlin_util_xzb(gml,t)
	
	//evaluate the anciliary parameters, ensuring correct scale
	gam1 	= exp(merlin_util_ap(gml,1))			//gamma of the Weibull model
	alpha 	= merlin_util_ap(gml,2)					//association parameter
	sd	= exp(merlin_util_ap(gml,3))			//scale parameter - standard deviation
	sd2	= sd^2									//variance
	
	//setup for numerical integration
	nq 	= 30
	gq 	= merlin_gq(nq,"legendre")			//extract basis nodes and weights
	Nobs	= rows(t)
	nodes 	= t :* J(Nobs,1,gq[,1]'):/2 :+ t:/2
	weights = t :* J(Nobs,1,gq[,2]'):/2
	
	//evaluate denominator of weight function -> constant 
	maxtime	= max(merlin_util_depvar(gml)[,1])		
	mvnmaxt = mvnormalcv(0,maxtime,0,sd2)
	
	//conduct integration
	//-> merlin_util_xzb_mod() evaluates the linear predictor of the biomarker, at specified timepoints
	intres 	= J(Nobs,1,0)
	for (q=1;q<=nq;q++) {
		intres = intres :+ weights[,q] :* 
                        normalden(t:-nodes[,q],sd) :* 
                        merlin_util_xzb_mod(gml,1,nodes[,q]) 
	}
	intres = intres :/ mvnmaxt
	
	//build and return hazard function
	logh = linpred :+ log(gam1) :+ (gam1:-1):*log(t)  :+ alpha :* intres
	return(logh)
}
end

mat svs = e(b)
mat svs2 = svs[1,1..7]		//logb
mat svs2 = svs2,svs[1,8]	//trt
mat svs2 = svs2,svs[1,10]	//intercept
mat svs2 = svs2,svs[1,9]	//alpha
mat svs2 = svs2,svs[1,11]	//loggamma
mat svs2 = svs2,-2		//sd
mat svs2 = svs2,svs[1,12..14]

merlin 	(logb rcs(time,df(3)) fp(time,pow(1))#M2[id]@1 M1[id]@1	, family(gaussian) timevar(time)) 		///
		(stime trt 				, family(user, loghfunc(userhaz) nap(3) failure(died)) timevar(stime))	///
		, covariance(unstructured) intmethod(gh) from(svs2)

do ./cert/predictions.do			
		
