
webuse catheter, clear

mata:
real matrix userhaz(gml,t)
{
	linpred = merlin_util_xzb(gml,t)
	gammas 	= merlin_util_ap(gml,1)\merlin_util_ap(gml,2)
	return(exp(linpred :+ merlin_fp(t,(0,1)) * gammas))
}
end

merlin (time age female, family(user, hfunction(userhaz) failure(infect) nap(2)))
do ./cert/predictions.do

merlin (time age female M1[patient]@1, family(user, hfunction(userhaz) failure(infect) nap(2)))
do ./cert/predictions.do

merlin (time age age#fp(time, pow(0)) female M1[patient]@1, ///
          family(user, hfunc(userhaz) failure(infect) nap(2)) timevar(time)) , zeros diff
do ./cert/predictions.do		  
