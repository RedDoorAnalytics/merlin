//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/forge"
adopath ++ "./forge"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 500
gen id = _n
gen u1 = rnormal()
expand 5
gen trt = runiform()>0.5

survsim stime dead , loghazard(-2.3 :+ 1:*log(#t) :+ 0.5:*trt :+ u1) maxt(5)
stset stime, f(dead)

mata:
real matrix loghaz(gml,t)
{
	xb 		= forge_util_xzb(gml,t)
	slope  	= forge_util_xb(gml,1)
	return(slope:*log(t) :+ xb)
}
end

forge (stime trt M1[id]@1, family(user, loghfunction(loghaz) nap(1)  failure(dead)))

range tvar 0 5 100

predict s1, survival fixedonly timevar(tvar)
predict s2, survival marginal timevar(tvar)
predict s3, hazard fixedonly timevar(tvar)
predict s4, hazard marginal timevar(tvar)
predict s5, rmst fixedonly timevar(tvar)
predict s6, rmst marginal timevar(tvar)
predict s7, chazard fixedonly timevar(tvar)
predict s8, chazard marginal timevar(tvar)
predict s9, eta fixedonly timevar(tvar)
predict s10, eta marginal timevar(tvar)
predict s11, mu fixedonly timevar(tvar)
predict s12, mu marginal timevar(tvar)
