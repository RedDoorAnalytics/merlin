//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/megenreg/Stata"
adopath ++ "./megenreg"
adopath ++ "`drive'/V1.3.1"

clear all

do ./build/buildmlib.do
mata mata clear

qui do `drive'/V1.3.1/stmt_matacode.mata"

clear
set seed 249857
set obs 100000
gen trt = runiform()>0.5
gen age = rnormal(50,5)
gen year = 1990 + floor(20*runiform())
gen yearc = year - 2000

survsim stime died , maxt(5)  cov(trt -0.5) ///
	hazard(exp(-2 :+ 0.1 :* log(#t) :+ (age :+ #t) :* 0.1 :- 0.01 :* (yearc :+ #t)))

mata:
real matrix rcs(transmorphic gml, 		///
					real matrix t)	
{
	linpred 	= megenreg_util_xzb(gml)
	xb 			= megenreg_util_xb(gml,1) //\megenreg_util_xb(gml,2)
	//knots 		= (-5,-1,1.609)
	//calculate and return the hazard function
	return(linpred :+ log(t) * xb)
}
end	
	
megenreg (stime trt rcs(df(1) offset(age))@ts1 		///
					rcs(df(1) offset(yearc))@ts2 	///
					, family(user, failure(died) loghfunc(rcs) nap(1)) timevar(stime))

stset stime, f(died)
stmt trt, time1(df(1)) time2(df(1) start(age) logtoff) time3(df(1) start(yearc) logtoff) noorthog nohr
