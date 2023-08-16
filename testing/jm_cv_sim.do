//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
adopath ++ "./stmerlin"
adopath ++ "./jm"
clear all

do ./build/buildmlib.do
mata mata clear

pr drop _all

set seed 7254
clear
set obs 500
gen id1 = _n
gen u0 = rnormal(0,1)
gen u1 = rnormal(0,0.3)
gen trt = runiform()>0.5

survsim stime died, hazard(0.1:*1.2:*{t}:^0.2 :* 		/// baseline hazard
        exp(0.3 :* (u0 :+ (0.1:+u1):*{t} :+ 0.2:*{t}:^2) 	///
        :- 0.1 :* (u0 :+ (0.1:+u1):*{t} :+ 0.2:*{t}:^2):^2) 	/// exp() of time-dependent linear predictor
        )						        ///
        covariates(trt -0.5) 					///	time-independent linear predictor
        maxt(10)						//	admin. censoring

expand 10
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u0 + (0.1+u1) * time + 0.2 * time^2
gen y = rnormal(xb,0.5)

bys id1 (time) : replace stime 	= . if _n>1
bys id1 (time) : replace died 	= . if _n>1

// gsem (y <- time M1[id]@1, family(gauss))
// cap drop test*
// predictnl test = predict(mu cond(ebmeans)), ci(tfest_lci test_uci) force

// jm 	(y fp(time, pow(1 2)) || id: time, family(gaussian) timevar(time))		///
// 	(stime trt EV[y], family(weibull, failure(died)) timevar(stime))
mata:
real matrix userelement(gml,t)
{
element = merlin_util_expval_mod(gml,1,t)
return(element)
}
end
	
	
merlin 	(y time fp(time, pow(2)) time#M2[id]@1 M1[id]@1, family(gaussian) timevar(time))	///
		(stime trt mf(userelement), family(w, failure(died)) timevar(stime))	        ///
		, cov(unstr)

// predict refs2*, reffects
// predict s0, surv fixedonly outcome(2)
// predict s1, surv fixedonly outcome(2) ci
// predict s2, surv fitted outcome(2)
range tvar1 0 5 20
predict s32, surv fitted outcome(2) panel(1) timevar(tvar1) devcode7(bc928r72crncfye98fyqc9r398ry)
predict s33, surv fitted outcome(2) panel(1) timevar(tvar1) //devcode7(bc928r72crncfye98fyqc9r398ry)

// //longitudinal
// range tvar1 0 5 20
// predict l1, mu outcome(1) fixedonly ci
// predict l2, mu outcome(1) fitted ci


forvalues i=2/5 {
	predict fel`i' if _n<=`i', mu outcome(1) fixedonly 
	predict l`i' if _n<=`i', mu outcome(1) fitted panel(1) ci
	gen lt`i' = `i'-1 if _n<=20
	range t`i' `i'-1 10 20
	predict s`i', blupif(_n<=`i') surv outcome(2) fitted panel(1) ltruncated(lt`i') timevar(t`i') ci at(trt 1)
	twoway 	(rarea l`i'_lci l`i'_uci time if _n<=`i', )(line l`i' time if _n<=`i')(scatter y time if _n<=`i')	///
                (rarea s`i'_lci s`i'_uci t`i', yaxis(2) color(ltblue))(line s`i' t`i', yaxis(2) color(gray)) 		///
                (function y=1, yaxis(2) range(0 `=`i'-1') col(gray))		///
                , 	name(g`i', replace)				///
                        legend(order(2 "Predicted long." 3 "Observed long." 5 "Predicted survival") cols(1))	///
                        ylabel(0(0.2)1, axis(2) format(%2.1f) angle(h)) ///
                        ylabel(-8(2)8, angle(h) axis(1)) 				///
                        ytitle("Survival",axis(2) orient(rvertical)) ytitle("Biomarker",axis(1)) 									///
                        title("Patient-specific conditional survival") ///
                        xtitle("Follow-up time") xlabel(0(2)12)
graph export g`i'.png, name(g`i') replace
// 	local all `all' g`i'
}
// grc1leg `all', rows(1) xcommon ycommon 



// local graph1 (rarea l1_lci l1_uci time if id==1)(line l1 time if id==1)
// local graph1 `graph1' (scatter y time if id==1)
// // twoway `graph1', name(g1, replace)
// local graph2 (rarea l2_lci l2_uci time if id==1)(line l2 time if id==1)
// local graph2 `graph2' (scatter y time if id==1)
// // twoway `graph2', name(g2, replace)

// //survival
// range tvar2 4 10 20
// gen t0var2 = 4 in 1/20
// predict cs1, surv outcome(2) fitted panel(1) timevar(tvar2) ltruncated(t0var2) at(trt 1)
// predict cs2, surv outcome(2) fitted panel(1) timevar(tvar2) ltruncated(t0var2) at(trt 1) ci

// local graph4 (rarea cs2_lci cs2_uci tvar2, yaxis(2))(line cs2 tvar2, yaxis(2))

// twoway `graph2' `graph4', name(cs1,replace)

twoway function y=3 + 0.4*x + 1*x^2, range(0 10)



