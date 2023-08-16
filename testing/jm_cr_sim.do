//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

clear
set seed 7254
set obs 500
gen id1 = _n
gen u1 = rnormal(0,1)
gen u2 = rnormal(0,0.5)
gen trt = runiform()>0.5

local maxt 10

//cause one
local l1 0.1
local g1 1.2
survsim stime1 d1, hazard(`l1':*`g1':*{t}:^(`g1'-1) :* exp( 0.1 :* (u1 :+ (0.1:+u2):*{t} :+ 0.1:*{t}:^2))) cov(trt -0.2) maxt(`maxt')

//cause two
local l2 0.05
local g2 1.2
survsim stime2 d2, hazard(`l2':*`g2':*{t}:^(`g2'-1) :* exp(0.05 :* (u1 :+ (0.1:+u2):*{t} :+ 0.1:*{t}:^2))) cov(trt -0.5) maxt(`maxt')
								
								
//find the minimum
egen double stime = rowmin(stime1 stime2)

replace d1 = 0
replace d1 = 1 if stime1==stime & stime!=`maxt'
replace d2 = 0
replace d2 = 1 if stime2==stime & stime!=`maxt'

expand 10
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u1 + (0.1+u2)*time + 0.1*time^2
gen y = rnormal(xb,0.5)

bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace d1 = . if _n>1
bys id1 (time) : replace d2 = . if _n>1

// merlin 	(y fp(time, pow(1 2)) fp(time, pow(1))#M2[id]@1 M1[id]@1, family(gaussian) timevar(time))	///
// 		(stime trt EV[y] , family(weibull, failure(d1)) timevar(stime))							///
// 		(stime trt EV[y], family(weibull, failure(d2)) timevar(stime))							///
// 		, 
// est save jmcr, replace
est use jmcr
		
// predict cif1, cif marginal outcome(2) at(trt 0) timevar(tvar)
// predict cif2, cif marginal outcome(2) at(trt 1) timevar(tvar)
// predict cif3, cif marginal outcome(3) at(trt 0) timevar(tvar)
// predict cif4, cif marginal outcome(3) at(trt 1) timevar(tvar)

// gen totalcif1 = cif1 + cif3
// gen totalcif2 = cif2 + cif4

// twoway (area totalcif1 tvar)(area cif1 tvar), name(g1,replace)
// twoway (area totalcif2 tvar)(area cif3 tvar), name(g2,replace)
		
// predictnl double difcif = predict(cif marginal outcome(2) at(trt 1) timevar(tvar)) 	///
// 						- predict(cif marginal outcome(2) at(trt 0) timevar(tvar))	///
// 						, ci(difcif_lci difcif_uci)
						
// twoway (rarea difcif_lci difcif_uci tvar)(line difcif tvar), name(g3,replace)
predict refs*, reffects
keep if id>=4
forvalues i=2(2)10 {
	predict l`i' if _n<=`i', mu outcome(1) fitted panel(4) at(trt 0)
	gen lt`i' = `=time[`i']' if _n<=20
	range t`i' `=time[`i']' 12 20
	predict cif1`i', blupif(_n<=`i') cif outcome(2) fitted panel(4) ltruncated(lt`i') timevar(t`i') at(trt 0)
	predict cif2`i', blupif(_n<=`i') cif outcome(3) fitted panel(4) ltruncated(lt`i') timevar(t`i') at(trt 0)
	gen totalcif`i' = cif1`i' + cif2`i'
	twoway 	(line l`i' time if _n<=`i')(scatter y time if _n<=`i')	///
			(area totalcif`i' t`i', yaxis(2) color(ltblue)) 		///
			(area cif2`i' t`i', yaxis(2) color(midblue)) 		///
			(function y = `=time[`i']', horizontal range(-1 12) lpattern(dash)) ///
			, 	name(g`i', replace)				///
				legend(order(4 "Prob(Death due to cancer)" 5 "Prob(Death due to other causes)") cols(1))	///
				ylabel(0(0.2)1, axis(2) format(%2.1f) angle(h)) ///
				ylabel(-1 "0" 0 "1" 1 "2" 2 "3" 3 "4" 4 "5" 5 "6" 6 "7" 7 "8" 8 "9" 9 "10" 10 "11" 11 "12" 12 "13", ///
				format(%2.1f) angle(h) axis(1)) 				///
				ytitle("Mortality",axis(2) orient(rvertical)) ytitle("Biomarker",axis(1)) 									///
				title("Dynamic conditional cause-specific CIFs") ///
				xtitle("Follow-up time") xlabel(0(2)12) 
graph display, xsize(20) ysize(10) //name(g`i',replace)
graph export cif`i'.png, name(g`i') replace

}


		
