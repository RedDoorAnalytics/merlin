
set seed 874
clear
set obs 300
gen id = _n

//treatment group with 0.5 prob of treatment vs placebo
gen trt = rbinomial(1,0.5)

//random intercept for each of the two biomarkers
gen b01 = rnormal(0,1)
gen b02 = rnormal(0,1)

//simulate survival times under current value assoc for both biomarkers
//exp baseline hazard with lambda = 0.1
//admin censoring at 5 years
//treatment effect of HR = 0.607
survsim stime died, hazard(	0.1 :* exp(							///	
								0.1 :* (b01 :- 0.5 :* #t) :+	///
								0.2 :* (b02 :+ 0.5 :* #t)		///
							)) cov(trt -0.5) maxtime(5)

// up to 5 observations each, at timepoints 0,1,2,3,4
expand 5
sort id
by id : gen time = _n - 1
//remove any after event/censoring time
drop if time>stime

//need only one survival time per patient (must do this!)
bys id: replace stime = . if _n>1
bys id: replace died = . if _n>1

//generate observed biomarkers, with residual error
gen y1 = b01 - 0.5 *time + rnormal(0,0.5)
gen y2 = b02 + 0.5 *time + rnormal(0,0.5)

//fit true model
merlin  (y1  fp(time, pow(1)) M1[id]@1  , family(gaussian) timevar(time))   ///
        (y2  fp(time, pow(1)) M2[id]@1  , family(gaussian) timevar(time))   ///
        (stime trt EV[y1] EV[y2]  		, family(weibull, failure(died)) timevar(stime))   ///
        , covariance(unstructured) 
