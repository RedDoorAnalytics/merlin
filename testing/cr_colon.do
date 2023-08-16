//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear


use ./data/colon, clear
 
gen female = (sex==2)

// Expand data
// This creates 2 copies of each observation
set seed 978
// keep if runiform()<0.1
expand 2

// Recode and set up data for competing risk analysis
bysort id: gen cause=_n		// cause =1 for cause 1, cause =2 for cause 2

gen cancer=(cause==1)		// indicator for observation for cancer
gen other=(cause==2)		// indicator for observation for other

// Create dummy variables for age
tab agegrp, gen(ag)

// Categorize age and create interactions with cause
forvalues i = 0/3 {
	gen age`i'can=(agegrp==`i' & cancer==1) 
	gen age`i'oth=(agegrp==`i' & other==1) 
}

// Allow different effect of sex for cancer and other */
gen fem_can = female*cancer
gen fem_other = female*other

// Event indicator
gen event=(cause==status)  // status=1 death due to cancer, =2 death due to other

stset surv_mm, failure(event) scale(12) exit(time 120.5)

timer clear
timer on 1
merlin 	(_t 	female ag2 ag3 ag4									///
				female#rcs(_t , df(3) log event orthog)				///
				ag2#rcs(_t , df(3) log event orthog)				///
				ag3#rcs(_t , df(3) log event orthog)				///
				ag4#rcs(_t , df(3) log event orthog)				///
				if cause==1											///
				, family(rp, failure(event) df(4)) timevar(_t))		///
		(_t 	female ag1 ag2 ag3									///
				if cause==2											///
				, family(rp, failure(event) df(4)))					///
				, 
timer off 1

// timer on 2
// merlin 	(_t 	female ag1 ag2 ag3									///
// 				/*female#rcs(_t , df(3) log event orthog)				///
// 				ag1#rcs(_t , df(3) log event orthog)				///
// 				ag2#rcs(_t , df(3) log event orthog)				///
// 				ag3#rcs(_t , df(3) log event orthog)*/				///
// 				if cause==1											///
// 				, family(rp, failure(event) df(4)))		///
// 		(_t 	female ag1 ag2 ag3									///
// 				if cause==2											///
// 				, family(rp, failure(event) df(4)))					///
// 				, evaltype(gf1)
// timer off 2

timer on 2
merlin 	(_t 	female ag2 ag3 ag4									///
				female#rcs(_t , df(3) log event orthog)				///
				ag2#rcs(_t , df(3) log event orthog)				///
				ag3#rcs(_t , df(3) log event orthog)				///
				ag4#rcs(_t , df(3) log event orthog)				///
				if cause==1											///
				, family(rp, failure(event) df(4)) timevar(_t))		///
		(_t 	female ag2 ag3 ag4									///
				if cause==2											///
				, family(rp, failure(event) df(4)))					///
				, evaltype(gf2) //hessian
timer off 2

// Fit a separate model for cancer and store the knot locations
qui stpm2 fem_can age1can age2can age3can if cancer == 1, ///
	df(4) scale(hazard) dftvc(3) tvc(fem_can age1can age2can age3can) eform nolog   
global knots_cancer `e(bhknots)'
global knots_cancer_tvc `e(tvcknots_age1can)'

// Fit a separate model for other and store the knot locations
qui stpm2 fem_oth age1oth age2oth age3oth if other == 1, ///
	df(4) scale(hazard) eform nolog   
global knots_other `e(bhknots)'


// timer on 4
// // Fit a single model using saved knot locations
// stpm2 cancer other fem_can fem_oth age1can age2can age3can age1oth age2oth age3oth ///
// 	, scale(hazard) rcsbaseoff nocons ///
// 	tvc(cancer other /*fem_can age1can age2can age3can*/) ///
// 	knotstvc(cancer $knots_cancer other $knots_other ///
// 	/*fem_can $knots_cancer_tvc ///
// 	age1can $knots_cancer_tvc ///
// 	age2can $knots_cancer_tvc ///
// 	age3can $knots_cancer_tvc*/)

// timer off 4
timer list
