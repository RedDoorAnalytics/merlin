
clear

set seed 934431

local N = 500

local a0 = 0
local a1 = 1
local a2 = 0.2

local sigmau2 = 1
local sigmav2 = 0.5
local sigmae2 = 1

local beta = 1
local eta = 1

local g1 = -1.5
local g2 = -0.5
local g3 = 1

set obs `N'
gen id = _n
gen trt = runiform()>0.5

gen double u = rnormal(0,`=sqrt(`sigmau2')')
gen double v = rnormal(0,`=sqrt(`sigmav2')')

//survival time
survsim stime died, dist(weibull) lambda(0.02) gamma(2) cov(trt `eta' u `g2' v `g3') maxt(12)

//recurrent event/observation times

gen double t0 = 0
local index = 0
count if t`index'==.
local NC = `r(N)'
while `NC'!=`N' {
	local index = `index' + 1
	gen double t`index'=.
	gen double _temp_t`index' = t`=`index'-1' + (-log(runiform()) / (0.08 * exp(trt * `beta' + u)))^(1/2) if t`=`index'-1'!=.
	su _temp_t`index'
	replace _temp_t`index' = . if _temp_t`index'>stime & _temp_t`index'!=.
	replace t`index' = _temp_t`index' if _temp_t`index'!=.
	count if t`index'==.
	local NC = `r(N)'
}

reshape long t, i(id) j(time)
drop if t==.  
drop _temp*
drop time
//generate longitudinal outcome

gen double xb = `a0' + trt*`a1' + t*`a2' + `g1'*u + v
gen double y = rnormal(xb,`=sqrt(`sigmae2')')

//for megenreg, need one row per id for survival model
bys id: replace stime = . if _n!=_N
bys id: replace died = . if _n!=_N

//get reset times for ob process
bys id: gen double obtime = t[_n+1]-t[_n]

gen ind = 1 if obtime!=.
replace ind=0 if obtime==.
replace obtime = stime-t if obtime==.

//fit true model
/*
gsem 	(stime <- trt M1[id] M2[id], family(weibull, failure(died)))	///
		(obtime <- trt M1[id]@1, family(weibull, failure(ind)))			///
		(y <- trt t M1[id] M2[id]@1, family(gaussian)), 				///
		covstruct(M1[id] M2[id], diag) 
*/		

timer clear
timer on 1
megenreg 	(stime trt M1[id]@g2 M2[id]@g3	, family(weibull, failure(died)))	///
			(obtime trt M1[id]				, family(weibull, failure(ind)))	///
			(y trt t M1[id]@g1 M2[id]		, family(gaussian)),
timer off 1
timer list




