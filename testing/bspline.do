//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

use https://www.mjcrowther.co.uk/data/jm_example.dta, clear

merlin 	(logb 	bs(time,df(5) d(3))							///
				M1[id]@1								///
				, family(gaussian)) 					///
		(stime trt EV[logb] 							///
				, family(weibull, failure(died)) timevar(stime))		///
		, 
		
predict l1, mu marginal 

// merlin 	(logb 	rcs(time,df(3))							///
// 				M1[id]@1								///
// 				, family(gaussian)) 					///
// 		(stime trt EV[logb] 							///
// 				, family(weibull, failure(died)) timevar(stime))		///
// 		, 
		
// predict l2, mu marginal 

scatter l1 time

