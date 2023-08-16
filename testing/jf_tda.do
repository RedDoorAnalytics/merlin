//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

do ./build/buildmlib.do
mata mata clear

use "./data/jointfrailty_example",clear

tab sex, gen(gender)
tab dukes, gen(dukes)
tab charlson, gen(ch)

merlin 	(time  gender1 /*dukes2 dukes3 ch2 ch3*/     						///
						, family(rp, df(3) failure(event)) timevar(time)) 						///
		(stime gender1 /*dukes2 dukes3 ch2 ch3*/  						 	///
						, family(rp, df(2) failure(death)) timevar(stime))

merlin 	(time  gender1 /*dukes2 dukes3 ch2 ch3*/ M1[id]@1     						///
						, family(rp, df(3) failure(event))) 						///
		(stime gender1 /*dukes2 dukes3 ch2 ch3*/ M1[id] 						 	///
						, family(rp, df(2) failure(death)) timevar(stime))

merlin 	(time  gender1 dukes2 dukes3 ch2 ch3 M1[id]@1 M2[id]@1 	///
			, family(rp, df(3) failure(event))) 							///
		(stime gender1 dukes2 dukes3 ch2 ch3 M2[id]@1 		///
			, family(rp, df(1) failure(death)))
			
predict refs*, reffects
predict serefs*, reses
						
// est store m1
range tvar 0 2000 100
predict s1, survival outcome(2) timevar(tvar)
						
						
// merlin 	(time  gender1 dukes2 dukes3 ch2 ch3 M1[id]@1     						///
// 						, family(rcs, df(3) failure(event))) 						///
// 		(stime gender1 dukes2 dukes3 ch2 ch3 M1[id] M1[id]#fp(stime, pow(0)) 	///
// 						, family(rcs, df(3) failure(death)) timevar(stime))
// est store m2
