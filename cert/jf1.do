use https://www.mjcrowther.co.uk/data/jointfrailty_example.dta, clear

tab sex, gen(gender)
tab dukes, gen(dukes)
tab charlson, gen(ch)

merlin 	(time  gender1 dukes2 dukes3 ch2 ch3 M1[id]@1 ///
			, family(rp, df(3) failure(event))) 					///
		(stime gender1 dukes2 dukes3 ch2 ch3 M1[id] ///
			, family(rp, df(1) failure(death)))
do ./cert/predictions.do

merlin 	(time  gender1 dukes2 dukes3 ch2 ch3 M1[id]@1 M2[id]@1 	///
			, family(rp, df(3) failure(event))) 							///
		(stime gender1 dukes2 dukes3 ch2 ch3 M2[id]@1 		///
			, family(rp, df(1) failure(death)))
do ./cert/predictions.do


// use https://www.mjcrowther.co.uk/data/jointfrailty_example.dta, clear
//
// tab sex, gen(gender)
// tab dukes, gen(dukes)
// tab charlson, gen(ch)
//
// merlin 	(tstop  M1[id]@1 ///
// 			, family(rp, df(3) failure(event) ltruncated(tstart))) 	///
// 	(stime  M1[id] ///
// 			, family(rp, df(1) failure(death))), intpoints(15)
// do ./cert/predictions.do
//
// merlin 	(tstop  gender1 dukes2 dukes3 ch2 ch3 M1[id]@1 M2[id]@1 	///
// 			, family(rp, df(3) failure(event) ltruncated(tstart))) 	///
// 		(stime gender1 dukes2 dukes3 ch2 ch3 M2[id]@1 		///
// 			, family(rp, df(1) failure(death))), intpoints(15)
// do ./cert/predictions.do
//
//
