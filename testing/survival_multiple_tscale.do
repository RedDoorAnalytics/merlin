
//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

do ./build/buildmlib.do
mata mata clear

//local drive Z:\
local drive /Users/Michael/Documents
//local drive c:
cd "`drive'/multistate/multistate"
adopath ++ "."
adopath ++ "./msset"
adopath ++ "./predictms"
adopath ++ "./stms"
adopath ++ "./msboxes"
adopath ++ "./msaj"
adopath ++ "./graphms"
clear all

tr:do ./build/buildmlib.do
mata mata clear

clear
set seed 249857
set obs 10000
gen trt = runiform()>0.5
gen age = rnormal(50,5)
gen agec = age - 50
gen year = 1990 + floor(20*runiform())
gen yearc = year - 2000

survsim stime died, maxtime(5) cov(trt -0.5) 					///
hazard(	0.1:*1.2:*{t}:^0.2 :*				///
        exp(								///
                        0.1 :* (agec :+ {t}) 		///
                        :- 0.1 :* (yearc :+ {t})	///
                )								///
        )

merlin (stime 	trt 																///
        fp(stime, pow(1) offset(agec)) 										///
        fp(stime, pow(1) offset(yearc)) 									///
        , family(weibull, failure(died)) timevar(stime))
est store m1
predictms, surv singleevent models(m1) timevar(stime) at1() devcode1(41bsjdh82e18ndu3)

merlin (stime 	trt 																///
        rcs(stime, df(1) offset(agec)) 										///
        rcs(stime, df(1) offset(yearc)) 									///
        , family(pwexp, knots(1) failure(died)) timevar(stime))
        
merlin (stime 	trt 														///
        rcs(stime, df(1) offset(agec)) 										///
        rcs(stime, df(1) offset(yearc)) 									///
        , family(rp, df(1) failure(died)) timevar(stime)), 

merlin (stime 	trt 														///
        rcs(stime, df(1) offset(agec)) 										///
        rcs(stime, df(1) offset(yearc)) 									///
        , family(cox, failure(died)) timevar(stime)), 

        
timer clear
// timer on 1				
// merlin (stime 	trt 																///
// 				rcs(stime, df(1) offset(agec)) 										///
// 				rcs(stime, df(1) offset(yearc)) 									///
// 				rcs(stime, df(1) offset(agec))#rcs(stime, df(1) offset(yearc))		///
// 				rcs(stime, df(1) log)												///
// 				, family(loghazard, failure(died)) timevar(stime)),					///
// 				evaltype(gf2)
// timer off 1
// timer on 2
// merlin (stime 	trt 																///
// 				rcs(stime, df(1) offset(agec)) 										///
// 				rcs(stime, df(1) offset(yearc)) 									///
// 				rcs(stime, df(1) offset(agec))#rcs(stime, df(1) offset(yearc))		///
// 				, family(cox, failure(died)) timevar(stime)),					///
// 				evaltype(gf2)
// timer off 2
timer on 3
merlin (stime 	trt 																///
        rcs(stime, df(2) offset(agec)) 										///
        rcs(stime, df(2) offset(yearc)) 									///
        /*rcs(stime, df(1) offset(agec))#rcs(stime, df(1) offset(yearc))*/		///
        , family(rp, df(3) failure(died)) timevar(stime)),					///
        evaltype(gf2)
timer off 3
        
stset stime, f(died)

timer on 4
stmt trt, 	time1(df(1)) 							///
time2(df(1) start(agec) logtoff) 		///
time3(df(1) start(yearc) logtoff) 		///
nohr noorthog
timer off 4
timer list
