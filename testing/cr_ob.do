//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

use "`drive'/multistate/multistate/data/multistate_example",clear
rename pid id

msset, id(id) states(rfi osi) times(rf os) cr

gen stime = _stop

bys id: gen d1 = _status[1]
bys id: gen d2 = _status[2]

bys id: drop if _n==2

set seed 35867
// keep if runiform()<0.3

merlin 	(stime chemo stime, family(exp, failure(d1)) )		///
		(stime chemo stime, family(exp, failure(d2)) )		///
		, 

predict cif1, cif outcome(1) //at(chemo 1)
predict cif2, cif outcome(2) //at(chemo 1)

// merlin 	(stime chemo , family(rp, df(5) failure(d1)))
// est store m1
// merlin 	(stime chemo , family(rp, df(5) failure(d2)))
// est store m2

// cap drop tvar
// range tvar 0 200 100
// galahad, models(m1 m2) cr at1(chemo 0) timevar(tvar)

		
// stset stime , f(d1)
// stcox chemo
		
merlin 	(stime chemo age , family(cox, failure(d1)) )		///
		(stime chemo age , family(cox, failure(d2)) )		///
		, 		
// merlin 	(stime chemo chemo#rcs(stime, df(1) log), family(cox, failure(d1)) timevar(stime))		///
// 		(stime chemo chemo#rcs(stime, df(1) log), family(cox, failure(d2)) timevar(stime))		///
// 		, 

predict cif3, cif outcome(1) standardise //ci reps(100) //at(chemo 1)
predict cif4, cif outcome(2) //at(chemo 1)
		
// merlin (stime chemo age, family(cox, failure(d1)))       ///
//        (stime chemo age, family(rp, df(3) failure(d2)))	

		
		
// gen tf = f1 + f2
// scatter f1 f2 stime]

cap range tvar 0 400 100

// timer clear
// timer on 1
predict cif1, cif outcome(1) at(age 50 chemo 0) timevar(tvar)
predict cif2, cif outcome(1) at(age 50 chemo 0) timevar(tvar) causes(1 2)

predict cif3, cif outcome(2) at(age 50 chemo 0) timevar(tvar)
predict cif4, cif outcome(2) at(age 50 chemo 0) timevar(tvar) causes(1 2)

predict rmst1, rmst outcome(1) at(age 50 chemo 0) timevar(tvar)
predict rmst2, rmst outcome(1) at(age 50 chemo 0) timevar(tvar) causes(1 2)

predict rmst3, rmst outcome(2) at(age 50 chemo 0) timevar(tvar)
predict rmst4, rmst outcome(2) at(age 50 chemo 0) timevar(tvar) causes(1 2)

predict tl1, timelost outcome(1) at(age 50 chemo 0) timevar(tvar)
predict tl2, timelost outcome(1) at(age 50 chemo 0) timevar(tvar) causes(1 2)

predict tl3, timelost outcome(2) at(age 50 chemo 0) timevar(tvar)
predict tl4, timelost outcome(2) at(age 50 chemo 0) timevar(tvar) causes(1 2)

predict ttl1, totaltimelost outcome(1) at(age 50 chemo 0) timevar(tvar)
predict ttl2, totaltimelost outcome(2) at(age 50 chemo 0) timevar(tvar) causes(1 2)


// predict rmst1, rmst fixedonly outcome(1) at(age 50 chemo 0) timevar(tvar)
// predict rmst2, rmst fixedonly outcome(2) at(age 50 chemo 0) timevar(tvar)
// timer off 1

// // quietly {
// stset stime, f(d1)
// stpm2 chemo age, scale(h) df(3)
// est store m1
// stset stime, f(d2)
// stpm2 chemo age, scale(h) df(3)
// est store m2
// // }

// timer on 2
// predictms , models(m1 m2) cr at1(age 50 chemo 0) timevar(tvar) n(100000) los
// timer off 2
// timer list
// twoway scatter _prob_at1_1_2 cif1 tvar, name(g1,replace) nodraw
// twoway scatter _prob_at1_1_3 cif2 tvar, name(g2,replace) nodraw
// twoway scatter _los_at1_1_1 rmst1 tvar, name(g3,replace) nodraw
// graph combine g1 g2 g3 
