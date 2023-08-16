//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear


set seed 7254
pr drop _all
clear
set obs 500
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal()


gen trt = runiform()>0.5
survsim stime died, hazard(0.1:*1.2:*{t}:^0.2 :* 	///
	exp(0.2:*(u1 :+ {t} :+0.5:*trt :+ 0.1:*(u2 :+ {t} :+0.2:*trt)) ///
	))   cov(trt -0.5) maxt(5)

expand 7
bys id1 : gen time = _n-1

drop if time>stime

gen xb2 = u1 + time + 0.5*trt
gen xb1 = u2 + time + 0.5*trt + 0.1 * xb2
gen y1 = rnormal(xb1,0.5)
gen y2 = rnormal(xb2,0.5)

gen start = time
bys id: gen stop = start[_n+1]
gen event = 0
bys id : replace event = died if _n==_N
bys id: replace stop = stime if _n==_N

stset stop, enter(start) id(id) f(event)

//
//stjm y trt, panel(id) survm(weib) ffp(1 2) gh(7) survcov(trt) 

bys id1 (time) : replace stime = . if _n<_N
bys id1 (time) : replace died = . if _n<_N

timer clear
timer on 1
merlin 	(stime trt EV[y1] , family(weib, failure(died)) timevar(stime))	///
			(y1 trt time M1[id1]@1 EV[y2] , timevar(time) family(gaussian))	 ///
			(y2 trt time M2[id1]@1 , timevar(time) family(gaussian)),	
timer off 1
timer list

// {#t}*M2[id1]
// trt*{#t}	
// trt*{#t}*M2[id1]	
