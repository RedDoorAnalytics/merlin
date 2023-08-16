//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 72549

pr drop _all
clear
set obs 250
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal(0,0.1)
gen trt = runiform()>0.5

survsim stime died, hazard(	0.1:*1.2:*#t:^0.2 :* 					///
							exp(0.2 :* (u1 :+ (0.1:+u2):*#t))) 		///
							cov(trt -0.5) maxt(5)

expand 5
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u1 + (0.1+u2)*time
gen y = rnormal(xb,0.5)

gen start = time
bys id: gen stop = start[_n+1]
gen event = 0
bys id : replace event = died if _n==_N
bys id: replace stop = stime if _n==_N
stset stop, enter(start) id(id) f(event)


bys id1: gen t0 = runiform() if _n==1
bys id1: replace t0 = t0[1]

drop if time<t0
bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1
bys id1 (time) : replace t0 = . if _n>1

timer clear
timer on 2
merlin	(stime trt EV[y], family(weib, failure(died) ltruncated(t0)) timevar(stime)) ///
		(y rcs(time, df(1))  M1[id1]@1,family(gaussian) timevar(time)) ///
		, adaptopts(log)
timer off 2
timer list
	
	
