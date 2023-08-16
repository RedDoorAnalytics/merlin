//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 2000
gen id1 = _n
gen u1 = rnormal()
gen u2 = rnormal(0,0.1)
gen trt = runiform()>0.5

survsim stime died, hazard(0.1:*1.2:*#t:^0.2 :* 					///
							exp(									///
			3 :* ( (u1 :+ (0.2:+u2):*#t) :- (u1 :+ (0.2:+u2):*(#t:-1))) :*(#t:>=1) 		///
								)) cov(trt -0.5) maxt(5)

expand 5
bys id1 : gen time = _n-1
drop if time>stime
drop if stime<1

gen xb = u1 + (0.2+u2)*time
gen y = rnormal(xb,0.5)

gen start = time
bys id: gen stop = start[_n+1]
gen event = 0
bys id : replace event = died if _n==_N
bys id: replace stop = stime if _n==_N
stset stop, enter(start) id(id) f(event)

bys id1 (time) : replace stime = . if _n>1
bys id1 (time) : replace died = . if _n>1
bys id1 (time) : gen t0 = 1 if _n==1

mata:

real matrix merlin_cv_diff(gml,t)
{
	diff = merlin_util_xzb_mod(gml,2,t) :- merlin_util_xzb_mod(gml,2,t:-1)	
	return(diff)
}
end

mata:
real matrix merlin_cv_diff3(gml,t)
{
	diff = merlin_util_xzb_mod(gml,2,t) :- merlin_util_xzb_mod(gml,2,t:-exp(merlin_util_ap(gml,1)))	
	return(diff)
}
end

merlin	(stime trt EV[y] , family(weib, failure(died) ltruncated(t0)) timevar(stime)) 		///
		(y fp(time,pow(1)) fp(time, pow(1))#M2[id1]@1 M1[id1]@1,family(gaussian) timevar(time)) 		///
			,  devcode1(385927) 
mat sv = e(b)

merlin	(stime trt mf(merlin_cv_diff) , family(weib, failure(died) ltruncated(t0)) timevar(stime)) 		///
		(y fp(time,pow(1)) fp(time, pow(1))#M2[id1]@1 M1[id1]@1,family(gaussian) timevar(time)) 		///
			,  devcode1(385927) from(sv)

mat sv = e(b)
mat sv = sv[1,1..4],0,sv[1,5..11]
merlin	(stime trt mf(merlin_cv_diff3) , family(weib, failure(died) ltruncated(t0) nap(1)) timevar(stime)) 		///
		(y fp(time,pow(1)) fp(time, pow(1))#M2[id1]@1 M1[id1]@1,family(gaussian) timevar(time)) 		///
			,  devcode1(385927) from(sv)	
	
