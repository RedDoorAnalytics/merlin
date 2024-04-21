//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"

clear all
tr:do ./build/buildmlib.do
mata mata clear

webuse jobhistory
stset tend, origin(tstart) fail(failure)
sort birthyear id

mestreg education njobs 		///
	|| birthyear: || id:		///
	, distribution(weibull) nohr //intmethod(gh) intpoints(25)
predict refs*, reffects

// stmixed education njobs 		///
// 	|| birthyear: || id:		///
// 	, distribution(weibull)	showmerlin
//	
	
merlin 	(_t education njobs 		///
	M2[birthyear>id]@1 		///
	M1[birthyear]@1 		///
	, family(weib, failure(_d))) 	///
	, intmethod(mvagh) intpoints(7) 

predict rfs*, reffects debug
