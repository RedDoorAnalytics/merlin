//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"

clear all
tr:do ./build/buildmlib.do
mata mata clear

clear
webuse jobhistory
stset tend, origin(tstart) fail(failure)

sort birthyear id
order birthyear id, first

mestreg education  || birthyear: || id:, 	///
	distribution(weib) adaptopts(log)
predict refs*, reffects

stmixed education || birthyear: || id: , dist(weib) adaptopts(log)

