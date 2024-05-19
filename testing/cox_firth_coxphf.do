//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"


//build mlib
clear all
tr:do ./build/buildmlib.do
mata mata clear

set seed 7254
pr drop _all
clear
set obs 3
input stime died x
1 1 1
2 1 1
3 1 0

merlin (stime x, family(cox, firth failure(died))) 

stset stime, f(died)
stmerlin x, dist(cox) firth 

// stset TIME, f(CENS)
// stmerlin T N G CD, dist(cox) firth 
 