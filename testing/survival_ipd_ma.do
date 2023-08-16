//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
clear all

do ./build/buildmlib.do
mata mata clear

use ./data/SurvivalIPD,clear

stset stime, failure(event==1)


tab trial, gen(trialvar)
stcox treat


stmixed  treat  || trial: treat ,  ///
        distribution(rp) df(3) cov(uns) 

