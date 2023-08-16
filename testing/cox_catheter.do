//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

webuse catheter,clear
replace time = time + 0.001*runiform()
stset time, f(infect) scale(365.25)

merlin (_t age female female#rcs(_t, df(1) log), family(cox, failure(_d)) timevar(_t))

predict s1, survival at(female 1)  ci

scatter s1* _t
