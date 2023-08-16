//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

webuse brcancer,clear
replace rectime = rectime + 0.001*runiform()
stset rectime, f(censrec) //scale(365.25)

merlin (_t hormon x1, family(cox, failure(_d)))

predict s1, logchazard ci

stcox hormon x1
// predictnl s2 = predict(basechazard) , ci(s2_lci s2_uci) force
