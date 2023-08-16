//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

use "/Users/Michael/Documents/stdev/gml/data/tvsfpors",clear

// meologit thk prethk cc tv || school:, trace startvalues(zero)
//meoprobit thk prethk cc tv || school:

// gsem (thk <-  , family(ordinal)), trace

// merlin (thk , family(ologit)), trace //restartvalues(M2 0.1)

sort school thk
merlin (thk prethk cc tv M2[school]@1, family(oprobit)), 
