//=================================================================================================================//
//
// development instructions for merlin
//
//=================================================================================================================//

// adding a new element

error checks in merlin_setup_check_clp()
extract varname (& possibly offsets) for merlin_setup_build_touses()
lots in merlin_setup_build_clp()
setup function in merlin_setup_elements()
needs new function in merlin_functions.mata
add to util_xzb()
add to eret
if possibly time-dependent then add to setup check for NI

// adding a new family

add to ParseDist in merlin_parse.ado
xb pointer in merlin_setup.mata
logl pointer in merlin_setup.mata
Ndistancp in merlin_setup.mata
nbetas for gf12
predict functions
help file
//->if survival 
predict, issurv
predictms will automatically work, add to doc
survsim needs error check removing for the family, add to doc
stmerlin needs error check removing for the family, add to doc
stmixed ...
