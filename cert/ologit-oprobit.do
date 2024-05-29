
//ologit 

use ./data/ologit.dta, clear

merlin (apply gpa, family(ologit))

assert `"`e(cmdline)'"'           == `"merlin (apply gpa, family(ologit))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"2"'
assert `"`e(constant1)'"'         == `"0"'
assert `"`e(cmplabels1)'"'        == `"gpa"'
assert `"`e(Nvars_1)'"'           == `"1"'
assert `"`e(response1)'"'         == `"apply"'
assert `"`e(family1)'"'           == `"ordinal"'
assert `"`e(allvars)'"'           == `"apply gpa"'
assert `"`e(title)'"'             == `"Fixed effects regression model"'
assert `"`e(cmd)'"'               == `"merlin"'
assert `"`e(hasopts)'"'           == `"0"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"merlin_p"'
assert `"`e(deriv_useminbound)'"' == `"off"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"merlin_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf0"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 3
assert         e(N)          == 400
assert         e(k)          == 3
assert         e(k_eq)       == 3
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -366.2999806662056) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,3,0)
mat T_b[1,1] =  .7248720108409579
mat T_b[1,2] =  2.374854794663344
mat T_b[1,3] =  4.399911382224572
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"_cmp_1_1_1:_cons dap1_1:_cons dap1_2:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(3,3,0)
mat T_V[1,1] =  .0621516126747953
mat T_V[1,2] =  .1870181106889627
mat T_V[1,3] =  .1903162268740929
mat T_V[2,1] =  .1870181106889627
mat T_V[2,2] =  .5730410479678641
mat T_V[2,3] =  .5788045799176031
mat T_V[3,1] =  .1903162268740929
mat T_V[3,2] =  .5788045799176031
mat T_V[3,3] =  .6107160221211283
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"_cmp_1_1_1:_cons dap1_1:_cons dap1_2:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"_cmp_1_1_1:_cons dap1_1:_cons dap1_2:_cons"'
mat drop C_V T_V

//oprobit

merlin (apply gpa, family(oprobit))

assert `"`e(cmdline)'"'           == `"merlin (apply gpa, family(oprobit))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"2"'
assert `"`e(constant1)'"'         == `"0"'
assert `"`e(cmplabels1)'"'        == `"gpa"'
assert `"`e(Nvars_1)'"'           == `"1"'
assert `"`e(response1)'"'         == `"apply"'
assert `"`e(family1)'"'           == `"ordinal"'
assert `"`e(allvars)'"'           == `"apply gpa"'
assert `"`e(title)'"'             == `"Fixed effects regression model"'
assert `"`e(cmd)'"'               == `"merlin"'
assert `"`e(hasopts)'"'           == `"0"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"merlin_p"'
assert `"`e(deriv_useminbound)'"' == `"off"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"merlin_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf0"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 3
assert         e(N)          == 400
assert         e(k)          == 3
assert         e(k_eq)       == 3
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -365.9264182022463) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,3,0)
mat T_b[1,1] =  .4568953461475451
mat T_b[1,2] =   1.49734204304666
mat T_b[1,3] =  2.672959346315235
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"_cmp_1_1_1:_cons dap1_1:_cons dap1_2:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(3,3,0)
mat T_V[1,1] =  .0225514712689514
mat T_V[1,2] =  .0678317962639689
mat T_V[1,3] =  .0688490049265496
mat T_V[2,1] =  .0678317962639689
mat T_V[2,2] =   .208021844519503
mat T_V[2,3] =  .2090506089016653
mat T_V[3,1] =  .0688490049265496
mat T_V[3,2] =  .2090506089016653
mat T_V[3,3] =  .2176390392091212
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"_cmp_1_1_1:_cons dap1_1:_cons dap1_2:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"_cmp_1_1_1:_cons dap1_1:_cons dap1_2:_cons"'
mat drop C_V T_V

