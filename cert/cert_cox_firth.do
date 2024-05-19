clear
set obs 3
input stime died x
1 1 1
2 1 1
3 1 0

merlin (stime x, family(cox, firth failure(died)))

assert `"`e(cmdline)'"'           == `"merlin (stime x, family(cox, firth failure(died)))"'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"0"'
assert `"`e(cmplabels1)'"'        == `"x"'
assert `"`e(Nvars_1)'"'           == `"1"'
assert `"`e(timevar1)'"'          == `"stime"'
assert `"`e(failure1)'"'          == `"died"'
assert `"`e(response1)'"'         == `"stime died"'
assert `"`e(family1)'"'           == `"cox"'
assert `"`e(allvars)'"'           == `"died stime x"'
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

assert         e(rank)       == 1
assert         e(N)          == 3
assert         e(k)          == 1
assert         e(k_eq)       == 1
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -1.709334009787867) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

mat T_b=  1.329103947653489
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"_cmp_1_1_1:_cons"'
mat drop C_b T_b

mat T_V=  3.723527300202224
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"_cmp_1_1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"_cmp_1_1_1:_cons"'
mat drop C_V T_V

mat T_gradient=  2.30488544917e-12
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"_cmp_1_1_1:_cons"'
mat drop C_gradient T_gradient

mat T_ml_scale=  .6835394680821688
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1"'
mat drop C_ml_scale T_ml_scale


stset stime, f(died)
stmerlin x, dist(cox) firth 

assert `"`e(cmdline2)'"'          == `"stmerlin x, dist(cox) firth"'
assert `"`e(cmd2)'"'              == `"stmerlin"'
_assert_streq `"`e(cmdline)'"' `"merlin (_t x    , family(cox, failure(_d)      firth)  ) if _st==1, bors  from(__000002)         "'
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"0"'
assert `"`e(constant1)'"'         == `"0"'
assert `"`e(cmplabels1)'"'        == `"x"'
assert `"`e(Nvars_1)'"'           == `"1"'
assert `"`e(timevar1)'"'          == `"_t"'
assert `"`e(failure1)'"'          == `"_d"'
assert `"`e(response1)'"'         == `"_t _d"'
assert `"`e(family1)'"'           == `"cox"'
assert `"`e(allvars)'"'           == `"_d _t x"'
assert `"`e(title)'"'             == `"Survival model"'
assert `"`e(cmd)'"'               == `"merlin"'
assert `"`e(hasopts)'"'           == `"1"'
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

assert         e(rank)       == 1
assert         e(N)          == 3
assert         e(k)          == 1
assert         e(k_eq)       == 1
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -1.709334009787866) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 1
assert         e(Nlevels)    == 1

mat T_b=  1.329103809473645
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"_cmp_1_1_1:_cons"'
mat drop C_b T_b

mat T_V=  3.722467235758117
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"_cmp_1_1_1:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"_cmp_1_1_1:_cons"'
mat drop C_V T_V

mat T_gradient=  3.03544876179e-08
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'
_assert_streq `"`: colfullnames C_gradient'"' `"_cmp_1_1_1:_cons"'
mat drop C_gradient T_gradient

mat T_ml_scale=  .3556074371471415
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1"'
mat drop C_ml_scale T_ml_scale
