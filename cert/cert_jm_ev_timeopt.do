set seed 7254
clear
set obs 1000
gen id1 = _n
gen u0 = rnormal(0,1)
gen u1 = rnormal(0,0.3)
gen trt = runiform()>0.5

survsim stime died, hazard(0.1:*1.2:*{t}:^0.2 :* /// baseline hazard
        exp(0.1 :* (0.5 :*trt :+ u0 :+ (0.1:+u1):*{t} :+ 0.2:*{t}:^2) ///
		:+ 0.3 :* (0.5 :*trt :+ u0))) 	///
        covariates(trt -0.5) 	///
        maxt(10)		//	admin. censoring

expand 10
bys id1 : gen time = _n-1
drop if time>stime

gen xb = u0 + (0.1+u1) * time + 0.2 * time^2 + 0.5 * trt
gen y = rnormal(xb,0.5)

bys id1 (time) : replace stime 	= . if _n>1
bys id1 (time) : replace died 	= . if _n>1
	
merlin 	(y trt time fp(time, pow(2)) time#M2[id]@1 M1[id]@1, 	///
		family(gaussian) timevar(time))			///
	(stime trt EV[y, time(0)] EV[y], 			///
		family(w, failure(died)) timevar(stime))	///
		, cov(unstr)
		
assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(intmethod1)'"'        == `"mvaghermite"'
assert `"`e(intpoints1)'"'        == `"7"'
assert `"`e(re_dist1)'"'          == `"normal"'
assert `"`e(re_label1)'"'         == `"sd(M1) sd(M2) corr(M2,M1)"'
assert `"`e(re_ivscale1)'"'       == `"exp exp tanh"'
assert `"`e(re_eqns1)'"'          == `"lns1_1 lns1_2 art1_1_2"'
assert `"`e(Nreparams1)'"'        == `"3"'
assert `"`e(Nres1)'"'             == `"2"'
assert `"`e(latents1)'"'          == `"M1 M2"'
assert `"`e(nap2)'"'              == `"0"'
assert `"`e(ndistap2)'"'          == `"1"'
assert `"`e(constant2)'"'         == `"1"'
assert `"`e(cmplabels2)'"'        == `"trt EV[,time()] EV[]"'
assert `"`e(Nvars_2)'"'           == `"1 1 1"'
assert `"`e(timevar2)'"'          == `"stime"'
assert `"`e(failure2)'"'          == `"died"'
assert `"`e(response2)'"'         == `"stime died"'
assert `"`e(family2)'"'           == `"weibull"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"1"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(cmplabels1)'"'        == `"trt time fp() time#M2[id] M1[id]"'
assert `"`e(Nvars_1)'"'           == `"1 1 1 1 1"'
assert `"`e(timevar1)'"'          == `"time"'
assert `"`e(response1)'"'         == `"y"'
assert `"`e(family1)'"'           == `"gaussian"'
assert `"`e(levelvars)'"'         == `"id"'
assert `"`e(allvars)'"'           == `"died id stime time trt y"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"merlin"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"merlin_p"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"merlin_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf0"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 13
assert         e(N)          == 5240
assert         e(k)          == 15
assert         e(k_eq)       == 15
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -8432.809734956487) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 2
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,15,0)
mat T_b[1,1] =  .5114156772091157
mat T_b[1,2] =   .108580961411121
mat T_b[1,3] =  .1995054706589502
mat T_b[1,4] =                  1
mat T_b[1,5] =                  1
mat T_b[1,6] =  .0642507537519986
mat T_b[1,7] = -.6831125041950871
mat T_b[1,8] = -.6472969262902526
mat T_b[1,9] =  .2844112418961536
mat T_b[1,10] =  .1040293170163463
mat T_b[1,11] =  -2.15991577833493
mat T_b[1,12] =  .1554578039447398
mat T_b[1,13] =  .0449009717273517
mat T_b[1,14] = -1.190829815879932
mat T_b[1,15] =  .0172445082444099
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'


mat drop C_b T_b

qui {
mat T_V = J(15,15,0)
mat T_V[1,1] =  .0049615839359173
mat T_V[1,2] = -.0000133005877004
mat T_V[1,3] =  3.24823892357e-06
mat T_V[1,6] = -.0026728170245046
mat T_V[1,7] = -3.52885174649e-08
mat T_V[1,8] = -.0001274364173335
mat T_V[1,9] = -.0000148737662817
mat T_V[1,10] = -4.55629332415e-08
mat T_V[1,11] =  .0000896819540141
mat T_V[1,12] = -7.34603021249e-06
mat T_V[1,13] = -5.34369331035e-06
mat T_V[1,14] = -2.12331425896e-06
mat T_V[1,15] =   .000030562140878
mat T_V[2,1] = -.0000133005877004
mat T_V[2,2] =  .0002185901723534
mat T_V[2,3] = -.0000131267955907
mat T_V[2,6] = -.0000895500738324
mat T_V[2,7] =  6.07907333081e-07
mat T_V[2,8] = -.0000136672844967
mat T_V[2,9] =  .0000245741404501
mat T_V[2,10] =  3.18065386883e-06
mat T_V[2,11] =  .0000174914192116
mat T_V[2,12] = -1.72072574970e-06
mat T_V[2,13] =  6.76849779888e-06
mat T_V[2,14] = -4.31061157187e-06
mat T_V[2,15] =  .0000455191431386
mat T_V[3,1] =  3.24823892357e-06
mat T_V[3,2] = -.0000131267955907
mat T_V[3,3] =  2.38861671569e-06
mat T_V[3,6] =  8.08966203382e-06
mat T_V[3,7] =  7.07742408844e-08
mat T_V[3,8] =  2.34460895697e-06
mat T_V[3,9] = -5.82026443422e-06
mat T_V[3,10] =  1.76310212014e-07
mat T_V[3,11] =  3.32293424032e-06
mat T_V[3,12] = -3.35956157350e-06
mat T_V[3,13] = -3.08489497311e-06
mat T_V[3,14] =  1.93777609059e-06
mat T_V[3,15] =  4.80047596649e-06
mat T_V[6,1] = -.0026728170245046
mat T_V[6,2] = -.0000895500738324
mat T_V[6,3] =  8.08966203382e-06
mat T_V[6,6] =  .0027246869684819
mat T_V[6,7] = -3.29033923222e-07
mat T_V[6,8] =    .00006572587437
mat T_V[6,9] =  .0000132765518059
mat T_V[6,10] = -1.36194824190e-06
mat T_V[6,11] = -.0001164338430396
mat T_V[6,12] =  .0000127011859026
mat T_V[6,13] =  6.11218457723e-06
mat T_V[6,14] =  1.55863894995e-06
mat T_V[6,15] = -.0000284890730886
mat T_V[7,1] = -3.52885174649e-08
mat T_V[7,2] =  6.07907333081e-07
mat T_V[7,3] =  7.07742408844e-08
mat T_V[7,6] = -3.29033923222e-07
mat T_V[7,7] =  .0001433684804711
mat T_V[7,8] = -5.36232141894e-06
mat T_V[7,9] =  9.06686502948e-06
mat T_V[7,10] =  1.13040441449e-08
mat T_V[7,11] = -2.21643123130e-06
mat T_V[7,12] =  1.19095324053e-06
mat T_V[7,13] = -.0000185037385701
mat T_V[7,14] = -.0000247689155763
mat T_V[7,15] =  .0000349188715619
mat T_V[8,1] = -.0001274364173335
mat T_V[8,2] = -.0000136672844967
mat T_V[8,3] =  2.34460895697e-06
mat T_V[8,6] =    .00006572587437
mat T_V[8,7] = -5.36232141894e-06
mat T_V[8,8] =  .0047824093356833
mat T_V[8,9] = -.0006759592745774
mat T_V[8,10] = -.0000710192949182
mat T_V[8,11] =   -.00155140758058
mat T_V[8,12] = -.0001100589082614
mat T_V[8,13] =  .0000149944406699
mat T_V[8,14] = -1.09077403602e-06
mat T_V[8,15] = -.0000177507682139
mat T_V[9,1] = -.0000148737662817
mat T_V[9,2] =  .0000245741404501
mat T_V[9,3] = -5.82026443422e-06
mat T_V[9,6] =  .0000132765518059
mat T_V[9,7] =  9.06686502948e-06
mat T_V[9,8] = -.0006759592745774
mat T_V[9,9] =  .0012428199180261
mat T_V[9,10] = -.0000211521871774
mat T_V[9,11] = -.0007872940863251
mat T_V[9,12] =  .0003449991284563
mat T_V[9,13] = -.0000278818654587
mat T_V[9,14] = -5.63525516869e-06
mat T_V[9,15] =  .0000337886025187
mat T_V[10,1] = -4.55629332415e-08
mat T_V[10,2] =  3.18065386883e-06
mat T_V[10,3] =  1.76310212014e-07
mat T_V[10,6] = -1.36194824190e-06
mat T_V[10,7] =  1.13040441449e-08
mat T_V[10,8] = -.0000710192949182
mat T_V[10,9] = -.0000211521871774
mat T_V[10,10] =  .0000747693250866
mat T_V[10,11] =  .0002202152816848
mat T_V[10,12] = -.0002682185482921
mat T_V[10,13] = -5.65741695403e-07
mat T_V[10,14] =  2.63832856143e-06
mat T_V[10,15] =  1.45862741649e-06
mat T_V[11,1] =  .0000896819540141
mat T_V[11,2] =  .0000174914192116
mat T_V[11,3] =  3.32293424032e-06
mat T_V[11,6] = -.0001164338430396
mat T_V[11,7] = -2.21643123130e-06
mat T_V[11,8] =   -.00155140758058
mat T_V[11,9] = -.0007872940863251
mat T_V[11,10] =  .0002202152816848
mat T_V[11,11] =   .008004062107746
mat T_V[11,12] = -.0029906701900228
mat T_V[11,13] = -.0000126145950315
mat T_V[11,14] =  .0000139506471219
mat T_V[11,15] =   .000016891597225
mat T_V[12,1] = -7.34603021249e-06
mat T_V[12,2] = -1.72072574970e-06
mat T_V[12,3] = -3.35956157350e-06
mat T_V[12,6] =  .0000127011859026
mat T_V[12,7] =  1.19095324053e-06
mat T_V[12,8] = -.0001100589082614
mat T_V[12,9] =  .0003449991284563
mat T_V[12,10] = -.0002682185482921
mat T_V[12,11] = -.0029906701900228
mat T_V[12,12] =  .0019121423695886
mat T_V[12,13] =  7.70735191729e-06
mat T_V[12,14] = -.0000120199214096
mat T_V[12,15] = -.0000122001502848
mat T_V[13,1] = -5.34369331035e-06
mat T_V[13,2] =  6.76849779888e-06
mat T_V[13,3] = -3.08489497311e-06
mat T_V[13,6] =  6.11218457723e-06
mat T_V[13,7] = -.0000185037385701
mat T_V[13,8] =  .0000149944406699
mat T_V[13,9] = -.0000278818654587
mat T_V[13,10] = -5.65741695403e-07
mat T_V[13,11] = -.0000126145950315
mat T_V[13,12] =  7.70735191729e-06
mat T_V[13,13] =   .000650595158777
mat T_V[13,14] =  5.82761134227e-06
mat T_V[13,15] = -.0001679820998379
mat T_V[14,1] = -2.12331425896e-06
mat T_V[14,2] = -4.31061157187e-06
mat T_V[14,3] =  1.93777609059e-06
mat T_V[14,6] =  1.55863894995e-06
mat T_V[14,7] = -.0000247689155763
mat T_V[14,8] = -1.09077403602e-06
mat T_V[14,9] = -5.63525516869e-06
mat T_V[14,10] =  2.63832856143e-06
mat T_V[14,11] =  .0000139506471219
mat T_V[14,12] = -.0000120199214096
mat T_V[14,13] =  5.82761134227e-06
mat T_V[14,14] =  .0009170104171091
mat T_V[14,15] = -.0001641122099889
mat T_V[15,1] =   .000030562140878
mat T_V[15,2] =  .0000455191431386
mat T_V[15,3] =  4.80047596649e-06
mat T_V[15,6] = -.0000284890730886
mat T_V[15,7] =  .0000349188715619
mat T_V[15,8] = -.0000177507682139
mat T_V[15,9] =  .0000337886025187
mat T_V[15,10] =  1.45862741649e-06
mat T_V[15,11] =   .000016891597225
mat T_V[15,12] = -.0000122001502848
mat T_V[15,13] = -.0001679820998379
mat T_V[15,14] = -.0001641122099889
mat T_V[15,15] =  .0018084195840078
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8

mat drop C_V T_V

qui {
mat T_gradient = J(1,15,0)
mat T_gradient[1,1] = -1.43050596075e-07
mat T_gradient[1,2] = -2.18772634291e-07
mat T_gradient[1,3] = -1.25245078688e-07
mat T_gradient[1,6] = -2.21441356179e-07
mat T_gradient[1,7] =  6.65646573133e-08
mat T_gradient[1,8] = -3.58825198569e-09
mat T_gradient[1,9] =  2.47694581637e-09
mat T_gradient[1,10] = -3.50254221676e-08
mat T_gradient[1,11] = -5.68156150252e-09
mat T_gradient[1,12] = -1.34509064928e-08
mat T_gradient[1,13] =  3.52859884513e-07
mat T_gradient[1,14] =  4.92458120109e-08
mat T_gradient[1,15] =  1.55451244089e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'

mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,13,0)
mat T_ml_scale[1,1] =  1.258212076276074
mat T_ml_scale[1,2] =  1.893638696603069
mat T_ml_scale[1,3] =   .091126243616412
mat T_ml_scale[1,4] =  15.40905702224557
mat T_ml_scale[1,5] =  .4270062568902656
mat T_ml_scale[1,6] =  1.179495995180331
mat T_ml_scale[1,7] =  3.822150112028886
mat T_ml_scale[1,8] =  .7851577992537245
mat T_ml_scale[1,9] =  .2027947304742919
mat T_ml_scale[1,10] =  1.140873927277686
mat T_ml_scale[1,11] =  15.44478765272994
mat T_ml_scale[1,12] =  .3617928167638673
mat T_ml_scale[1,13] =  54.87817693695079
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13"'
mat drop C_ml_scale T_ml_scale

do ./cert/predictions.do

merlin 	(y trt time fp(time, pow(2)) time#M2[id]@1 M1[id]@1, 	///
		family(gaussian) timevar(time))			///
	(stime trt EV[y, time(0)], 			///
		family(w, failure(died)) timevar(stime))	///
		, cov(unstr)


assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(intmethod1)'"'        == `"mvaghermite"'
assert `"`e(intpoints1)'"'        == `"7"'
assert `"`e(re_dist1)'"'          == `"normal"'
assert `"`e(re_label1)'"'         == `"sd(M1) sd(M2) corr(M2,M1)"'
assert `"`e(re_ivscale1)'"'       == `"exp exp tanh"'
assert `"`e(re_eqns1)'"'          == `"lns1_1 lns1_2 art1_1_2"'
assert `"`e(Nreparams1)'"'        == `"3"'
assert `"`e(Nres1)'"'             == `"2"'
assert `"`e(latents1)'"'          == `"M1 M2"'
assert `"`e(nap2)'"'              == `"0"'
assert `"`e(ndistap2)'"'          == `"1"'
assert `"`e(constant2)'"'         == `"1"'
assert `"`e(cmplabels2)'"'        == `"trt EV[,time()]"'
assert `"`e(Nvars_2)'"'           == `"1 1"'
assert `"`e(timevar2)'"'          == `"stime"'
assert `"`e(failure2)'"'          == `"died"'
assert `"`e(response2)'"'         == `"stime died"'
assert `"`e(family2)'"'           == `"weibull"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"1"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(cmplabels1)'"'        == `"trt time fp() time#M2[id] M1[id]"'
assert `"`e(Nvars_1)'"'           == `"1 1 1 1 1"'
assert `"`e(timevar1)'"'          == `"time"'
assert `"`e(response1)'"'         == `"y"'
assert `"`e(family1)'"'           == `"gaussian"'
assert `"`e(levelvars)'"'         == `"id"'
assert `"`e(allvars)'"'           == `"died id stime time trt y"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"merlin"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"merlin_p"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"merlin_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf0"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 12
assert         e(N)          == 5240
assert         e(k)          == 14
assert         e(k_eq)       == 14
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -8494.96747302738 ) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 2
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,14,0)
mat T_b[1,1] =   .511120740050131
mat T_b[1,2] =  .1071137635683485
mat T_b[1,3] =  .1987233625280199
mat T_b[1,4] =                  1
mat T_b[1,5] =                  1
mat T_b[1,6] =  .0648787047363927
mat T_b[1,7] = -.6826454716659203
mat T_b[1,8] =  -.572562083505552
mat T_b[1,9] =  .3357212173732097
mat T_b[1,10] = -2.502720599851075
mat T_b[1,11] =  .4834238255999958
mat T_b[1,12] =   .045714107604905
mat T_b[1,13] = -1.196302265798748
mat T_b[1,14] =  .0161976726083015
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"_cmp_1_1_1:_cons _cmp_1_2_1:_cons _cmp_1_3_1:_cons _cmp_1_4_1:_cons _cmp_1_5_1:_cons cons1:_cons dap1_1:_cons _cmp_2_1_1:_cons _cmp_2_2_1:_cons cons2:_cons dap2_1:_cons lns1_1:_cons lns1_2:_cons art1_1_2:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(14,14,0)
mat T_V[1,1] =  .0049678097227934
mat T_V[1,2] = -.0000135969132577
mat T_V[1,3] =  3.20234732218e-06
mat T_V[1,6] = -.0026760987341562
mat T_V[1,7] = -9.84839214618e-08
mat T_V[1,8] = -.0001467751280891
mat T_V[1,9] = -.0000140833196023
mat T_V[1,10] =   .000093798262829
mat T_V[1,11] = -4.35531441395e-06
mat T_V[1,12] = -4.42120576354e-06
mat T_V[1,13] = -2.32789583780e-06
mat T_V[1,14] =  .0000205919708522
mat T_V[2,1] = -.0000135969132577
mat T_V[2,2] =  .0002173278627755
mat T_V[2,3] = -.0000131345879869
mat T_V[2,6] = -.0000897939291995
mat T_V[2,7] =  8.47901457058e-07
mat T_V[2,8] = -.0000105846896318
mat T_V[2,9] =  .0000236898521857
mat T_V[2,10] =  .0000181779326472
mat T_V[2,11] =  6.61821660407e-06
mat T_V[2,12] =  6.68864648478e-06
mat T_V[2,13] = -8.19107231537e-06
mat T_V[2,14] =  .0000457011845493
mat T_V[3,1] =  3.20234732218e-06
mat T_V[3,2] = -.0000131345879869
mat T_V[3,3] =  2.38204966204e-06
mat T_V[3,6] =  8.13838455264e-06
mat T_V[3,7] = -9.01813236083e-08
mat T_V[3,8] =  1.95114852836e-06
mat T_V[3,9] = -4.71061803126e-06
mat T_V[3,10] =  1.00706969764e-06
mat T_V[3,11] = -1.71439792056e-06
mat T_V[3,12] = -3.07844732878e-06
mat T_V[3,13] =  1.67268867437e-06
mat T_V[3,14] =  4.94933243949e-06
mat T_V[6,1] = -.0026760987341562
mat T_V[6,2] = -.0000897939291995
mat T_V[6,3] =  8.13838455264e-06
mat T_V[6,6] =  .0027283773639294
mat T_V[6,7] = -3.72582950066e-07
mat T_V[6,8] =  .0000749487554227
mat T_V[6,9] =  .0000130665337119
mat T_V[6,10] = -.0001144094736133
mat T_V[6,11] =  2.91793340935e-06
mat T_V[6,12] =  5.51142914814e-06
mat T_V[6,13] =  3.01972048619e-06
mat T_V[6,14] =  -.000023482312296
mat T_V[7,1] = -9.84839214618e-08
mat T_V[7,2] =  8.47901457058e-07
mat T_V[7,3] = -9.01813236083e-08
mat T_V[7,6] = -3.72582950066e-07
mat T_V[7,7] =  .0001435956868135
mat T_V[7,8] = -6.45307987715e-06
mat T_V[7,9] =  .0000108454586581
mat T_V[7,10] = -3.45134816282e-06
mat T_V[7,11] =  1.48453484305e-06
mat T_V[7,12] = -.0000184852394196
mat T_V[7,13] = -.0000264774053787
mat T_V[7,14] =  .0000360985416701
mat T_V[8,1] = -.0001467751280891
mat T_V[8,2] = -.0000105846896318
mat T_V[8,3] =  1.95114852836e-06
mat T_V[8,6] =  .0000749487554227
mat T_V[8,7] = -6.45307987715e-06
mat T_V[8,8] =  .0047148130072164
mat T_V[8,9] =  -.000688763230362
mat T_V[8,10] = -.0014085385296659
mat T_V[8,11] = -.0002755300910122
mat T_V[8,12] =  .0000185710239294
mat T_V[8,13] =  9.86431114601e-06
mat T_V[8,14] = -.0000378340635647
mat T_V[9,1] = -.0000140833196023
mat T_V[9,2] =  .0000236898521857
mat T_V[9,3] = -4.71061803126e-06
mat T_V[9,6] =  .0000130665337119
mat T_V[9,7] =  .0000108454586581
mat T_V[9,8] =  -.000688763230362
mat T_V[9,9] =  .0012032284569466
mat T_V[9,10] =  -.000652086241621
mat T_V[9,11] =  .0001967542816783
mat T_V[9,12] = -.0000358731138076
mat T_V[9,13] = -.0000184939598075
mat T_V[9,14] =  .0000712451662442
mat T_V[10,1] =   .000093798262829
mat T_V[10,2] =  .0000181779326472
mat T_V[10,3] =  1.00706969764e-06
mat T_V[10,6] = -.0001144094736133
mat T_V[10,7] = -3.45134816282e-06
mat T_V[10,8] = -.0014085385296659
mat T_V[10,9] =  -.000652086241621
mat T_V[10,10] =  .0084413481927368
mat T_V[10,11] = -.0021532862640791
mat T_V[10,12] = -4.40539108448e-06
mat T_V[10,13] =  .0000123593942539
mat T_V[10,14] = -2.98253514670e-06
mat T_V[11,1] = -4.35531441395e-06
mat T_V[11,2] =  6.61821660407e-06
mat T_V[11,3] = -1.71439792056e-06
mat T_V[11,6] =  2.91793340935e-06
mat T_V[11,7] =  1.48453484305e-06
mat T_V[11,8] = -.0002755300910122
mat T_V[11,9] =  .0001967542816783
mat T_V[11,10] = -.0021532862640791
mat T_V[11,11] =  .0007505240377315
mat T_V[11,12] =  2.68019483602e-06
mat T_V[11,13] = -5.01617616897e-06
mat T_V[11,14] =  1.73276169756e-06
mat T_V[12,1] = -4.42120576354e-06
mat T_V[12,2] =  6.68864648478e-06
mat T_V[12,3] = -3.07844732878e-06
mat T_V[12,6] =  5.51142914814e-06
mat T_V[12,7] = -.0000184852394196
mat T_V[12,8] =  .0000185710239294
mat T_V[12,9] = -.0000358731138076
mat T_V[12,10] = -4.40539108448e-06
mat T_V[12,11] =  2.68019483602e-06
mat T_V[12,12] =  .0006504832268623
mat T_V[12,13] =  7.66429138397e-06
mat T_V[12,14] = -.0001706208395836
mat T_V[13,1] = -2.32789583780e-06
mat T_V[13,2] = -8.19107231537e-06
mat T_V[13,3] =  1.67268867437e-06
mat T_V[13,6] =  3.01972048619e-06
mat T_V[13,7] = -.0000264774053787
mat T_V[13,8] =  9.86431114601e-06
mat T_V[13,9] = -.0000184939598075
mat T_V[13,10] =  .0000123593942539
mat T_V[13,11] = -5.01617616897e-06
mat T_V[13,12] =  7.66429138397e-06
mat T_V[13,13] =  .0009195773569102
mat T_V[13,14] = -.0001737214663782
mat T_V[14,1] =  .0000205919708522
mat T_V[14,2] =  .0000457011845493
mat T_V[14,3] =  4.94933243949e-06
mat T_V[14,6] =  -.000023482312296
mat T_V[14,7] =  .0000360985416701
mat T_V[14,8] = -.0000378340635647
mat T_V[14,9] =  .0000712451662442
mat T_V[14,10] = -2.98253514670e-06
mat T_V[14,11] =  1.73276169756e-06
mat T_V[14,12] = -.0001706208395836
mat T_V[14,13] = -.0001737214663782
mat T_V[14,14] =  .0018166979923904
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"_cmp_1_1_1:_cons _cmp_1_2_1:_cons _cmp_1_3_1:_cons _cmp_1_4_1:_cons _cmp_1_5_1:_cons cons1:_cons dap1_1:_cons _cmp_2_1_1:_cons _cmp_2_2_1:_cons cons2:_cons dap2_1:_cons lns1_1:_cons lns1_2:_cons art1_1_2:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"_cmp_1_1_1:_cons _cmp_1_2_1:_cons _cmp_1_3_1:_cons _cmp_1_4_1:_cons _cmp_1_5_1:_cons cons1:_cons dap1_1:_cons _cmp_2_1_1:_cons _cmp_2_2_1:_cons cons2:_cons dap2_1:_cons lns1_1:_cons lns1_2:_cons art1_1_2:_cons"'
mat drop C_V T_V

qui {
mat T_gradient = J(1,14,0)
mat T_gradient[1,1] = -1.72059398097e-07
mat T_gradient[1,2] = -1.26058667524e-06
mat T_gradient[1,3] = -6.18858304966e-06
mat T_gradient[1,6] = -2.98800812083e-07
mat T_gradient[1,7] =  1.15703596655e-07
mat T_gradient[1,8] = -7.96377785115e-10
mat T_gradient[1,9] = -6.89108814933e-10
mat T_gradient[1,10] = -9.24535280365e-10
mat T_gradient[1,11] = -6.50917848354e-09
mat T_gradient[1,12] =  3.88081913493e-07
mat T_gradient[1,13] =  5.57636743205e-07
mat T_gradient[1,14] =  6.95455728910e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'

mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,12,0)
mat T_ml_scale[1,1] =  1.578455007202345
mat T_ml_scale[1,2] =  1.877325760188616
mat T_ml_scale[1,3] =  .1534647611563707
mat T_ml_scale[1,4] =  8.229089752610335
mat T_ml_scale[1,5] =  .2317341415675373
mat T_ml_scale[1,6] =  1.862799765586885
mat T_ml_scale[1,7] =  1.643052461458864
mat T_ml_scale[1,8] =  .2430738524519503
mat T_ml_scale[1,9] =  .4228650067359527
mat T_ml_scale[1,10] =  12.91835391197532
mat T_ml_scale[1,11] =  .3402894922670585
mat T_ml_scale[1,12] =   48.7597637874115
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12"'
mat drop C_ml_scale T_ml_scale

do ./cert/predictions.do		
		
	
		
merlin 	(y trt time fp(time, pow(2)) time#M2[id]@1 M1[id]@1, 	///
		family(gaussian) timevar(time))			///
	(stime trt EV[y, time(0)], 			///
		family(w, failure(died)))	///
		, cov(unstr)


assert `"`e(chintpoints)'"'       == `"30"'
assert `"`e(intmethod1)'"'        == `"mvaghermite"'
assert `"`e(intpoints1)'"'        == `"7"'
assert `"`e(re_dist1)'"'          == `"normal"'
assert `"`e(re_label1)'"'         == `"sd(M1) sd(M2) corr(M2,M1)"'
assert `"`e(re_ivscale1)'"'       == `"exp exp tanh"'
assert `"`e(re_eqns1)'"'          == `"lns1_1 lns1_2 art1_1_2"'
assert `"`e(Nreparams1)'"'        == `"3"'
assert `"`e(Nres1)'"'             == `"2"'
assert `"`e(latents1)'"'          == `"M1 M2"'
assert `"`e(nap2)'"'              == `"0"'
assert `"`e(ndistap2)'"'          == `"1"'
assert `"`e(constant2)'"'         == `"1"'
assert `"`e(cmplabels2)'"'        == `"trt EV[,time()]"'
assert `"`e(Nvars_2)'"'           == `"1 1"'
assert `"`e(timevar2)'"'          == `"stime"'
assert `"`e(failure2)'"'          == `"died"'
assert `"`e(response2)'"'         == `"stime died"'
assert `"`e(family2)'"'           == `"weibull"'
assert `"`e(nap1)'"'              == `"0"'
assert `"`e(ndistap1)'"'          == `"1"'
assert `"`e(constant1)'"'         == `"1"'
assert `"`e(cmplabels1)'"'        == `"trt time fp() time#M2[id] M1[id]"'
assert `"`e(Nvars_1)'"'           == `"1 1 1 1 1"'
assert `"`e(timevar1)'"'          == `"time"'
assert `"`e(response1)'"'         == `"y"'
assert `"`e(family1)'"'           == `"gaussian"'
assert `"`e(levelvars)'"'         == `"id"'
assert `"`e(allvars)'"'           == `"died id stime time trt y"'
assert `"`e(title)'"'             == `"Mixed effects regression model"'
assert `"`e(cmd)'"'               == `"merlin"'
assert `"`e(hasopts)'"'           == `"1"'
assert `"`e(from)'"'              == `"1"'
assert `"`e(predict)'"'           == `"merlin_p"'
assert `"`e(opt)'"'               == `"moptimize"'
assert `"`e(vce)'"'               == `"oim"'
assert `"`e(user)'"'              == `"merlin_gf()"'
assert `"`e(crittype)'"'          == `"log likelihood"'
assert `"`e(ml_method)'"'         == `"gf0"'
assert `"`e(singularHmethod)'"'   == `"m-marquardt"'
assert `"`e(technique)'"'         == `"nr"'
assert `"`e(which)'"'             == `"max"'
assert `"`e(properties)'"'        == `"b V"'

assert         e(rank)       == 12
assert         e(N)          == 5240
assert         e(k)          == 14
assert         e(k_eq)       == 14
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -8494.965406368865) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 2
assert         e(Nlevels)    == 2

qui {
mat T_b = J(1,14,0)
mat T_b[1,1] =   .511120859429521
mat T_b[1,2] =  .1071135829300313
mat T_b[1,3] =  .1987234095222355
mat T_b[1,4] =                  1
mat T_b[1,5] =                  1
mat T_b[1,6] =  .0648786245268762
mat T_b[1,7] = -.6826454910436669
mat T_b[1,8] =  -.572554551798982
mat T_b[1,9] =  .3357158382739179
mat T_b[1,10] = -2.502659578224493
mat T_b[1,11] =  .4834033121537687
mat T_b[1,12] =  .0457140129975802
mat T_b[1,13] = -1.196302134309318
mat T_b[1,14] =  .0161976519962929
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
_assert_streq `"`: rowfullnames C_b'"' `"y1"'
_assert_streq `"`: colfullnames C_b'"' `"_cmp_1_1_1:_cons _cmp_1_2_1:_cons _cmp_1_3_1:_cons _cmp_1_4_1:_cons _cmp_1_5_1:_cons cons1:_cons dap1_1:_cons _cmp_2_1_1:_cons _cmp_2_2_1:_cons cons2:_cons dap2_1:_cons lns1_1:_cons lns1_2:_cons art1_1_2:_cons"'
mat drop C_b T_b

qui {
mat T_V = J(14,14,0)
mat T_V[1,1] =  .0049678088173964
mat T_V[1,2] = -.0000135969065016
mat T_V[1,3] =  3.20234666164e-06
mat T_V[1,6] = -.0026760982291437
mat T_V[1,7] = -9.84662101593e-08
mat T_V[1,8] = -.0001467726750005
mat T_V[1,9] = -.0000140833261389
mat T_V[1,10] =   .000093799609804
mat T_V[1,11] = -4.35623208996e-06
mat T_V[1,12] = -4.42119627068e-06
mat T_V[1,13] = -2.32789396727e-06
mat T_V[1,14] =  .0000205920180481
mat T_V[2,1] = -.0000135969065016
mat T_V[2,2] =  .0002173278824854
mat T_V[2,3] = -.0000131345893473
mat T_V[2,6] = -.0000897939327689
mat T_V[2,7] =  8.47839487509e-07
mat T_V[2,8] = -.0000105850465873
mat T_V[2,9] =  .0000236898822879
mat T_V[2,10] =  .0000181735943255
mat T_V[2,11] =  6.61957578799e-06
mat T_V[2,12] =  6.68865987636e-06
mat T_V[2,13] = -8.19099703537e-06
mat T_V[2,14] =  .0000457010063735
mat T_V[3,1] =  3.20234666164e-06
mat T_V[3,2] = -.0000131345893473
mat T_V[3,3] =  2.38204987169e-06
mat T_V[3,6] =  8.13838486952e-06
mat T_V[3,7] = -9.01721019102e-08
mat T_V[3,8] =  1.95124274223e-06
mat T_V[3,9] = -4.71062889209e-06
mat T_V[3,10] =  1.00809248241e-06
mat T_V[3,11] = -1.71474659345e-06
mat T_V[3,12] = -3.07844554046e-06
mat T_V[3,13] =  1.67270626199e-06
mat T_V[3,14] =  4.94931227687e-06
mat T_V[6,1] = -.0026760982291437
mat T_V[6,2] = -.0000897939327689
mat T_V[6,3] =  8.13838486952e-06
mat T_V[6,6] =  .0027283768636239
mat T_V[6,7] = -3.72596378761e-07
mat T_V[6,8] =  .0000749474361289
mat T_V[6,9] =  .0000130665285374
mat T_V[6,10] =  -.000114409633673
mat T_V[6,11] =  2.91858413803e-06
mat T_V[6,12] =  5.51143231751e-06
mat T_V[6,13] =  3.01968588263e-06
mat T_V[6,14] = -.0000234823055246
mat T_V[7,1] = -9.84662101593e-08
mat T_V[7,2] =  8.47839487509e-07
mat T_V[7,3] = -9.01721019102e-08
mat T_V[7,6] = -3.72596378761e-07
mat T_V[7,7] =   .000143595670095
mat T_V[7,8] = -6.45288941295e-06
mat T_V[7,9] =  .0000108450887424
mat T_V[7,10] = -3.45206354793e-06
mat T_V[7,11] =  1.48476786265e-06
mat T_V[7,12] =  -.000018485243703
mat T_V[7,13] = -.0000264773724177
mat T_V[7,14] =  .0000360984831132
mat T_V[8,1] = -.0001467726750005
mat T_V[8,2] = -.0000105850465873
mat T_V[8,3] =  1.95124274223e-06
mat T_V[8,6] =  .0000749474361289
mat T_V[8,7] = -6.45288941295e-06
mat T_V[8,8] =  .0047148283257708
mat T_V[8,9] = -.0006887737847991
mat T_V[8,10] = -.0014083733493745
mat T_V[8,11] = -.0002755905901253
mat T_V[8,12] =  .0000185704966079
mat T_V[8,13] =  9.86456407227e-06
mat T_V[8,14] = -.0000378336662439
mat T_V[9,1] = -.0000140833261389
mat T_V[9,2] =  .0000236898822879
mat T_V[9,3] = -4.71062889209e-06
mat T_V[9,6] =  .0000130665285374
mat T_V[9,7] =  .0000108450887424
mat T_V[9,8] = -.0006887737847991
mat T_V[9,9] =  .0012032322954856
mat T_V[9,10] = -.0006522001557568
mat T_V[9,11] =  .0001967962654472
mat T_V[9,12] = -.0000358723939263
mat T_V[9,13] = -.0000184939584628
mat T_V[9,14] =  .0000712443775129
mat T_V[10,1] =   .000093799609804
mat T_V[10,2] =  .0000181735943255
mat T_V[10,3] =  1.00809248241e-06
mat T_V[10,6] =  -.000114409633673
mat T_V[10,7] = -3.45206354793e-06
mat T_V[10,8] = -.0014083733493745
mat T_V[10,9] = -.0006522001557568
mat T_V[10,10] =  .0084426211288257
mat T_V[10,11] =  -.002153759960756
mat T_V[10,12] = -4.40695740479e-06
mat T_V[10,13] =  .0000123623128704
mat T_V[10,14] = -2.98375171108e-06
mat T_V[11,1] = -4.35623208996e-06
mat T_V[11,2] =  6.61957578799e-06
mat T_V[11,3] = -1.71474659345e-06
mat T_V[11,6] =  2.91858413803e-06
mat T_V[11,7] =  1.48476786265e-06
mat T_V[11,8] = -.0002755905901253
mat T_V[11,9] =  .0001967962654472
mat T_V[11,10] =  -.002153759960756
mat T_V[11,11] =  .0007506993611602
mat T_V[11,12] =  2.68075834525e-06
mat T_V[11,13] = -5.01722685858e-06
mat T_V[11,14] =  1.73316809704e-06
mat T_V[12,1] = -4.42119627068e-06
mat T_V[12,2] =  6.68865987636e-06
mat T_V[12,3] = -3.07844554046e-06
mat T_V[12,6] =  5.51143231751e-06
mat T_V[12,7] =  -.000018485243703
mat T_V[12,8] =  .0000185704966079
mat T_V[12,9] = -.0000358723939263
mat T_V[12,10] = -4.40695740479e-06
mat T_V[12,11] =  2.68075834525e-06
mat T_V[12,12] =  .0006504832155188
mat T_V[12,13] =  7.66423615670e-06
mat T_V[12,14] = -.0001706207300378
mat T_V[13,1] = -2.32789396727e-06
mat T_V[13,2] = -8.19099703537e-06
mat T_V[13,3] =  1.67270626199e-06
mat T_V[13,6] =  3.01968588263e-06
mat T_V[13,7] = -.0000264773724177
mat T_V[13,8] =  9.86456407227e-06
mat T_V[13,9] = -.0000184939584628
mat T_V[13,10] =  .0000123623128704
mat T_V[13,11] = -5.01722685858e-06
mat T_V[13,12] =  7.66423615670e-06
mat T_V[13,13] =  .0009195773278022
mat T_V[13,14] = -.0001737212720606
mat T_V[14,1] =  .0000205920180481
mat T_V[14,2] =  .0000457010063735
mat T_V[14,3] =  4.94931227687e-06
mat T_V[14,6] = -.0000234823055246
mat T_V[14,7] =  .0000360984831132
mat T_V[14,8] = -.0000378336662439
mat T_V[14,9] =  .0000712443775129
mat T_V[14,10] = -2.98375171108e-06
mat T_V[14,11] =  1.73316809704e-06
mat T_V[14,12] = -.0001706207300378
mat T_V[14,13] = -.0001737212720606
mat T_V[14,14] =  .0018166970330895
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
_assert_streq `"`: rowfullnames C_V'"' `"_cmp_1_1_1:_cons _cmp_1_2_1:_cons _cmp_1_3_1:_cons _cmp_1_4_1:_cons _cmp_1_5_1:_cons cons1:_cons dap1_1:_cons _cmp_2_1_1:_cons _cmp_2_2_1:_cons cons2:_cons dap2_1:_cons lns1_1:_cons lns1_2:_cons art1_1_2:_cons"'
_assert_streq `"`: colfullnames C_V'"' `"_cmp_1_1_1:_cons _cmp_1_2_1:_cons _cmp_1_3_1:_cons _cmp_1_4_1:_cons _cmp_1_5_1:_cons cons1:_cons dap1_1:_cons _cmp_2_1_1:_cons _cmp_2_2_1:_cons cons2:_cons dap2_1:_cons lns1_1:_cons lns1_2:_cons art1_1_2:_cons"'
mat drop C_V T_V

qui {
mat T_gradient = J(1,14,0)
mat T_gradient[1,1] = -1.53429502287e-07
mat T_gradient[1,2] = -1.16090597228e-06
mat T_gradient[1,3] = -5.75550353488e-06
mat T_gradient[1,6] = -2.65974799908e-07
mat T_gradient[1,7] =  1.01799998384e-07
mat T_gradient[1,8] = -3.78964530329e-10
mat T_gradient[1,9] = -6.40658107294e-10
mat T_gradient[1,10] = -6.77756678905e-11
mat T_gradient[1,11] = -3.93046109201e-09
mat T_gradient[1,12] =  3.37206414830e-07
mat T_gradient[1,13] =  4.94162709461e-07
mat T_gradient[1,14] =  6.39131792852e-07
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
_assert_streq `"`: rowfullnames C_gradient'"' `"r1"'

mat drop C_gradient T_gradient

qui {
mat T_ml_scale = J(1,12,0)
mat T_ml_scale[1,1] =  1.595129290694157
mat T_ml_scale[1,2] =  1.938786061478532
mat T_ml_scale[1,3] =  .1543220402065401
mat T_ml_scale[1,4] =  8.251690122227629
mat T_ml_scale[1,5] =  .4381732882676909
mat T_ml_scale[1,6] =  1.866399970431286
mat T_ml_scale[1,7] =  1.750875888095688
mat T_ml_scale[1,8] =  .2446997991905092
mat T_ml_scale[1,9] =  .4273124678177102
mat T_ml_scale[1,10] =  11.27314178929863
mat T_ml_scale[1,11] =  .3395728970663057
mat T_ml_scale[1,12] =  52.08230198519526
}
matrix C_ml_scale = e(ml_scale)
assert mreldif( C_ml_scale , T_ml_scale ) < 1E-8
_assert_streq `"`: rowfullnames C_ml_scale'"' `"r1"'
_assert_streq `"`: colfullnames C_ml_scale'"' `"c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12"'
mat drop C_ml_scale T_ml_scale


do ./cert/predictions.do
