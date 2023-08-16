use https://www.mjcrowther.co.uk/data/jm_example.dta, clear

merlin	(logb  time M1[id]@1              	, family(gaussian))                 ///
		(logp  time M2[id]@1            , family(gaussian))                ///
		(stime trt M1[id] M2[id]	, family(weibull, failure(died)))   ///
		, covariance(unstructured) 
do ./cert/predictions.do
		   
mat startvals = e(b)
set seed 9741308
merlin 	(logb  time M1[id]@1              , family(gaussian))                 	///
		(logp  time M2[id]@1              , family(gaussian))                 	///
		(stime trt M1[id] M2[id]		  , family(weibull, failure(died)))   	///
		, covariance(unstructured) intmethod(mcarlo) intpoints(300)         	///
		redist(t) df(3) from(startvals) 
do ./cert/predictions.do

merlin  (logb  fp(time, pow(1)) M1[id]@1            , family(gaussian) timevar(time))   ///
        (logp  fp(time, pow(1)) M2[id]@1            , family(gaussian) timevar(time))   ///
		(stime trt EV[logb] EV[logp]  				, family(weibull, failure(died)) timevar(stime))   ///
		, covariance(unstructured) 
do ./cert/predictions.do

merlin 	(logb  fp(time, pow(1)) M1[id]@1              , family(gaussian) timevar(time))   	///
		(logp  fp(time, pow(1)) M2[id]@1              , family(gaussian) timevar(time))   	///
		(stime trt EV[logb] EV[logp] EV[logb]#EV[logp] , 					              	///
                                                family(weibull, failure(died)) timevar(stime))  			///
		, covariance(unstructured) 
do ./cert/predictions.do

replace logp = . in 2
merlin	(logb  time M1[id]@1              	, family(gaussian))                 ///
	(logp  time M2[id]@1              	, family(gaussian))                	///
	(stime trt M1[id] M2[id]			, family(weibull, failure(died)))   ///
	, covariance(unstructured) 
do ./cert/predictions.do
