*! version 0.0.0  ??????2017

/*
//predictions
-> mu										- expected value
-> eta 										- expected value of linear predictor	
-> survival									- survival function, with ltruncation
-> hazard									- hazard function
-> chazard									- cumulative hazard function
-> rmst										- restricted mean survival time
-> marginal									- marginal over random effects
// -> std										- standardised over all covariates
*/

program merlin_p, sortpreserve
        version 14.1
        local vv : display "version " string(_caller()) ":"

        tempname tname
        capture noisily `vv' Predict `tname' `0'
        local rc = c(rc)
        capture drop `tname'*
        capture mata: rmexternal("`tname'")
        exit `rc'
end

program Predict
        version 14.1
        gettoken GML 0 : 0
        syntax  newvarname      ///
                [if] [in] [,                            ///
														///
				OUTcome(string)							///
														/// statistics
                mu                                      ///
                eta                                     ///
                SURVival                                ///
				/*LTruncation(string)*/						///
				Hazard									///
				CHazard									///
				RMST									///
				USER(string)							///
														///
														///
                FIXEDonly                               ///
                /*EBMEANs*/                                 ///
                MARGinal                                ///
														///
				CI										///
				TIMEvar(varname)						///
				AT(string)								///
														///
                INTPoints(numlist int max=1 >0) 		///
														///
				DEBUG									///	NOTDOC
        ]

		local anything `varlist'
		
        if "`intpoints'" == "1" {
                di as err "invalid intpoints() option;"
                di as err "intpoints(1) is not allowed by predict"
                exit 198
        }

		if "`outcome'"=="" {
			local outcome = 1
		}
		
		if "`debug'"!="" {
			local noisily noisily
		}
		
        // parse statistics

        local STAT      `mu'            ///
                        `eta'           ///
                        `survival'      ///
                        `hazard'      	///
						`chazard'      	///
                        `rmst'		    
        opts_exclusive "`STAT'"

        if "`STAT'" == "" {
                di as txt "(option {bf:mu} assumed)"
                local STAT mu
        }

		if "`STAT'"=="marginal" & "`e(levelvars)'"=="" {
			di as error "No random effects to marginalise over"
			exit 1986
		}
		
		if "`ltruncation'"!="" & ("`STAT'"=="survival" | "`timevar'"=="") {
			di as error "ltruncation() only allowed with survival and timevar()"
			exit 198
		}
		
        // parse options
        
        opts_exclusive "`fixedonly' `ebmeans' `marginal'"
		local xbtype "`fixedonly'`ebmeans'`marginal'"
		
		// postestimation sample

        tempname touse
        mark `touse' `if' `in'

		if "`timevar'"!="" {
			markout `touse' `timevar'
			qui count if `touse'==1
			local ptvar ptvar(`timevar')					//merlin_build_touses() updated on this
		}
		
		//integration method
		
		if "`e(levelvars)'"!="" {
			local Nrelevels = e(Nlevels) - 1
			forval k=1/`Nrelevels' {
				local ims `ims' `e(intmethod`k')'
				if "`intpoints'"=="" {
					local ips `ips' `e(intpoints`k')'
				}
				else {
					local ips `ips' `intpoints'
				}
			}
			mata: st_local("ims",subinstr(st_local("ims"),"mvaghermite","ghermite"))
			local intmethods intmethod(`ims')
			local intpoints intpoints(`imp')
		}
		
		
		//====================================================================================================================//
		
		//Preserve data for out of sample prediction etc.
		tempfile newvars 
		preserve	
		
		//Out of sample predictions using at()
		if "`at'" != "" {
			tokenize `at'
			while "`1'"!="" {
				unab 1: `1'
				cap confirm var `2'
				if _rc {
					cap confirm num `2'
					if _rc {
						di in red "invalid at(... `1' `2' ...)"
						exit 198
					}
				}
				qui replace `1' = `2' if `touse'
				mac shift 2
			}
		}	
		
		if "`ci'"=="" {
		
			//get coefficients and refill struct
			tempname best
			mat `best' = e(b)
			
			//remove any options
			local cmd `e(cmdline)'
			gettoken merlin cmd : cmd
			gettoken cmd rhs : cmd, parse(",") bind
			
			//recall merlin
			tempname tousem
			quietly `noisily' merlin_parse `GML'	, touse(`tousem') : `cmd'		///
												, 								///
													predict 					///
													nogen 						///
													from(`best') 				///
													`intmethods' 				///
													`intpoints' 				///
													`ptvar'
			
			if "`fixedonly'" != "" | "`marginal'"!="" {

				mata: merlin_predict("`GML'","`anything'","`touse'","`STAT'","`xbtype'")
				MISSMSG `anything'

			}

		}
		else {
			predictnl double `anything' = predict(`STAT' `xbtype' timevar(`timevar')) if `touse', ci(`anything'_lci `anything'_uci)
		}
		
		
        // The following predictions use empirical Bayes' estimates.

//         local log = "`log'" != ""

//         if _caller() < 14.2 {
//                 capture checkestimationsample
//                 if c(rc) {
//                         di as err "{p 0 0 2}"
//                         di as err "data have changed since estimation;{break}"
//                         di as err ///
// "prediction of empirical Bayes `EBTYPE' requires the original " ///
// "estimation data"
//                         di as err "{p_end}"
//                         exit 459
//                 }
//                 tempname esample
//                 quietly gen byte `esample' = e(sample)
//         }
//         else {
//                 local esample : copy local touse
//         }

//         mata: _gsem_predict_latent(     "`TNAME'",              ///
//                                         "double `TNAME'_*",     ///
//                                         "",                     ///
//                                         "`esample'",            ///
//                                         "`EBTYPE'",             ///
//                                         "",                     ///
//                                         "`intpoints'",          ///
//                                         "`iterate'",            ///
//                                         "`tolerance'",          ///
//                                         `log')
//         if _caller() < 14.2 {
// 			capture assert `touse' == `esample'
// 			if c(rc) {
// 				FILL `touse' `VARLIST'
// 			}
//         }
//         local LATENT : copy local VARLIST

//         if "`sortlist'" != "" {
//                 sort `sortlist'
//         }


//         if "`STAT'" == "survival" {
//                 mata: _gsem_predict_lsurv("`TNAME'",            ///
//                                         "`anything'",           ///
//                                         "`touse'",              ///
//                                         `offset',               ///
//                                         "`outcome'",            ///
//                                         "`LATENT'")
//                 MISSMSG `VARLIST'
//                 exit
//         }

	//====================================================================================================================================================//
	// Restore original data and merge in new variables 
	
		
		local keep `anything'

		if "`ci'" != "" { 
			local keep `keep' `anything'_lci `anything'_uci 
		}
		
		keep `keep'
		qui save `newvars'
		restore
		merge 1:1 _n using `newvars', nogenerate noreport

end

program MISSMSG
        tempname touse
        quietly gen byte `touse' = 1
        markout `touse' `0'
        quietly count if !`touse'
        if r(N) {
                di as txt "(`r(N)' missing values generated)"
        }
end

program FILL
        gettoken touse 0 : 0
        foreach var of local 0 {
                local gvars : char `var'[gvars]
                if "`gvars'" != "" {
                        quietly bysort `gvars' (`var') : ///
                        replace `var' = `var'[1] if `touse'
                }
                quietly replace `var' = . if !`touse'
        }
end

exit
