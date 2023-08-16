*! version 0.1.8 23apr2018 MJC

/*
History
23apr2018: version 1.0.0 - merlin released to SSC
						 - megenreg -> merlin
						 - bug fix with if statements - merlin_modeltouses_post_xb() creating tousei was giving missings not 0s. Now fixed.
						 - offset() added to fp() and rcs()
						 - error check added for . notation or ## in CLP
						 - handling of elements re-written -> @ only for constraints
						 - predict function added - fixedonly and margina allowed
						 - results help file
						 - EV now expected value of response, added XB for expected value of linear predictor
						 - random effects must be M#[]
						 - documented all utility functions
						 - results table greatly improved
						 - interactions between fp() and rcs() etc now allowed
16dec2017: version 0.1.8 - family(gamma) added
						 - bug fix in inverse links for EV functions
20nov2017: version 0.1.7 - bug fix in not documented links when using family(user)
						 - mf() element added = user-defined mata function
						 - family(rcs) added
						 - startval for @s changed from 0 to 0.0001
						 - ereturn list additions
27oct2017: version 0.1.6 - fixed output of bernoulli model, labelled equation cloglog when it should be logit (was still using logit)
						 - random effects in results table now transformed to sds and corrs, display much improved
						 - replay fixed
23oct2017: version 0.1.5 - added family(lquantile [, quantile(#)]) for linear quantile regression
						 - stable added to sort operations
20oct2017: version 0.1.4 - default cumulative hazard integration changed to 15-point Gauss-Kronrod
						 - loghfunction() added to family(user)
						 - improved error checks
19oct2017: version 0.1.3 - starting value for @ parameters changed to 0.1 in prev. release; put back to 0s
						 - minor improvements to the mlib
15oct2017: version 0.1.2 - constraints() fixed, family(null) now working with ob. level models
12oct2017: version 0.1.1 - restartvalues() added
08oct2017: version 0.1.0 - dev. release
*/

/*
Familys

exponential
Weibull
Gompertz
RP
rcs
Gaussian
Bernoulli
Poisson
beta
negative binomial
ordered logit
ordered probit
lquantile
gamma
*/

program merlin, eclass 
        version 14.2

        if replay() {
			if "`e(cmd)'" != "merlin" {
				error 301
			}
			Display `0'
			exit
        }
		
		tempname GML
        capture noisily Estimate `GML' `0'
        local rc = c(rc)
		capture mata: rmexternal("`GML'")
        capture drop `GML'*
        if (`rc') exit `rc'
        ereturn local cmdline `"merlin `0'"'
end

program Estimate, eclass
        version 14.2
        gettoken GML : 0

        `vv' `BY' Fit `0'		//!! should leave behind diopts

		mata: merlin_ereturn("`GML'")
		
        Display, `diopts'
end

program Fit, eclass sortpreserve
        version 14.2
        gettoken GML 0 : 0

        tempname touse b
        merlin_parse `GML', touse(`touse') : `0'
		if "`r(predict)'"!="" {
			exit
		}
        local mltype    `"`r(mltype)'"'
        local mleval    `"`r(mleval)'"'
        local mlspec    `"`r(mlspec)'"'
        local mlopts    `"`r(mlopts)'"'
        local mlvce     `"`r(mlvce)'"'
        local mlwgt     `"`r(mlwgt)'"'
		local mlcns		`"`r(constr)'"'
		local mlinitcns	`"`r(initconstr)'"'
        local nolog     `"`r(nolog)'"'
        local mlprolog  `"`r(mlprolog)'"'
		local mftodrop  `r(mftodrop)'
        c_local diopts  `"`r(diopts)'"'
		local mlfrom  `"`r(mlfrom)'"'
		local mlzeros	`"`r(mlzeros)'"'
		if "`mlcns'" != "" {
			local cnsopt constraint(`mlcns')
        }
		
        if "`mlfrom'"!="" {
			matrix `b' = r(b)
			local mlinit init(`b',copy)
		}
		
		di
		di as txt "Fitting full model:"

		cap n ml model `mltype' `mleval'              			///
								`mlspec'                        ///
								`mlwgt'                         ///
								if `touse',                     ///
								`mlopts'                        ///
								`mlvce'                         ///
								`mlprolog'                      ///
								`cnsopt'                        ///
								`mlinit'						///
								collinear                       ///
								maximize     	                ///
								missing    		                ///
								nopreserve   	                ///
								search(off)                     ///
								userinfo(`GML')	                ///
								wald(0)   	                    
														 
			if _rc>0 {
				if _rc==1400 & "`mlzeros'"=="" {
				
				di ""
				di as text "-> Starting values failed - trying zero vector"
			
				ml model `mltype' 	`mleval'              			///
									`mlspec'                        ///
									`mlwgt'                         ///
									if `touse',                     ///
									`mlopts'                        ///
									`mlvce'                         ///
									`mlprolog'                      ///
									`cnsopt'                        ///
									collinear                       ///
									maximize     	                ///
									missing    		                ///
									nopreserve   	                ///
									search(off)                     ///
									userinfo(`GML')	                ///
									wald(0)   	                    
				}
				else {
					exit _rc
				}
			}
		
        ereturn local title "Mixed effects regression model"
        ereturn local predict   merlin_p
		ereturn local cmd merlin
		cap mata: mata drop `mftodrop'
end

program Display
        syntax [,       noHeader        ///
                        noDVHeader      ///
                        noLegend        ///
                        notable         ///
                        *               ///
        ]

        _get_diopts diopts, `options'
        if e(estimates) == 0 {
                local coefl coeflegend selegend
                local coefl : list diopts & coefl
                if `"`coefl'"' == "" {
                        local diopts `diopts' coeflegend
                }
        }
		
		local Nrelevels = `e(Nlevels)'-1
		
		local plus
		if "`e(Nres1)'"!="" {
			local plus plus
			local neq = `e(k_eq)'
			forval l=1/`Nrelevels' {
				local neq = `neq' - `e(Nreparams`l')'
			}
			local neq neq(`neq')
		}
		
		_coef_table_header
		_coef_table, neq(0) plus nocnsreport
		
		forval mod=1/`e(Nmodels)' {
		
			local y : word 1 of `e(response`mod')'
			if "`y'"=="" {
				local y null
			}
			_diparm __lab__, label("`y':") eqlabel
			
			local Ncmps : word count `e(Nvars_`mod')'
			
			forval c = 1/`Ncmps' {
			
				local clab : word `c' of `e(cmplabels`mod')'
				local np : word `c' of `e(Nvars_`mod')'
				if `np'>1 {
					forval el=1/`np' {
						_diparm _cmp_`mod'_`c'_`el', label("`clab':`el'")
					}
				}
				else {
					_diparm _cmp_`mod'_`c'_1, label("`clab'")
				}
				
			}
			
			if `e(constant`mod')' {
				_diparm cons`mod', label("_cons")
			}
		
			if `e(ndistap`mod')'>0 & "`e(family`mod')'"!="rcs" & "`e(family`mod')'"!="rp" {
				
				if "`e(family`mod')'"=="weibull" {
					_diparm dap`mod'_1, label("log(gamma)") 
				}
				else if "`e(family`mod')'"=="gompertz" {
					_diparm dap`mod'_1, label("gamma") 
				}
				else if "`e(family`mod')'"=="beta" {
					_diparm dap`mod'_1, label("log(s)") 
				}
				else if "`e(family`mod')'"=="negbinomial" {
					_diparm dap`mod'_1, label("log(alpha)") 
				}
				else if "`e(family`mod')'"=="gamma" {
					_diparm dap`mod'_1, label("log(s)") 
				}
				
				else if "`e(family`mod')'"=="gaussian" | "`e(family`mod')'"=="lquantile" {
					_diparm dap`mod'_1, label("sd(resid.)") exp
				}
				else if "`e(family`mod')'"=="ordinal" {
					forval a=1/`e(ndistap`mod')' {
						_diparm dap`mod'_`a', label("cut`a'")
					}
				}
				else {
					forval a=1/`e(ndistap`mod')' {
						_diparm dap`mod'_`a', label("dap:`a'")
					}
				}
			}
		
			if `e(nap`mod')'>0 {
				forval a=1/`e(nap`mod')' {
					_diparm ap`mod'_`a', label("ap:`a'")
				}
			}
		
			if `mod'==`e(Nmodels)' & `e(Nlevels)'==1 {
				_diparm __bot__	
			}
			else {
				_diparm __sep__	
			}

		}		
		
		
		//VCV display
		if "`e(Nres1)'"!="" {
			
			forval i=1/`Nrelevels' {
				local lev : word `i' of `e(levelvars)'
				_diparm __lab__ , label("`lev':") eqlabel
				
				forval j=1/`e(Nreparams`i')' {
					local param : word `j' of `e(re_eqns`i')'
					local scale : word `j' of `e(re_ivscale`i')'
					local label : word `j' of `e(re_label`i')'
					_diparm `param', `scale' label(`label')
				}
				if `i'<`Nrelevels' {
					_diparm __sep__
				}
			}
			
			_diparm __bot__
		}
		
end

exit
