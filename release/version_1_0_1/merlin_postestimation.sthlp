{smcl}
{* *! version 1.0.0  ?????2018}{...}
{vieweralsosee "merlin" "help merlin"}{...}
{viewerjumpto "Syntax" "merlin_postestimation##syntax"}{...}
{viewerjumpto "Description" "merlin_postestimation##description"}{...}
{viewerjumpto "Options" "merlin_postestimation##options"}{...}
{viewerjumpto "Remarks" "merlin_postestimation##remarks"}{...}
{viewerjumpto "Examples" "merlin_postestimation##examples"}{...}

{marker syntax}{...}
{title:Syntax for predict}

{pstd}
Syntax for predictions following a {helpb merlin:merlin} model

{p 8 16 2}
{cmd:predict}
{it:newvarname}
{ifin} [{cmd:,}
{it:{help merlin_postestimation##statistic:statistic}}
{it:{help merlin_postestimation##opts_table:options}}]


{phang}
The default is to make predictions based only on the fixed portion of the 
model.  

{marker statistic}{...}
{synoptset 22 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{syntab:Main}
{synopt :{opt mu}}expected value of {depvar}; the default{p_end}
{synopt :{opt eta}}expected value of linear prediction of {depvar}{p_end}
{synopt :{opt surv:ival}}survivor function at {depvar}{p_end}
{synopt :{opt f:ailure}}failure probability at {depvar}{p_end}
{synopt :{opt h:azard}}hazard function at {depvar}{p_end}
{synopt :{opt ch:azard}}cumulative hazard function at {depvar}{p_end}
{synopt :{opt rmst}}restricted mean survival time at {depvar}{p_end}
{synopt :{opt rmft}}restricted mean failure/event time at {depvar}{p_end}
{synoptline}

{marker opts_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr:options}
{synoptline}
{syntab:Main}
{synopt :{opt fixedonly}}compute {it:statistic} based only on the fixed portion of the model; the default{p_end}
{synopt :{opt marginal}}compute {it:statistic} marginally with respect to the latent variables{p_end}
{synopt :{cmd:outcome(}{it:#}{cmd:)}}specify observed response variable (default 1){p_end}
{synopt :{opt at(at_spec)}}specify covariate values for prediction{p_end}
{synopt :{opt ci}}calculate confidence intervals{p_end}
{synopt :{cmd:timevar(}{varname}{cmd:)}}calculate predictions at specified time-points{p_end}

{syntab :Integration}
{synopt :{opt intp:oints(#)}}use
        {it:#} integration points to compute marginal predictions {p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:predict} is a standard postestimation command of Stata.
This entry concerns use of {cmd:predict} after {cmd:merlin}.

{pstd}
{cmd:predict} after {cmd:merlin} creates new variables containing
observation-by-observation values of estimated observed response variables,
linear predictions of observed response variables, or other such functions.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{cmd:mu}, the default, calculates the expected value of the outcomes.

{phang} 
{cmd:eta} calculates the fitted linear prediction.

{phang} 
{cmd:survival} calculates the survival function.

{phang} 
{cmd:failure} calculates the failure probability, which is 1 - survival.

{phang} 
{cmd:hazard} calculates the hazard function.

{phang}
{cmd:chazard} calculates the cumulative hazard function.

{phang} 
{cmd:rmst} calculates the restricted mean survival time, which is the integral of the survival functions.

{phang} 
{cmd:rmft} calculates the restricted mean failure/event time, which is the integral of 1 minus the survival functions.

{phang}
{cmd:fixedonly} specifies that the predicted {it:statistic} be computed
based only on the fixed portion of the model. This is the default.

{phang}
{cmd:marginal} specifies that the predicted {it:statistic} be computed
marginally with respect to the latent variables.

{phang2}
Although this is not the default, marginal predictions are often very useful
in applied analysis.  They produce what are commonly called
population-averaged estimates. 

{phang2}
For models with continuous latent variables, the {it:statistic} is calculated
by integrating the prediction function with respect to all the latent
variables over their entire support.

{phang}
{cmd:outcome(}{it:#}{cmd:)} specifies that predictions for
outcome {it:#} be calculated.

{phang}
{cmd:ci} specifies that confidence intervals are calculated for the predicted {it:statistic}. They will 
stored in {it:newvarname_lci} and {it:newvarname_uci}.

{phang}
{cmd:timevar(}{varname}{cmd:)}calculate predictions at specified time-points. 
For survival models, the default is to calculate predictions at the response times. 
For a {cmd:merlin} model where a {cmd:timevar()} was specified, then the default will use the original 
{cmd:timevar()}. This option overides it.{p_end}

{dlgtab:Integration}

{phang}
{opt intpoints(#)} specifies the number of integration points used to
compute marginal predictions; the default is the value from estimation.


{marker remarks}{...}
{title:Remarks}

{pstd}
Out-of-sample prediction is allowed for all {cmd:predict} options.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. use http://fmwww.bc.edu/repec/bocode/s/stjm_pbc_example_data, clear}{p_end}

{pstd}Linear mixed effects model with {cmd:merlin}{p_end}
{phang2}{cmd:. merlin (logb time age trt time#M1[id]@1 M2[id]@1, family(gaussian))}{p_end}

{pstd}Predict the expected value of {cmd:logb} marginalised over the random effects{p_end}
{phang2}{cmd:. predict ev1, eta marginal}{p_end}

