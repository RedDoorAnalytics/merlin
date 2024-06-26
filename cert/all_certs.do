//============================================================================//
// cert. script for merlin

//source paths
local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/stmerlin/stmerlin"


//build mlib
clear all
do ./build/buildmlib.do
mata mata clear

//============================================================================//

do ./cert/ologit-oprobit

//survival
do ./cert/streg.do
do ./cert/strcs.do
do ./cert/stpm2.do
do ./cert/pwexp.do
do ./cert/cert_cox.do
do ./cert/cert_cox_firth.do
do ./cert/relsurverror.do
do ./cert/userhaz.do
do ./cert/cert_bhazard_comptostpm2.do

//gf2
do ./cert/cert_streg_cox_gf2.do
do ./cert/cert_rp_rcs_gf2.do
do ./cert/cert_gauss_gf2.do
do ./cert/cert_poisson_gf2.do
do ./cert/cert_rp_intcens.do
do ./cert/cert_haz_loghaz_rp_bhazard.do

//multilevel survival
do ./cert/cert_mestreg.do
do ./cert/weibull_3level.do
do ./cert/jf1.do

//longitudinal
do ./cert/level1.do

//joint longitudinal-survival
do ./cert/bivariate_jm.do
do ./cert/weighted_cumulative_jm.do
do ./cert/jm_ltruncated.do


//============================================================================//

di "certs complete - no errors"

//============================================================================//
