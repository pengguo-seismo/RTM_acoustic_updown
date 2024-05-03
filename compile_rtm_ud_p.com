#!/bin/bash

rm *.o
F90=mpif90
mp=-qopenmp
#module load fftw

opt1='-axCORE-AVX2'

module load openmpi/4.1.4-ofed54-intel20

DUG='-g -traceback -check all -check bounds -stand f03 -assume realloc_lhs -fstack-protector  -assume protect_parens'

$F90 -c precision_m.f90 $opt1
$F90 -c main_module.f90 $opt1
$F90 -c wave_module.f90 $opt1
$F90 -c fftw_module_f.f90  $opt1
$F90 -c parameter_input.f90 $opt1
$F90 -c wave_module.f90  $opt1
$F90 -c wave_drec_sep_module.f90 $opt1
$F90 -c wave_hilbert_module.f90 $opt1
$F90 -c wave_propagator.f90 $mp $opt1
$F90 -c pad.f90 $opt1
$F90 -c utili.f90 $opt1
$F90 -c vp_f0_2_vp.f90 $opt1
$F90 -c assign_parameter.f90 $opt1
$F90 -c adding_abs_val.f90  $opt1
$F90 -c wavlet.f90 $opt1
$F90 -c ricker_w_hb.f90 $opt1
$F90 -c snapshot.f90 $opt1
$F90 -c cpml_coef.f90 $opt1
$F90 -c imaging_condition.f90  $opt1
$F90 -c imaging_condition_udlr.f90  $opt1
$F90 -c imaging_condition_udlr_im.f90 $opt1
$F90 -c imaging_condition_udlr_eff.f90  $opt1
$F90 -c bvr_record_wave_p.f90 $opt1
$F90 -c bvr_record_wave_vxz.f90 $opt1
$F90 -c hilbert_seismo.f90 $opt1
$F90 -c wavefield_decomp.f90 $opt1
$F90 -c visco_acous_ud_hilbert.f90 $mp $opt1
$F90 *.o -o rtm_ud_p $mp  -lfftw3f $opt1 -Wl,-rpath=$LD_RUN_PATHi 



