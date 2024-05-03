#!/bin/bash

rm *.o
F90=mpif90
mp=-qopenmp

DUG='-g -traceback -check all -check bounds'

module load openmpi
module load fftw

$F90 -c precision_m.f90
$F90 -c main_module.f90
$F90 -c wave_module.f90
$F90 -c parameter_input.f90
$F90 -c wave_module.f90
$F90 -c wave_drec_sep_module.f90
$F90 -c wave_hilbert_module.f90
$F90 -c fftw_module_f.f90
$F90 -c wave_propagator.f90 $mp
$F90 -c pad.f90
$F90 -c utili.f90
$F90 -c vp_f0_2_vp.f90
$F90 -c assign_parameter.f90
$F90 -c adding_abs_val.f90 
$F90 -c wavlet.f90
$F90 -c ricker_w_hb.f90
$F90 -c snapshot.f90
$F90 -c cpml_coef.f90
$F90 -c imaging_condition.f90 
$F90 -c imaging_condition_udlr_eff.f90 
$F90 -c bvr_record_wave_p.f90
$F90 -c bvr_record_wave_vxz.f90
$F90 -c hilbert_seismo.f90
$F90 -c visco_acous_ud_eff_hilbert_t.f90 $mp 
$F90 *.o -o rtm_ud_eff_t_p $mp  -lfftw3f -Wl,-rpath=$LD_RUN_PATH
