
rm *.o
ftn="mpif90 -O2 -xHOST -fPIC"
mp=-fopenmp
module load openmpi

$ftn -c precision_m.f90
$ftn -c main_module.f90
$ftn -c wave_module.f90
$ftn -c parameter_input.f90
$ftn -c wave_propagator.f90 $mp
$ftn -c adding_abs_val.f90
$ftn -c vp_f0_2_vp.f90
$ftn -c assign_parameter.f90
$ftn -c wavlet.f90
$ftn -c snapshot.f90
$ftn -c cpml_coef.f90
$ftn -c imaging_condition.f90 
$ftn -c bvr_record_wave_p.f90
$ftn -c bvr_record_wave_vxz.f90
$ftn -c visco_acous_recons.f90 $mp
$ftn *.o -o recons_p $mp
