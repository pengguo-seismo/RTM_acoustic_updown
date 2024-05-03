
rm *.o
ftn="mpif90"

mp=-qopenmp

opt1='-axCORE-AVX2'

#opt1='-axCORE-AVX2 -g -traceback -check all -check bounds'

#module load openmpi/4.1.1-ofed51-intel20

module load openmpi/4.1.4-ofed54-intel20

$ftn -c precision_m.f90 $opt1
$ftn -c main_module.f90 $opt1
$ftn -c wave_module.f90 $opt1
$ftn -c parameter_input.f90 $opt1
$ftn -c wave_propagator.f90 $mp $opt1
$ftn -c adding_abs_val.f90 $opt1
$ftn -c vp_f0_2_vp.f90 $opt1
$ftn -c assign_parameter.f90 $opt1
$ftn -c wavlet.f90 $opt1
$ftn -c snapshot.f90 $opt1
$ftn -c cpml_coef.f90 $opt1
$ftn -c imaging_condition.f90 $opt1 
$ftn -c bvr_record_wave_p.f90  $opt1
$ftn -c bvr_record_wave_vxz.f90 $opt1
$ftn -c visco_acous_pstm.f90 $mp $opt1
$ftn *.o -o rtm_pstm_p $mp $opt1
