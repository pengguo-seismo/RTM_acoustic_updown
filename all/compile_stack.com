
rm *.o
ftn="mpif90 -O2 -xHOST -fPIC"
mp=-qopenmp
mp=-fopenmp

module load openmpi

$ftn -c stack_image.f90
$ftn -c filter.f90
$ftn *.o -o stack
