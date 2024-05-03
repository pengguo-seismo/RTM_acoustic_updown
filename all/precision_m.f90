
! This module provides a simple facility for changing between single
! and double precision
 
     module precision_m
       integer, parameter, public :: singlePrecision = kind(0.0) ! Single precision
       integer, parameter, public :: doublePrecision = kind(0.0d0) ! Double precision

! Comment out one of the lines below
       integer, parameter, public :: fp_kind = singlePrecision
!integer, parameter, public :: fp_kind = doublePrecision
     end module precision_m
