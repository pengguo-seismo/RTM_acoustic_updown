program adder

INTERFACE
  SUBROUTINE addnums(a, b) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
    IMPLICIT NONE
    INTEGER(C_INT) :: a, b
  END SUBROUTINE addnums
END INTERFACE

integer a,b
a=1
b=2
call addnums(a,b)
stop    
end program
