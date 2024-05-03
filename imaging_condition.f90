
!! apply cross-correlation imaging condition

! --------------------------------------------------------------------

       subroutine imaging_condition(sou_p,rcv_p,image,mx,mz,lx,lz)

          implicit none
          integer::mx,mz,lx,lz
          real::sou_p(-lx:mx+lx,-lz:mz+lz),rcv_p(-lx:mx+lx,-lz:mz+lz),&
                image(1:mx,1:mz)
          integer::i,j

          do i=1,mx
             do j=1,mz
                image(i,j)=image(i,j)+sou_p(i,j)*rcv_p(i,j)
             enddo 
          enddo 

        end subroutine imaging_condition
          
