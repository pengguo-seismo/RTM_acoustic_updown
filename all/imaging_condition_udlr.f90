
!! apply cross-correlation imaging condition

! --------------------------------------------------------------------

       subroutine imaging_condition_udlr(sou_p_up,sou_p_dw,&
        rcv_p_up,rcv_p_dw,image_du,image_ud,image_dd,image_uu,&
        image_duud,image_dduu,nx,nz,lv)

          implicit none
          integer::nx,nz,lv
          real,dimension(-lv:nx+lv,-lv:nz+lv)::sou_p_up,sou_p_dw,&
                         rcv_p_up,rcv_p_dw
          real,dimension(1:nx,1:nz)::image_ud,image_du,image_duud
          real,dimension(1:nx,1:nz)::image_uu,image_dd,image_dduu
          integer::i,j

          do j=1,nz
             do i=1,nx
                image_du(i,j)=image_du(i,j)+sou_p_dw(i,j)*rcv_p_up(i,j)
                image_ud(i,j)=image_ud(i,j)+sou_p_up(i,j)*rcv_p_dw(i,j)
                image_uu(i,j)=image_uu(i,j)+sou_p_up(i,j)*rcv_p_up(i,j)
                image_dd(i,j)=image_dd(i,j)+sou_p_dw(i,j)*rcv_p_dw(i,j)
             enddo 
          enddo 

          do j=1,nz
             do i=1,nx
                image_duud(i,j)=image_duud(i,j)+sou_p_dw(i,j)*rcv_p_up(i,j)&
                               +sou_p_up(i,j)*rcv_p_dw(i,j)
                image_dduu(i,j)=image_dduu(i,j)+sou_p_up(i,j)*rcv_p_up(i,j)&
                               +sou_p_dw(i,j)*rcv_p_dw(i,j)
             enddo 
          enddo 

        end subroutine imaging_condition_udlr
          
