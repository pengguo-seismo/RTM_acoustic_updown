
!!  define cpml coefs   !!

! Coded by : Peng Guo
! Date : January 2014
! Language : Fortran 90
! Copyright: Center for Lithospheric Studies
!            The University of Texas at Dallas, 2014
!            TOTAL E&P USA, 2014
! --------------------------------------------------------------------

      subroutine cpml_coef(a_x,b_x,k_x,a_z,b_z,k_z,&
      a_x_half,b_x_half,k_x_half,a_z_half,b_z_half,k_z_half,&
      npoints_pml,mx,mz,dx,dz,deltat,c0,f0)

       implicit none

       real::a_x(1:mx,1:mz),b_x(1:mx,1:mz),k_x(1:mx),a_z(1:mx,1:mz),b_z(1:mx,1:mz),&
             k_z(1:mz)

       real::a_x_half(1:mx,1:mz),b_x_half(1:mx,1:mz),k_x_half(1:mx),&
             a_z_half(1:mx,1:mz),b_z_half(1:mx,1:mz),k_z_half(1:mz)

       real::c0(1:mx,1:mz)

       real::thickness_PML_x,thickness_PML_z,xoriginleft,xoriginright,&
             zoriginbottom,zorigintop,zero

       real::Rcoef,xval,zval,abscissa_in_PML,&
             abscissa_normalized,alpha_max_pml,quasi_cp_max

       real::deltat,dx,dz,pi

       integer::mx,mz,npoints_pml

       real::d_x(1:mx),d_z(1:mz),alpha_x(1:mx),alpha_z(1:mz)
       real::d_x_half(1:mx),d_z_half(1:mz),alpha_x_half(1:mx),alpha_z_half(1:mz)
       real::d0_x(1:mx,1:mz),d0_z(1:mx,1:mz)

       real::npower,k_max_pml,f0

       integer::i,j,nxpmls,nzpmls

       nxpmls=npoints_pml
       nzpmls=npoints_pml

       a_x(:,:)=0.0
       b_x(:,:)=0.0
       k_x(:)=1.0
       a_z(:,:)=0.0
       b_z(:,:)=0.0
       k_z(:)=1.0

       a_x_half(:,:)=0.0
       b_x_half(:,:)=0.0
       k_x_half(:)=1.0
       a_z_half(:,:)=0.0
       b_z_half(:,:)=0.0
       k_z_half(:)=1.0

       d_x(:)=0.0
       d_z(:)=0.0
       alpha_x(:)=0.0
       alpha_z(:)=0.0
       d_x_half(:)=0.0
       d_z_half(:)=0.0
       alpha_x_half(:)=0.0
       alpha_z_half(:)=0.0

       zero=0.0

       npower=2.0
       k_max_pml=1.0

       quasi_cp_max=c0(1,1)

       thickness_PML_x = nxpmls * DX
       thickness_PML_z = nzpmls * Dz

       Rcoef = 0.000001

       pi=4.0*atan(1.0)
       alpha_max_pml=2.0*pi*f0/2.0  
 
       do i=1,mx
          do j=1,mz
            d0_x(i,j) = - (NPOWER + 1) * c0(i,j) * log(Rcoef) / (2.0 * thickness_PML_x)
            d0_z(i,j) = - (NPOWER + 1) * c0(i,j) * log(Rcoef) / (2.0 * thickness_PML_z)
          enddo 
       enddo 

       xoriginleft = thickness_PML_x
       xoriginright = (mx-1)*DX - thickness_PML_x

       do j=1,mz
         do i = 1,mX

! abscissa of current grid point along the damping profile
             xval = DX * dble(i-1)

!---------- left edge
!          if(USE_PML_XMIN) then

! define damping profile at the grid points
             abscissa_in_PML = xoriginleft - xval
             if(abscissa_in_PML >= ZERO) then
               abscissa_normalized = abscissa_in_PML / thickness_PML_x
               d_x(i) = d0_x(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
               K_x(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
               alpha_x(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML
             endif

             abscissa_in_PML = xoriginleft - (xval+dx/2.e0)
             if(abscissa_in_PML >= ZERO) then
               abscissa_normalized = abscissa_in_PML / thickness_PML_x
               d_x_half(i) = d0_x(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
               K_x_half(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
               alpha_x_half(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML
             endif
!          endif

!---------- right edge
!          if(USE_PML_XMAX) then

! define damping profile at the grid points
             abscissa_in_PML = xval - xoriginright
             if(abscissa_in_PML >= ZERO) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x(i) = d0_x(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                K_x(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                alpha_x(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML
             endif
!          endif

             abscissa_in_PML = xval + dx/2.e0 - xoriginright
             if(abscissa_in_PML >= ZERO) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x
                d_x_half(i) = d0_x(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                K_x_half(i) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                alpha_x_half(i) = ALPHA_MAX_PML * (1.0 - abscissa_normalized) + 0.1 * ALPHA_MAX_PML
             endif

! just in case, for -5 at the end
             if(alpha_x(i) < ZERO) alpha_x(i) = ZERO
             if(alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

             b_x(i,j) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
             b_x_half(i,j) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT)

! this to avoid division by zero outside the PML
             if(abs(d_x(i)) > 1.e-6) a_x(i,j) = d_x(i) * (b_x(i,j) - 1.0) / (K_x(i) * (d_x(i) &
                                           + K_x(i) * alpha_x(i)))
             if(abs(d_x_half(i)) > 1.e-6) a_x_half(i,j) = d_x_half(i) * &
                      (b_x_half(i,j) - 1.0) / (K_x_half(i) * (d_x_half(i) &
                      + K_x_half(i) * alpha_x_half(i)))

        enddo
        enddo 

! damping in the z direction

! origin of the PML layer (position of right edge minus thickness, in meters)
       zoriginbottom = thickness_PML_z
       zorigintop = (mz-1)*Dz - thickness_PML_z

       do i=1,mx
          do j = 1,mz

! abscissa of current grid point along the damping profile
             zval = Dz * real(j-1)
!---------- top edge
!           if(USE_PML_zMIN) then
! define damping profile at the grid points
             abscissa_in_PML = zoriginbottom - zval
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z(j) = d0_z(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized) &
                            + 0.1 * ALPHA_MAX_PML
             endif
             abscissa_in_PML = zoriginbottom - (zval+dz/2.e0)
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z_half(j) = d0_z(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z_half(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z_half(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized) &
                                 + 0.1 * ALPHA_MAX_PML
             endif

!            endif

!---------- bottom edge
!           if(USE_PML_zMAX) then

! define damping profile at the grid points
             abscissa_in_PML = zval - zorigintop
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z(j) = d0_z(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized)&
                            + 0.1 * ALPHA_MAX_PML
             endif

             abscissa_in_PML = zval+dz/2.e0 - zorigintop
             if(abscissa_in_PML >= ZERO) then
                 abscissa_normalized = abscissa_in_PML / thickness_PML_z
                 d_z_half(j) = d0_z(i,j) * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
                 K_z_half(j) = 1.0 + (K_MAX_PML - 1.0) * abscissa_normalized**NPOWER
                 alpha_z_half(j) = ALPHA_MAX_PML * (1.0 - abscissa_normalized) &
                                 + 0.1 * ALPHA_MAX_PML
             endif
!           endif

             b_z(i,j) = exp(- (d_z(j) / K_z(j) + alpha_z(j)) * DELTAT)
             b_z_half(i,j) = exp(- (d_z_half(j) / K_z_half(j) &
                           + alpha_z_half(j)) * DELTAT)

! this to avoid division by zero outside the PML
             if(abs(d_z(j)) > 1.e-6) a_z(i,j) = d_z(j) * (b_z(i,j) - 1.0) &
                                     / (K_z(j) * (d_z(j) + K_z(j) * alpha_z(j)))
             if(abs(d_z_half(j)) > 1.e-6) a_z_half(i,j) = d_z_half(j) * &
                                   (b_z_half(i,j) - 1.0) / (K_z_half(j) &
                               * (d_z_half(j) + K_z_half(j) * alpha_z_half(j)))

         enddo
        enddo 

       end subroutine cpml_coef


