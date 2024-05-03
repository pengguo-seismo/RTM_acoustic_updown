
        real,allocatable::syn_p(:,:),noise(:,:)
        integer::seed,seed_tmp
        real::mu,sigma

        integer::num_shot,nrel,nx,nz,nt,nx_rec,nt_snap
        integer::nx_p,nz_p
        integer::nx_p1,nx_p2,nz_p1,nz_p2
        real(fp_kind)::deltax,deltaz,deltat
        real(fp_kind)::f0,f_ref
        logical::lossless        

        character*90::vpname,rhoname,tau_epsname,tau_sigmaname,&
                      file_p_part,file_p

        seed_tmp=20171123

        open(12,file='parameter.txt')
          read(12,*) num_shot
!!    number of relaxation mechanisms
          read(12,*) nrel
          print *,nrel
!!    grid number -x
          read(12,*) nx_p
!!    grid number -z
          read(12,*) nz_p
!!    time step 
          read(12,*) nt
          read(12,*) nx_rec
          read(12,*) lossless
!!    snapshot output interval (time step)
          read(12,*) nt_snap
!!    grid increment - x 
          read(12,*) deltax
          print *,deltax
!!    grid increment - z
          read(12,*) deltaz
          print *,deltaz
!!    time increment
          read(12,*) deltat
          print *,deltat
!!    dominant frequency
          read(12,*) f0
          print *,f0
!!    reference frequency
!!    for velocity
          read(12,*) f_ref
!!      velocity file name
          read(12,*) vpname
          print *,vpname,'vpname'
!!     density file name
          read(12,*) rhoname
          print *,rhoname,'rhoname'
!!     strain relaxation time file name (Q)
          read(12,*) tau_epsname
          print *,tau_epsname,'tau_epsname'
!!     stress relaxation times file name (Q)
          read(12,*) tau_sigmaname
          print *,tau_sigmaname,'tau_epsname'
!!     seismogram file name
          read(12,*) file_p_part
          print *,file_p_part,'seis_part1'

        close(12)

        allocate(syn_p(1:nx_rec,1:nt),noise(1:nx_rec,1:nt))

        do i_shot=1,num_shot

           file_p=trim(file_p_part)//char(i_shot/100+48)//&
                  char(mod(i_shot,100)/10+48)//char(mod(i_shot,10)+48)

           open(24,file=file_p,access='direct',form='unformatted',&
                   recl=ii_bi*(nx_rec))
             do it=1,nt
                read(24,rec=it) (syn_p(i,it),i=1,nx_rec)
             enddo
           close(24)

           seed=seed_tmp+(i_shot-1)*10

           do j=1,nt
              do i=1,nx_rec
                 noise(i,j)=r4_normal_01 ( seed )
             enddo 
           enddo 

           syn_p(:,:)=syn_p(:,:)+noise(:,:)

         enddo 

        end


