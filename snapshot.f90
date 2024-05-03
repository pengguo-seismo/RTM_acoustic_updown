
!!  output wavefield snapshot every nt_snap time step  !!

! --------------------------------------------------------------------

        subroutine snap_shot(p,it,mx,mz,lv,nt_snap,snap_part1,file_snap_no)

           implicit none
           integer::it,nt_snap,iit,file_snap_no
           integer::mx,mz,lv
           real::p(-lv:mx+lv,-lv:mz+lv)
           character *1::snap1,snap2,snap3,snap4,snap5
           character *30::file_snap
           character *30::snap_part1
           integer::i,j
           integer::ii_bi

           ii_bi=1

           if(mod(it,nt_snap).eq.0) then
             iit=int(it/10000)
             snap1=char(iit+48)
             iit=int(mod(it,10000)/1000)
             snap2=char(iit+48)
             iit=int(mod(it,1000)/100)
             snap3=char(iit+48)
             iit=int(mod(it,100)/10)
             snap4=char(iit+48)
             iit=int(mod(it,10))
             snap5=char(iit+48)
             file_snap=trim(snap_part1)//snap1//snap2//snap3//snap4//snap5

             print *,'mx',mx,'mz',mz

             open(file_snap_no,file=file_snap,access='direct',&
                  form='unformatted',recl=ii_bi*mz)
              do i=1,mx
                  write(file_snap_no,rec=i) (p(i,j),j=1,mz)
              enddo
             close(file_snap_no)

           endif

        end subroutine snap_shot

