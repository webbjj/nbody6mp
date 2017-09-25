program zpars
!A program to check to make sure each star keeps its zpars

implicit none

!       Variables for timestep header
        real :: tphys, rc, rbar, rtide, xc(1:3), zmbar, turn, rscale, tphys2
        integer :: ns, nc, nb

!       Variables for stars
        integer :: nm1, kw1, nmb1, nmb2, kwb1, kwb2, kcm
        real :: ms1,logl1,logrl,xd(1:3),vd(1:3),dum1,dum2
        real :: ecc,pb,semi,msb1,msb2,loglb1,loglb2,logrb1,logrb2,xdb(1:3),vdb(1:3),dum3,dum4
        real zpars11,zpars11b

        integer :: i, j, k, Nin
        integer,parameter :: N1=250, N2=250
        real :: RTIDE0,R,MFRAC
        real :: zpars0(N1+N2)


        open (unit=1,file="fort.83",status="old")
        open (unit=2,file="fort.82",status="old")

        do i=1,10000

                !Read in the timestep header
                read(1,*,end=100) ns, tphys
                print *,i,tphys,ns
                read(1,*) nc, rc, rbar, rtide, xc(1:3)
                read(1,*) zmbar, turn, rscale

                do j=1,2*ns !Make sure you get all the stars

                        !Read in the star info. Sort it according to its id and
                        !print into the right fort.83
        
                        read(1,*) nm1,kw1,ms1,logl1,logrl,xd(1),xd(2),xd(3),vd(1),vd(2),vd(3),dum1,dum2,zpars11

                        if (nm1.EQ.-1000) then
!                                print *,'FOUND IT: ',j,nm1
                                go to 10

                        end if

                        if (i.EQ.1)then
                          zpars0(nm1)=zpars11
!                          print *,j,nm1,zpars0(nm1)
                        else
                          if(zpars11.ne.zpars0(nm1)) then
                            print *,'SINGLE ERROR ',i,j,nm1,ms1,zpars11,zpars0(nm1)
                          end if
                        end if

!                        if (i.eq.1 .and. j.ne.nm1) print *,'LOST STAR ',j,nm1

5                       continue

                end do

10              continue

                !Now repeat the process for the fort.82

                !Read in the timestep header

                read(2,*) nb,tphys2
                print *,i,tphys2,nb

                if (nb.EQ.0) then

                        read(2,*)nmb1,nmb2,kwb1,kwb2,kcm,ecc,pb,semi,msb1,msb2,loglb1,loglb2,logrb1,logrb2,xdb(1),xdb(2),xdb(3),vdb(1),vdb(2),vdb(3),dum1,dum2,dum3,dum4,zpars11,zpars11b

                else

                        do j=1,2*nb

                                read(2,*) nmb1,nmb2,kwb1,kwb2,kcm,ecc,pb,semi,msb1,msb2,loglb1,loglb2,logrb1,logrb2,xdb(1),xdb(2),xdb(3),vdb(1),vdb(2),vdb(3),dum1,dum2,dum3,dum4,zpars11,zpars11b

                                if (nmb1.EQ.-1000) then
!                                        print *,'FOUND THE BINARY: ',j,nmb1
                                        go to 20
                                end if

                                if (i.EQ.1)then
                                  zpars0(nmb1)=zpars11
                                  zpars0(nmb2)=zpars11b
                                  print *,nmb1,zpars0(nmb1)
                                  print *,nmb1,zpars0(nmb1)
                                else
                                  if(zpars11.ne.zpars0(nmb1)) then
                                    print *,'BINARY ERROR ',i,j,nmb1,msb1,zpars11,zpars0(nmb1)
                                  end if
                                  if(zpars11b.ne.zpars0(nmb2)) then
                                    print *,'BINARY ERROR ',i,j,nmb2,msb2,zpars11b,zpars0(nmb2)
                                  end if
                                end if

15                      continue
        
                        end do
               
                 end if

20              continue

        end do

        close(1)
        close(2)

100     continue

end program
