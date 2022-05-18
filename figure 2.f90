       !      p vs beta in mutation free model

      integer t, i, j

      real x, y, z, s1, beta, delta, xi, k, h, k11, k12, k13, &
      k21, k22, k23, k31, k32, k33, k41, k42, k43, x1, y1, z1, tol, &
      x0, y0, z0, alpha, miu , c0, c1, c2, c3, c4, c5, c6, c7, fp1, &
      fp2,fp3,p,eta,s2,s3,sig1, sig2




       open(11,file='x8.dat')
       open(22,file='y8.dat')
       open(33,file='xy8.dat')
       open(44,file='origin8.dat')
       open(100,file='ecoevo_ts8.dat')



	xi = 0.15
        !p=0.30 !0.42
        miu=0.0 !0.35
        sig1 = 0.75
        sig2 = 0.25
        !beta= 1.1
        eta=0.85



	do i = 0, 100
	   p = 0.0+(1.0/100.0)*i
	do j = 1, 99
           beta = 1.0+(1.0/100.0)*j

        write(*,*) i,j



	x = 0.35
	y = 0.35


	h = 0.01


	do 200 t = 1, 2000000

	k = x*((1-sig1)*x*(1-x-y)-sig1*y*(1-x-y)+sig1*(1-x-y)-xi)+miu*(y-x)
	k11 = k*h
	k = y*((beta-sig2)*x*(1-x-y)+(2*p*eta-eta-sig2)*y*(1-x-y)+sig2*(1-x-y)-xi)+miu*(x-y)
	k12 = k*h


     	x1 = x + k11/2.0
     	y1 = y + k12/2.0



    	k = x*((1-sig1)*x*(1-x-y)-sig1*y*(1-x-y)+sig1*(1-x-y)-xi)+miu*(y-x)
	k21 = k*h
	k = y*((beta-sig2)*x*(1-x-y)+(2*p*eta-eta-sig2)*y*(1-x-y)+sig2*(1-x-y)-xi)+miu*(x-y)
	k22 = k*h


     	x1 = x + k21/2.0
     	y1 = y + k22/2.0


    	k = x*((1-sig1)*x*(1-x-y)-sig1*y*(1-x-y)+sig1*(1-x-y)-xi)+miu*(y-x)
	k31 = k*h
	k = y*((beta-sig2)*x*(1-x-y)+(2*p*eta-eta-sig2)*y*(1-x-y)+sig2*(1-x-y)-xi)+miu*(x-y)
	k32 = k*h


     	x1 = x + k31
     	y1 = y + k32



   	k = x*((1-sig1)*x*(1-x-y)-sig1*y*(1-x-y)+sig1*(1-x-y)-xi)+miu*(y-x)
	k41 = k*h
	k = y*((beta-sig2)*x*(1-x-y)+(2*p*eta-eta-sig2)*y*(1-x-y)+sig2*(1-x-y)-xi)+miu*(x-y)
	k42 = k*h


     	x = x + k11/6.0 + k21/3.0 + k31/3.0 + k41/6.0
     	y = y + k12/6.0 + k22/3.0 + k32/3.0 + k42/6.0



       if(t.eq.2000000) then
       print*,'******************************************'
       write(*,*) p,beta,x,y,x+y
       print*,'******************************************'
       write(100,*) p,beta,x,y,x+y


	if(x+y .le. 1.0 .and. x+y .ge. 0.0) then

       if(x .lt. 0.0001 .and. y .lt. 0.0001)then
       if(x .ge. 0.0 .and. y .ge. 0.0) then
       write(44,*) p,beta
       end if
       end if


       if(y .lt. 0.0001 .and. y .ge. 0.0)then
       if(x .gt. 0.0001) then
       write(11,*) p,beta
       end if
       end if


       if(x .lt. 0.0001 .and. x .ge. 0.0)then
       if(y .gt. 0.0001) then
       write(22,*) p,beta
       end if
       end if


       if(y .gt. 0.0001 .and. x .gt. 0.0001) then
       write(33,*) p,beta
       end if




	end if



       end if




200   continue

     end do
      end do



      stop
      end
