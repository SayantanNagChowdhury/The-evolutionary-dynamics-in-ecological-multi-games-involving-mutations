        ! Bifurcation codes

      integer t, i, j

      real x, y, z, s1, beta, delta, xi, k, h, k11, k12, k13, &
      k21, k22, k23, k31, k32, k33, k41, k42, k43, x1, y1, z1, tol, &
      x0, y0, z0, alpha, miu , c0, c1, c2, c3, c4, c5, c6, c7, fp1, &
      fp2,fp3,p,eta,s2,s3, sig1, sig2


       open(100,file='bifurcation.dat')


	xi = 0.25
        p = 0.65
        !miu=0.5
        beta = 1.1
        eta=0.85
        sig1=0.75
        sig2=0.25



	do j = 0, 100
	   miu = 0.0+(1.0/100.0)*j

        write(*,*) j



	x = 0.35
	y = 0.35



	h = 0.01



	do 200 t = 1, 2000000


	k = x*((((1-sig1)*x)-(sig1*y)+sig1)*(1-x-y)-xi)+miu*(y-x)
	k11 = k*h
	k = y*((((beta-sig2)*x)+((2*p*eta-eta-sig2)*y)+sig2)*(1-x-y)-xi)+miu*(x-y)
	k12 = k*h


     	x1 = x + k11/2.0
     	y1 = y + k12/2.0




    	k = x*((((1-sig1)*x)-(sig1*y)+sig1)*(1-x-y)-xi)+miu*(y-x)
	k21 = k*h
	k = y*((((beta-sig2)*x)+((2*p*eta-eta-sig2)*y)+sig2)*(1-x-y)-xi)+miu*(x-y)
	k22 = k*h


     	x1 = x + k21/2.0
     	y1 = y + k22/2.0


    	k = x*((((1-sig1)*x)-(sig1*y)+sig1)*(1-x-y)-xi)+miu*(y-x)
	k31 = k*h
	k = y*((((beta-sig2)*x)+((2*p*eta-eta-sig2)*y)+sig2)*(1-x-y)-xi)+miu*(x-y)
	k32 = k*h

     	x1 = x + k31
     	y1 = y + k32



   	k = x*((((1-sig1)*x)-(sig1*y)+sig1)*(1-x-y)-xi)+miu*(y-x)
	k41 = k*h
	k = y*((((beta-sig2)*x)+((2*p*eta-eta-sig2)*y)+sig2)*(1-x-y)-xi)+miu*(x-y)
	k42 = k*h


     	x = x + k11/6.0 + k21/3.0 + k31/3.0 + k41/6.0
     	y = y + k12/6.0 + k22/3.0 + k32/3.0 + k42/6.0



       if(t.eq.2000000) then
       print*,'******************************************'
       write(*,*) miu, x, y, x+y
       print*,'******************************************'
       write(100,*) miu, x, y, x+y



       end if

200   continue

      end do





      stop
      end
