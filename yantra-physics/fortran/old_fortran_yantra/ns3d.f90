Subroutine  compute_macro_var(f, rho0, rho, p, u, nodetype, lz, ly, lx)
!computation of macroscopic variables from distribution function 
!for D3Q19 lattice
	Implicit None
	Real(kind = 8), Intent(In) :: f(lz,ly,lx,19)	!Distribution function
	Real(kind = 8), Intent(Inout) :: rho(lz,ly,lx)	!density
	Real(kind = 8), Intent(Inout) :: p(lz,ly,lx)	!pressure
	Real(kind = 8), Intent(Inout) :: u(3,lz,ly,lx)	!velocity
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: rho0	!density of fluid
	Integer, Intent(In) :: lz, ly, lx	!length in z, y and x direction respectively
	Integer :: i, j, k
	Real(kind = 8) :: M1x, M1y, M1z
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,M1x,M1y,M1z)
	!$OMP DO SCHEDULE(static)
	Do k= 1, lx
	  Do j= 1, ly
	    Do i= 1, lz
	      If (nodetype(i,j,k) <= 0) Then
	        M1x= 0.d0; M1y= 0.d0; M1z= 0.d0
	        rho(i, j, k)= sum(f(i,j,k,:))
                M1x= f(i,j,k,2) + f(i,j,k,5) + f(i,j,k,6) +f(i,j,k,7) + f(i,j,k,8)-f(i,j,k,11) &
&                -f(i,j,k,14) - f(i,j,k,15) - f(i,j,k,16) - f(i,j,k,17)  
                M1y= f(i,j,k,3) + f(i,j,k,7) - f(i,j,k,8) + f(i,j,k,9) - f(i,j,k,10) - &
&                 f(i,j,k,12) - f(i,j,k,16) + f(i,j,k,17) - f(i,j,k,18) + f(i,j,k,19)   
                M1z= f(i,j,k,4) + f(i,j,k,5) - f(i,j,k,6) + f(i,j,k,9) + f(i,j,k,10) - &
&                 f(i,j,k,13) - f(i,j,k,14) + f(i,j,k,15) - f(i,j,k,18) - f(i,j,k,19)   
	        u(1,i, j, k)= (1.d0/rho(i,j,k))*M1x
	        u(2,i, j, k)= (1.d0/rho(i,j,k))*M1y
	        u(3,i, j, k)= (1.d0/rho(i,j,k))*M1z
	        p(i, j, k)= 1/3d0*(rho(i,j,k) - rho0)
	      Else
	        rho(i,j,k)=0.d0; p(i, j, k)=0.d0; u(1,i, j, k)=0.d0; u(2,i, j, k)=0.d0; u(3,i, j, k)=0.d0
	      End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine  get_feq(feq, nodetype, rho, u, lz, ly, lx)
!subroutine to compute the equilibrium distribution function
!for D3Q19 lattice with order 2 terms of u
	Implicit None
	Real(kind = 8), Intent(out) :: feq(lz,ly,lx,19) !Distribution function
	Real(kind = 8), Intent(In) :: rho(lz,ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)!nodetype
	Real(kind = 8), Intent(In) :: u(3,lz,ly,lx)	!velocity 
	Integer, Intent(In) :: lz,ly,lx	!length in z,y, and x direction respectively
	Integer :: i, j, k
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j, k)
	!$OMP DO SCHEDULE(static)
	Do k= 1, lx
	  Do j= 1, ly
	    Do i= 1, lz
           If (nodetype(i,j,k) <= 0) Then
	      feq(i, j, k, 1)= rho(i, j, k)*(-1.0/2.0*u(1,i, j, k)**2 - 1.0/2.0*u(2,i, j, k)**2 - 1.0/2.0*u(3,i, j, k)**2 + &
      1.0/3.0)
	      feq(i, j, k, 2)= rho(i, j, k)*((1.0/6.0)*u(1,i, j, k)**2 + (1.0/6.0)*u(1,i, j, k) - 1.0/12.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 3)= rho(i, j, k)*(-1.0/12.0*u(1,i, j, k)**2 + (1.0/6.0)*u(2,i, j, k)**2 + (1.0/6.0)*u(2,i, j, k) - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 4)= rho(i, j, k)*(-1.0/12.0*u(1,i, j, k)**2 - 1.0/12.0*u(2,i, j, k)**2 + (1.0/6.0)*u(3,i, j, k)**2 + &
      (1.0/6.0)*u(3,i, j, k) + 1.0/18.0)
	      feq(i, j, k, 5)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(u(1,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 6)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(u(1,i, j, k) - u(3,i, j, k) &
      )**2 + 1.0/36.0)
	      feq(i, j, k, 7)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 + ( &
      1.0/12.0)*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(u(1,i, j, k) + &
      u(2,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 8)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 + (1.0/12.0)*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(u(1,i, j, k) - u(2,i, j, k) &
      )**2 + 1.0/36.0)
	      feq(i, j, k, 9)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 + (1.0/12.0)*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(u(2,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 10)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 - 1.0/12.0*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(-u(2,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 11)= rho(i, j, k)*((1.0/6.0)*u(1,i, j, k)**2 - 1.0/6.0*u(1,i, j, k) - 1.0/12.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 12)= rho(i, j, k)*(-1.0/12.0*u(1,i, j, k)**2 + (1.0/6.0)*u(2,i, j, k)**2 - 1.0/6.0*u(2,i, j, k) - &
      1.0/12.0*u(3,i, j, k)**2 + 1.0/18.0)
	      feq(i, j, k, 13)= rho(i, j, k)*(-1.0/12.0*u(1,i, j, k)**2 - 1.0/12.0*u(2,i, j, k)**2 + (1.0/6.0)*u(3,i, j, k)**2 - &
      1.0/6.0*u(3,i, j, k) + 1.0/18.0)
	      feq(i, j, k, 14)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(-u(1,i, j, k) - &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 15)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/24.0*u(3,i, j, k)**2 + (1.0/12.0)*u(3,i, j, k) + (1.0/8.0)*(-u(1,i, j, k) + &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 16)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 - &
      1.0/12.0*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(-u(1,i, j, k) - &
      u(2,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 17)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/12.0*u(1,i, j, k) - 1.0/24.0*u(2,i, j, k)**2 + ( &
      1.0/12.0)*u(2,i, j, k) - 1.0/24.0*u(3,i, j, k)**2 + (1.0/8.0)*(-u(1,i, j, k) + &
      u(2,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 18)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 - 1.0/12.0*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(-u(2,i, j, k) - &
      u(3,i, j, k))**2 + 1.0/36.0)
	      feq(i, j, k, 19)= rho(i, j, k)*(-1.0/24.0*u(1,i, j, k)**2 - 1.0/24.0*u(2,i, j, k)**2 + (1.0/12.0)*u(2,i, j, k) - &
      1.0/24.0*u(3,i, j, k)**2 - 1.0/12.0*u(3,i, j, k) + (1.0/8.0)*(u(2,i, j, k) - u(3,i, j, k) &
      )**2 + 1.0/36.0)
          End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine get_feq
!
Subroutine  collide_srt(f, rho, u, Fv, nodetype, tau, lz, ly, lx)
!subroutine to compute collision step for single relaxation time scheme
!for D3Q19 lattice with order 2 terms of u considered in feq
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(lz, ly,lx,19)	!Distribution function
	Real(kind = 8), Intent(In) :: rho(lz,ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(InOut) :: u(3,lz,ly,lx)	!velocity in z direction
     Real(kind = 8), Intent(In) :: Fv(3,lz,ly,lx) !Volumetric Forcing 
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: tau(lz,ly,lx)	!nodetype
	Integer, Intent(In) :: lz,ly,lx	!length in z, y and x direction respectively
	Integer :: i, j, k, n
	Real(kind = 8):: omega, omega1, feq(19), fswap,ss(19)
    ss = 0.d0
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, n, feq, fswap, ss, omega, omega1)
	!$OMP DO SCHEDULE(static)
	Do k= 1, lx
	   Do j= 1, ly
	    Do i= 1, lz
	      If (nodetype(i,j,k) <= 0) Then
	        omega = 1.d0/tau(i,j,k)
        	omega1 = 1.d0 - omega
  		    feq(1)=  rho(i,j,k)*(-1.0/2.0*u(1,i,j,k)**2 - 1.0/2.0*u(2,i,j,k)**2 - 1.0/2.0*u(3,i,j,k)**2 + &
      1.0/3.0)
	        feq(2)=  rho(i,j,k)*((1.0/6.0)*u(1,i,j,k)**2 + (1.0/6.0)*u(1,i,j,k) - 1.0/12.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(3)=  rho(i,j,k)*(-1.0/12.0*u(1,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k) - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(4)=  rho(i,j,k)*(-1.0/12.0*u(1,i,j,k)**2 - 1.0/12.0*u(2,i,j,k)**2 + (1.0/6.0)*u(3,i,j,k)**2 + &
      (1.0/6.0)*u(3,i,j,k) + 1.0/18.0)
	        feq(5)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(u(1,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(6)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(u(1,i,j,k) - u(3,i,j,k) &
      )**2 + 1.0/36.0)
	        feq(7)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 + ( &
      1.0/12.0)*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(u(1,i,j,k) + &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(8)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 + (1.0/12.0)*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(u(1,i,j,k) - u(2,i,j,k) &
      )**2 + 1.0/36.0)
	        feq(9)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 + (1.0/12.0)*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(u(2,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(10)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 - 1.0/12.0*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(-u(2,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(11)=  rho(i,j,k)*((1.0/6.0)*u(1,i,j,k)**2 - 1.0/6.0*u(1,i,j,k) - 1.0/12.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(12)=  rho(i,j,k)*(-1.0/12.0*u(1,i,j,k)**2 + (1.0/6.0)*u(2,i,j,k)**2 - 1.0/6.0*u(2,i,j,k) - &
      1.0/12.0*u(3,i,j,k)**2 + 1.0/18.0)
	        feq(13)=  rho(i,j,k)*(-1.0/12.0*u(1,i,j,k)**2 - 1.0/12.0*u(2,i,j,k)**2 + (1.0/6.0)*u(3,i,j,k)**2 - &
      1.0/6.0*u(3,i,j,k) + 1.0/18.0)
	        feq(14)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(-u(1,i,j,k) - &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(15)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/24.0*u(3,i,j,k)**2 + (1.0/12.0)*u(3,i,j,k) + (1.0/8.0)*(-u(1,i,j,k) + &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(16)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 - &
      1.0/12.0*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(-u(1,i,j,k) - &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(17)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/12.0*u(1,i,j,k) - 1.0/24.0*u(2,i,j,k)**2 + ( &
      1.0/12.0)*u(2,i,j,k) - 1.0/24.0*u(3,i,j,k)**2 + (1.0/8.0)*(-u(1,i,j,k) + &
      u(2,i,j,k))**2 + 1.0/36.0)
	        feq(18)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 - 1.0/12.0*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(-u(2,i,j,k) - &
      u(3,i,j,k))**2 + 1.0/36.0)
	        feq(19)=  rho(i,j,k)*(-1.0/24.0*u(1,i,j,k)**2 - 1.0/24.0*u(2,i,j,k)**2 + (1.0/12.0)*u(2,i,j,k) - &
      1.0/24.0*u(3,i,j,k)**2 - 1.0/12.0*u(3,i,j,k) + (1.0/8.0)*(u(2,i,j,k) - u(3,i,j,k) &
      )**2 + 1.0/36.0)
             call get_forcing(ss,Fv(1,i,j,k),Fv(2,i,j,k),Fv(3,i,j,k))
             f(i,j,k,:)= omega1*f(i,j,k,:) + omega*feq + ss
             u(1,i,j,k) = u(1,i,j,k) + Fv(1,i,j,k)/(2.0*rho(i,j,k))
             u(2,i,j,k) = u(2,i,j,k) + Fv(2,i,j,k)/(2.0*rho(i,j,k))
             u(3,i,j,k) = u(3,i,j,k) + Fv(3,i,j,k)/(2.0*rho(i,j,k))
	        Do n=2,10
	          fswap = f(i,j,k,n)
	          f(i,j,k,n)=f(i,j,k,n+9)
	          f(i,j,k,n+9)=fswap
	        End Do
	      End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
     Contains
        Subroutine get_forcing(SS,Fx,Fy,Fz)
          Real(kind = 8), Intent(InOut):: SS(19)
          Real(kind = 8), Intent(In):: Fx, Fy, Fz
          SS(1) = 0
          SS(2) = (1.0/6.0)*Fx
          SS(3) = (1.0/6.0)*Fy
          SS(4) = (1.0/6.0)*Fz
          SS(5) = (1.0/12.0)*Fx + (1.0/12.0)*Fz
          SS(6) = (1.0/12.0)*Fx - 1.0/12.0*Fz
          SS(7) = (1.0/12.0)*Fx + (1.0/12.0)*Fy
          SS(8) = (1.0/12.0)*Fx - 1.0/12.0*Fy
          SS(9) = (1.0/12.0)*Fy + (1.0/12.0)*Fz
          SS(10) = -1.0/12.0*Fy + (1.0/12.0)*Fz
          SS(11) = -1.0/6.0*Fx
          SS(12) = -1.0/6.0*Fy
          SS(13) = -1.0/6.0*Fz
          SS(14) = -1.0/12.0*Fx - 1.0/12.0*Fz
          SS(15) = -1.0/12.0*Fx + (1.0/12.0)*Fz
          SS(16) = -1.0/12.0*Fx - 1.0/12.0*Fy
          SS(17) = -1.0/12.0*Fx + (1.0/12.0)*Fy
          SS(18) = -1.0/12.0*Fy - 1.0/12.0*Fz
          SS(19) = (1.0/12.0)*Fy - 1.0/12.0*Fz
        End Subroutine
End Subroutine collide_srt
!
Subroutine  stream_and_bounce(f, nodetype, lz, ly, lx)
!subroutine for propogation step using latt(2007) swap algorithim
!for D3Q19 lattices
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(lz,ly,lx,19)	!Distribution function
	Real(kind = 8), Intent(In) :: nodetype(lz,ly,lx)	!nodetype
	Integer, Intent(In) :: lz, ly, lx	!length in y and x direction respectively
	Integer :: ex(19)	!lattice velocity in x direction
	Integer :: ey(19)	!lattice velocity in y direction
	Integer :: ez(19)	!lattice velocity in z direction
	Integer :: i, j, k, n, half, nextX, nextY, nextZ
	Real(kind = 8) :: ftemp
	ex = (/0, 1, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, -1, -1, -1, -1, 0, 0/)
	ey = (/0, 0, 1, 0, 0, 0, 1, -1, 1, -1, 0, -1, 0, 0, 0, -1, 1, -1, 1/)
	ez = (/0, 0, 0, 1, 1, -1, 0, 0, 1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1/)
	half= 9
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, n, nextX, nextY, nextZ, ftemp)
	!$OMP DO SCHEDULE(static)
	Do k= 1, lx
	  Do j= 1, ly
	    Do i= 1, lz
	      If (nodetype(i,j,k) <= 0) Then
	        Do  n= 2,10
	          nextX = k + ex (n)
	          nextY = j - ey (n)
	          nextz = i - ez (n)
	          !Apply periodicity
                  If (nextX < 1) nextX = lx
                  If (nextY < 1) nextY = ly
                  If (nextZ < 1) nextZ = lz
                  If (nextX > lx) nextX = 1
                  If (nextY > ly) nextY = 1
                  If (nextZ > lz) nextZ = 1
	          If (nextX > 0 .And. nextX < lx+1 &
	          & .And. nextY > 0  .And. nextY < ly+1 &
	          & .And. nextZ > 0  .And. nextY < lz+1) Then
	            If (nodetype(nextZ, nextY, nextX) <= 0 .And. &
	            & nodetype(nextZ,j,k)<=0 .And. nodetype(i,nextY,k)<=0 .And. nodetype(i,j,nextX) <=0) Then
	            !If (nodetype(nextZ, nextY, nextX) <= 0) Then
	              ftemp = f (nextZ, nextY, nextX, n)
	              f(nextZ, nextY, nextX, n) = f (i, j, k, n+half)
	              f(i, j, k, n+half) = ftemp
	            End If
	          End If
	        End Do
	      End If
	    End Do
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine stream_and_bounce
!
Subroutine apply_bc(f,rho,rho0,u,nodetype,interp,left_type,left_val, &
&     right_type,right_val,top_type,top_val,bottom_type,bottom_val,front_type,front_val, &
&     back_type,back_val,lz,ly,lx)
    !Subroutine to apply boundary condition for navier stokes equation
    !for d3q19 lattice. p boundary implemented as extrapolated scheme and 
    !u boundary implemented as anti-bounceback scheme
    Use omp_lib
    Implicit none
    Real(kind = 8), Intent(InOut):: f(lz,ly,lx,19) ! distribution function
    Real(kind = 8), Intent(In):: rho(lz,ly,lx) ! density matrix
    Real(kind = 8), Intent(In):: rho0 ! density of fluid
    Real(kind = 8), Intent(In):: u(3,lz,ly,lx) ! u
    Real(kind = 8), Intent(In):: nodetype(lz,ly,lx) ! nodetype indicator
    Character(Len=*), Intent(In):: left_type ! left boundary type e.g., 'p','u','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: right_type ! right boundary type e.g., 'p','u','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: top_type ! top boundary type e.g., 'p','u','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: bottom_type ! bottom boundary type e.g., 'p','u','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: front_type ! front boundary type e.g., 'p','u','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: back_type ! back boundary type e.g., 'p','u','normal_vel_grdient','copy_neighbour'
    Real(kind = 8), Intent(In):: left_val(3) ! value for left boundary
    Real(kind = 8), Intent(In):: right_val(3) ! value for right boundary
    Real(kind = 8), Intent(In):: top_val(3) ! value for top boundary
    Real(kind = 8), Intent(In):: bottom_val(3) ! value for bottom boundary
    Real(kind = 8), Intent(In):: front_val(3) ! value for front boundary
    Real(kind = 8), Intent(In):: back_val(3) ! value for back boundary    
    Integer, Intent(In):: interp ! 1 to allow interpolation to set boundary condtion 0 to not interpolate at the boundary node
    Integer, Intent(In):: lz ! number of nodes in z direction
    Integer, Intent(In):: ly ! number of nodes in y direction
    Integer, Intent(In):: lx ! number of nodes in x direction
    !other variables...
    Integer:: i, j, k
    Real(kind = 8):: bc_val1,bc_val2,bc_val3,l,ubx,uby,ubz,es2
    es2= 0.333333333333
    !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(i,j,k,bc_val1,bc_val2,bc_val3,ubx,uby,ubz,l)    
    !left bc
    select case (left_type)
     case ('p')
        bc_val1 = left_val(1)/es2 + rho0
       !$OMP DO SCHEDULE (static)
        Do j= 1,ly
          Do i=1,lz
            If (nodetype(i,j,1)<=0) Then
              !convert p to density
              If (interp==1) Then
                bc_val1 = 2.0 *bc_val1 - sum(f(i,j,2,:))
              End If 
              call left_p_bc(f(i,j,1,:),f(i,j,2,:),bc_val1, &
&               rho(i,j,2),u(1,i,j,2),u(2,i,j,2),u(3,i,j,2))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = left_val(1)
        bc_val2 = left_val(2)
        bc_val3 = left_val(3)
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,1)<=0) Then
              !If (interp==1) Then
              !  bc_val1 = 2.0*bc_val1 - u(1,i,j,2)
              !  bc_val2 = 2.0*bc_val2 - u(2,i,j,2)
              !  bc_val3 = 2.0*bc_val3 - u(3,i,j,2)
              !End If
              call left_u_bc(f(i,j,1,:),rho0,bc_val1,bc_val2,bc_val3)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = left_val(1)
        bc_val2 = left_val(2)
        bc_val3 = left_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(ly*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(ly*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,1)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                uby= 0.0d0
                ubz= 0.0d0
              Else
                ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                uby= 0.0d0
                ubz= 0.0d0
              End If
              call left_u_bc(f(i,j,1,:),rho0,ubx,uby,ubz)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i= 1,lz
            If (nodetype(i,j,1)<=0) Then
              f(i,j,1,:) = f(i,j,2,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (right_type)
     case ('p')
        bc_val1 = right_val(1)/es2 + rho0
       !$OMP DO SCHEDULE (static)
        Do j= 1,ly
          Do i=1,lz
            If (nodetype(i,j,lx)<=0) Then
              !convert p to density
              If (interp==1) Then
                bc_val1 = 2.0 *bc_val1 - sum(f(i,j,lx-1,:))
              End If 
              call right_p_bc(f(i,j,lx,:),f(i,j,lx-1,:),bc_val1, &
&               rho(i,j,lx-1),u(1,i,j,lx-1),u(2,i,j,lx-1),u(3,i,j,lx-1))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = right_val(1)
        bc_val2 = right_val(2)
        bc_val3 = right_val(3)
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,lx)<=0) Then
              !If (interp==1) Then
              !  bc_val1 = 2.0*bc_val1 - u(1,i,j,lx-1)
              !  bc_val2 = 2.0*bc_val2 - u(2,i,j,lx-1)
              !  bc_val3 = 2.0*bc_val3 - u(3,i,j,lx-1)
              !End If
              call right_u_bc(f(i,j,lx,:),rho0,bc_val1,bc_val2,bc_val3)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = right_val(1)
        bc_val2 = right_val(2)
        bc_val3 = right_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(ly*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(ly*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, lz
            If (nodetype(i,j,lx)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                uby= 0.0d0
                ubz= 0.0d0
              Else
                ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                uby= 0.0d0
                ubz= 0.0d0
              End If
              call right_u_bc(f(i,j,lx,:),rho0,ubx,uby,ubz)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i= 1,lz
            If (nodetype(i,j,lx)<=0) Then
              f(i,j,lx,:) = f(i,j,lx-1,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (top_type)
     case ('p')
        bc_val1 = top_val(1)/es2 + rho0
       !$OMP DO SCHEDULE (static)
        Do j= 1,lx
          Do i=1,lz
            If (nodetype(i,1,j)<=0) Then
              !convert p to density
              If (interp==1) Then
                bc_val1 = 2.0 *bc_val1 - sum(f(i,2,j,:))
              End If 
              call top_p_bc(f(i,1,j,:),f(i,2,j,:),bc_val1, &
&               rho(i,2,j),u(1,i,2,j),u(2,i,2,j),u(3,i,2,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = top_val(1)
        bc_val2 = top_val(2)
        bc_val3 = top_val(3)
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,1,j)<=0) Then
              !If (interp==1) Then
              !  bc_val1 = 2.0*bc_val1 - u(1,i,2,j)
              !  bc_val2 = 2.0*bc_val2 - u(2,i,2,j)
              !  bc_val3 = 2.0*bc_val3 - u(3,i,2,j)
              !End If
              call top_u_bc(f(i,1,j,:),rho0,bc_val1,bc_val2,bc_val3)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = top_val(1)
        bc_val2 = top_val(2)
        bc_val3 = top_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(lx*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(lx*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,1,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                ubz= 0.0d0
              Else
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                ubz= 0.0d0
              End If
              call top_u_bc(f(i,1,j,:),rho0,ubx,uby,ubz)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i= 1,lz
            If (nodetype(i,1,j)<=0) Then
              f(i,1,j,:) = f(i,2,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (bottom_type)
     case ('p')
        bc_val1 = bottom_val(1)/es2 + rho0
       !$OMP DO SCHEDULE (static)
        Do j= 1,lx
          Do i=1,lz
            If (nodetype(i,ly,j)<=0) Then
              !convert p to density
              If (interp==1) Then
                bc_val1 = 2.0 *bc_val1 - sum(f(i,ly-1,j,:))
              End If 
              call bottom_p_bc(f(i,ly,j,:),f(i,ly-1,j,:),bc_val1, &
&               rho(i,ly-1,j),u(1,i,ly-1,j),u(2,i,ly-1,j),u(3,i,ly-1,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = bottom_val(1)
        bc_val2 = bottom_val(2)
        bc_val3 = bottom_val(3)
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,ly,j)<=0) Then
              !If (interp==1) Then
              !  bc_val1 = 2.0*bc_val1 - u(1,i,ly-1,j)
              !  bc_val2 = 2.0*bc_val2 - u(2,i,ly-1,j)
              !  bc_val3 = 2.0*bc_val3 - u(3,i,ly-1,j)
              !End If
              call bottom_u_bc(f(i,ly,j,:),rho0,bc_val1,bc_val2,bc_val3)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = bottom_val(1)
        bc_val2 = bottom_val(2)
        bc_val3 = bottom_val(3)
        if ((nint(lz*bc_val2)-lz*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lz*bc_val2) + 1
          else
            l =nint(lx*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lz*bc_val2)
          else
            l=nint(lx*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lz
            If (nodetype(i,ly,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
                ubz= 0.0d0
              Else
                ubx= 0.0d0
                uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
                ubz= 0.0d0
              End If
              call bottom_u_bc(f(i,ly,j,:),rho0,ubx,uby,ubz)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i= 1,lz
            If (nodetype(i,ly,j)<=0) Then
              f(i,ly,j,:) = f(i,ly-1,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (front_type)
     case ('p')
        bc_val1 = front_val(1)/es2 + rho0
       !$OMP DO SCHEDULE (static)
        Do j= 1,ly
          Do i=1,ly
            If (nodetype(1,i,j)<=0) Then
              !convert p to density
              If (interp==1) Then
                bc_val1 = 2.0 *bc_val1 - sum(f(2,i,j,:))
              End If 
              call front_p_bc(f(1,i,j,:),f(2,i,j,:),bc_val1, &
&               rho(2,i,j),u(1,2,i,j),u(2,2,i,j),u(3,2,i,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = front_val(1)
        bc_val2 = front_val(2)
        bc_val3 = front_val(3)
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, ly
            If (nodetype(1,i,j)<=0) Then
              !If (interp==1) Then
              !  bc_val1 = 2.0*bc_val1 - u(1,2,i,j)
              !  bc_val2 = 2.0*bc_val2 - u(2,2,i,j)
              !  bc_val3 = 2.0*bc_val3 - u(3,2,i,j)
              !End If
              call front_u_bc(f(1,i,j,:),rho0,bc_val1,bc_val2,bc_val3)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = front_val(1)
        bc_val2 = front_val(2)
        bc_val3 = front_val(3)
        if ((nint(ly*bc_val2)-ly*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(ly*bc_val2) + 1
          else
            l =nint(ly*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(ly*bc_val2)
          else
            l=nint(ly*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i=1, ly
            If (nodetype(1,i,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
              Else
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
              End If
              call front_u_bc(f(1,i,j,:),rho0,ubx,uby,ubz)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, ly
          Do i= 1,ly
            If (nodetype(1,i,j)<=0) Then
              f(1,i,j,:) = f(2,i,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    select case (back_type)
     case ('p')
        bc_val1 = back_val(1)/es2 + rho0
       !$OMP DO SCHEDULE (static)
        Do j= 1,lx
          Do i=1,lx
            If (nodetype(lz,i,j)<=0) Then
              !convert p to density
              If (interp==1) Then
                bc_val1 = 2.0 *bc_val1 - sum(f(lz-1,i,j,:))
              End If 
              call back_p_bc(f(lz,i,j,:),f(lz-1,i,j,:),bc_val1, &
&               rho(lz-1,i,j),u(1,lz-1,i,j),u(2,lz-1,i,j),u(3,lz-1,i,j))
            End If
          End Do
        End Do
        !$OMP END DO
      case ('u')
        bc_val1 = back_val(1)
        bc_val2 = back_val(2)
        bc_val3 = back_val(3)
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lx
            If (nodetype(lz,i,j)<=0) Then
              !If (interp==1) Then
              !  bc_val1 = 2.0*bc_val1 - u(1,lz-1,i,j)
              !  bc_val2 = 2.0*bc_val2 - u(2,lz-1,i,j)
              !  bc_val3 = 2.0*bc_val3 - u(3,lz-1,i,j)
              !End If
              call back_u_bc(f(lz,i,j,:),rho0,bc_val1,bc_val2,bc_val3)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = back_val(1)
        bc_val2 = back_val(2)
        bc_val3 = back_val(3)
        if ((nint(lx*bc_val2)-lx*bc_val2)<0) Then
          if (bc_val3 == 1) Then
            l =nint(lx*bc_val2) + 1
          else
            l =nint(lx*bc_val2) + 1
          end if
        else
          if (bc_val3 == 1) Then 
            l=nint(lx*bc_val2)
          else
            l=nint(lx*bc_val2)
          end if
        endif
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i=1, lx
            If (nodetype(lz,i,j)<=0) Then
              If (interp==1) Then
                if (bc_val3==1) Then
                  k=i
                else
                  k=j
                end if
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-2.0d0)**2*(k-1.5d0)*(l-k-0.5d0)
              Else
                ubx= 0.0d0
                uby= 0.0d0
                ubz= (4.0d0*bc_val1)/(l-1.0d0)**2*(k-1.0d0)*(l-k)
              End If
              call back_u_bc(f(lz,i,j,:),rho0,ubx,uby,ubz)
            End If
          End Do
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do j= 1, lx
          Do i= 1,lx
            If (nodetype(lz,i,j)<=0) Then
              f(lz,i,j,:) = f(lz-1,i,j,:)
            End If
          End Do
        End Do
        !$OMP END DO
    End select
    !$OMP END PARALLEL
    Contains
    Subroutine left_p_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh,uznbh)
      Real(kind = 8), Intent(In):: fnbh(19),m0,m0nbh,uxnbh,uynbh,uznbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*((1.0d0/6.0d0)*uxnbh**2 + (1.0d0/6.0d0)*uxnbh - 1.0d0/12.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uznbh**2 + 1.0d0/18.0d0)
      feqnbh= m0nbh*((1.0d0/6.0d0)*uxnbh**2 + (1.0d0/6.0d0)*uxnbh - 1.0d0/12.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uznbh**2 + 1.0d0/18.0d0)
      f(2)= feq + (fnbh(2) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(5)= feq + (fnbh(5) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(6)= feq + (fnbh(6) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh + uynbh)** &
      2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh + uynbh)** &
      2 + 1.0d0/36.0d0)
      f(7)= feq + (fnbh(7) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      f(8)= feq + (fnbh(8) - feqnbh)
    End Subroutine left_p_bc
    Subroutine right_p_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh,uznbh)
      Real(kind = 8), Intent(In):: fnbh(19),m0,m0nbh,uxnbh,uynbh,uznbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*((1.0d0/6.0d0)*uxnbh**2 - 1.0d0/6.0d0*uxnbh - 1.0d0/12.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uznbh**2 + 1.0d0/18.0d0)
      feqnbh= m0nbh*((1.0d0/6.0d0)*uxnbh**2 - 1.0d0/6.0d0*uxnbh - 1.0d0/12.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uznbh**2 + 1.0d0/18.0d0)
      f(11)= feq + (fnbh(11) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(14)= feq + (fnbh(14) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(15)= feq + (fnbh(15) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      f(16)= feq + (fnbh(16) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + &
      1.0d0/36.0d0)
      f(17)= feq + (fnbh(17) - feqnbh)
    End Subroutine right_p_bc
    Subroutine top_p_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh,uznbh)
      Real(kind = 8), Intent(In):: fnbh(19),m0,m0nbh,uxnbh,uynbh,uznbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      f(8)= feq + (fnbh(8) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(10)= feq + (fnbh(10) - feqnbh)
      feq= m0*(-1.0d0/12.0d0*uxnbh**2 + (1.0d0/6.0d0)*uynbh**2 - 1.0d0/6.0d0*uynbh - 1.0d0/ &
      12.0d0*uznbh**2 + 1.0d0/18.0d0)
      feqnbh= m0nbh*(-1.0d0/12.0d0*uxnbh**2 + (1.0d0/6.0d0)*uynbh**2 - 1.0d0/6.0d0*uynbh - 1.0d0/ &
      12.0d0*uznbh**2 + 1.0d0/18.0d0)
      f(12)= feq + (fnbh(12) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + &
      1.0d0/36.0d0)
      f(16)= feq + (fnbh(16) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(18)= feq + (fnbh(18) - feqnbh)
    End Subroutine top_p_bc
    Subroutine bottom_p_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh,uznbh)
      Real(kind = 8), Intent(In):: fnbh(19),m0,m0nbh,uxnbh,uynbh,uznbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*(-1.0d0/12.0d0*uxnbh**2 + (1.0d0/6.0d0)*uynbh**2 + (1.0d0/6.0d0)*uynbh - 1.0d0 &
      /12.0d0*uznbh**2 + 1.0d0/18.0d0)
      feqnbh= m0nbh*(-1.0d0/12.0d0*uxnbh**2 + (1.0d0/6.0d0)*uynbh**2 + (1.0d0/6.0d0)*uynbh - 1.0d0 &
      /12.0d0*uznbh**2 + 1.0d0/18.0d0)
      f(3)= feq + (fnbh(3) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh + uynbh)** &
      2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(uxnbh + uynbh)** &
      2 + 1.0d0/36.0d0)
      f(7)= feq + (fnbh(7) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(9)= feq + (fnbh(9) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh - 1.0d0/24.0d0*uznbh**2 + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + &
      1.0d0/36.0d0)
      f(17)= feq + (fnbh(17) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(19)= feq + (fnbh(19) - feqnbh)
    End Subroutine bottom_p_bc
    Subroutine front_p_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh,uznbh)
      Real(kind = 8), Intent(In):: fnbh(19),m0,m0nbh,uxnbh,uynbh,uznbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(6)= feq + (fnbh(6) - feqnbh)
      feq= m0*(-1.0d0/12.0d0*uxnbh**2 - 1.0d0/12.0d0*uynbh**2 + (1.0d0/6.0d0)*uznbh**2 - &
      1.0d0/6.0d0*uznbh + 1.0d0/18.0d0)
      feqnbh= m0nbh*(-1.0d0/12.0d0*uxnbh**2 - 1.0d0/12.0d0*uynbh**2 + (1.0d0/6.0d0)*uznbh**2 - &
      1.0d0/6.0d0*uznbh + 1.0d0/18.0d0)
      f(13)= feq + (fnbh(13) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uxnbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(14)= feq + (fnbh(14) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(-uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(18)= feq + (fnbh(18) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 - 1.0d0/12.0d0*uznbh + (1.0d0/8.0d0)*(uynbh - uznbh)**2 + &
      1.0d0/36.0d0)
      f(19)= feq + (fnbh(19) - feqnbh)
    End Subroutine front_p_bc
    Subroutine back_p_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh,uznbh)
      Real(kind = 8), Intent(In):: fnbh(19),m0,m0nbh,uxnbh,uynbh,uznbh
      Real(kind = 8), Intent(InOut):: f(19)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*(-1.0d0/12.0d0*uxnbh**2 - 1.0d0/12.0d0*uynbh**2 + (1.0d0/6.0d0)*uznbh**2 + ( &
      1.0d0/6.0d0)*uznbh + 1.0d0/18.0d0)
      feqnbh= m0nbh*(-1.0d0/12.0d0*uxnbh**2 - 1.0d0/12.0d0*uynbh**2 + (1.0d0/6.0d0)*uznbh**2 + ( &
      1.0d0/6.0d0)*uznbh + 1.0d0/18.0d0)
      f(4)= feq + (fnbh(4) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(5)= feq + (fnbh(5) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 + (1.0d0/12.0d0)*uynbh - 1.0d0 &
      /24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(9)= feq + (fnbh(9) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/24.0d0*uynbh**2 - 1.0d0/12.0d0*uynbh - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uynbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(10)= feq + (fnbh(10) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      24.0d0*uznbh**2 + (1.0d0/12.0d0)*uznbh + (1.0d0/8.0d0)*(-uxnbh + uznbh)**2 + &
      1.0d0/36.0d0)
      f(15)= feq + (fnbh(15) - feqnbh)
    End Subroutine back_p_bc
    Subroutine left_u_bc(f,m0,ux,uy,uz)
      Real(kind = 8), Intent(In):: m0,ux,uy,uz
      Real(kind = 8), Intent(InOut):: f(19)
      f(2)= f(11)-(-1.0d0/3.0d0*m0*ux)
      f(5)= f(14)-((1.0d0/6.0d0)*m0*(-ux - uz))
      f(6)= f(15)-((1.0d0/6.0d0)*m0*(-ux + uz))
      f(7)= f(16)-((1.0d0/6.0d0)*m0*(-ux - uy))
      f(8)= f(17)-((1.0d0/6.0d0)*m0*(-ux + uy))
    End Subroutine left_u_bc
    Subroutine right_u_bc(f,m0,ux,uy,uz)
      Real(kind = 8), Intent(In):: m0,ux,uy,uz
      Real(kind = 8), Intent(InOut):: f(19)
      f(11)= f(2)-((1.0d0/3.0d0)*m0*ux)
      f(14)= f(5)-((1.0d0/6.0d0)*m0*(ux + uz))
      f(15)= f(6)-((1.0d0/6.0d0)*m0*(ux - uz))
      f(16)= f(7)-((1.0d0/6.0d0)*m0*(ux + uy))
      f(17)= f(8)-((1.0d0/6.0d0)*m0*(ux - uy))
    End Subroutine right_u_bc
    Subroutine top_u_bc(f,m0,ux,uy,uz)
      Real(kind = 8), Intent(In):: m0,ux,uy,uz
      Real(kind = 8), Intent(InOut):: f(19)
      f(8)= f(17)-((1.0d0/6.0d0)*m0*(-ux + uy))
      f(10)= f(19)-((1.0d0/6.0d0)*m0*(uy - uz))
      f(12)= f(3)-((1.0d0/3.0d0)*m0*uy)
      f(16)= f(7)-((1.0d0/6.0d0)*m0*(ux + uy))
      f(18)= f(9)-((1.0d0/6.0d0)*m0*(uy + uz))
    End Subroutine top_u_bc
    Subroutine bottom_u_bc(f,m0,ux,uy,uz)
      Real(kind = 8), Intent(In):: m0,ux,uy,uz
      Real(kind = 8), Intent(InOut):: f(19)
      f(3)= f(12)-(-1.0d0/3.0d0*m0*uy)
      f(7)= f(16)-((1.0d0/6.0d0)*m0*(-ux - uy))
      f(9)= f(18)-((1.0d0/6.0d0)*m0*(-uy - uz))
      f(17)= f(8)-((1.0d0/6.0d0)*m0*(ux - uy))
      f(19)= f(10)-((1.0d0/6.0d0)*m0*(-uy + uz))
    End Subroutine bottom_u_bc
    Subroutine front_u_bc(f,m0,ux,uy,uz)
      Real(kind = 8), Intent(In):: m0,ux,uy,uz
      Real(kind = 8), Intent(InOut):: f(19)
      f(6)= f(15)-((1.0d0/6.0d0)*m0*(-ux + uz))
      f(13)= f(4)-((1.0d0/3.0d0)*m0*uz)
      f(14)= f(5)-((1.0d0/6.0d0)*m0*(ux + uz))
      f(18)= f(9)-((1.0d0/6.0d0)*m0*(uy + uz))
      f(19)= f(10)-((1.0d0/6.0d0)*m0*(-uy + uz))
    End Subroutine front_u_bc
    Subroutine back_u_bc(f,m0,ux,uy,uz)
      Real(kind = 8), Intent(In):: m0,ux,uy,uz
      Real(kind = 8), Intent(InOut):: f(19)
      f(4)= f(13)-(-1.0d0/3.0d0*m0*uz)
      f(5)= f(14)-((1.0d0/6.0d0)*m0*(-ux - uz))
      f(9)= f(18)-((1.0d0/6.0d0)*m0*(-uy - uz))
      f(10)= f(19)-((1.0d0/6.0d0)*m0*(uy - uz))
      f(15)= f(6)-((1.0d0/6.0d0)*m0*(ux - uz))
    End Subroutine back_u_bc
End subroutine apply_bc
