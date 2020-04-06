Subroutine  compute_macro_var(f, rho0, rho, p, u, nodetype,ly, lx)
!computation of macroscopic variables from distribution function 
!for D2Q9 lattice
	Implicit None
	Real(kind = 8), Intent(In) :: f(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(Inout) :: rho(ly,lx)	!density
	Real(kind = 8), Intent(Inout) :: p(ly,lx)	!pressure
	Real(kind = 8), Intent(Inout) :: u(2,ly,lx)	!velocity
	Real(kind = 8), Intent(In) :: nodetype(ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: rho0	!density of fluid
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: i,j
	Real(kind= 8) :: M1x,M1y
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,M1x,M1y)
	!$OMP DO SCHEDULE(static)
	Do j= 1, lx
	  Do i= 1, ly
	    If (nodetype(i,j) <= 0) Then
	      rho(i, j)= sum(f(i,j,:)) 
	      M1x= f(i,j,2) + f(i,j,4) + f(i,j,5) - f(i,j,6) - f(i,j,8) - f(i,j,9)
	      M1y= f(i,j,3) + f(i,j,4) - f(i,j,5) - f(i,j,7) - f(i,j,8) + f(i,j,9)
	      u(1,i, j)= (1.d0/rho(i,j))*M1x
	      u(2,i, j)= (1.d0/rho(i,j))*M1y
	      p(i, j)=  1.0d0/3.0d0*(rho(i,j) - rho0)
	    Else
	      p(i,j)=0.d0; u(1,i,j)=0.d0; u(2,i,j)=0.d0; rho(i,j)=0.d0
	    End If
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine compute_macro_var
!
Subroutine  get_feq(feq, nodetype, rho, u, ly, lx)
!subroutine to compute the equilibrium distribution function
!for D2Q9 lattice with order 2 terms of u
	Implicit None
	Real(kind = 8), Intent(Out) :: feq(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(InOut) :: nodetype(ly,lx)	!nodetype
	Real(kind = 8), Intent(In) :: rho(ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(In) :: u(2,ly,lx)	!velocity 
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: i, j
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
	!$OMP DO SCHEDULE(static)
	Do j= 1, lx
	  Do i= 1, ly
	   If (nodetype(i,j) <= 0) Then
	    feq(i, j, 1)= rho(i, j)*(-2.0d0/3.0d0*u(1,i, j)**2 - 2.0d0/3.0d0*u(2,i, j)**2 + 4.0d0/9.0d0)
	    feq(i, j, 2)= rho(i, j)*((1.0d0/3.0d0)*u(1,i, j)**2 + (1.0d0/3.0d0)*u(1,i, j) - 1.0d0/6.0d0*u(2,i, j)**2 + &
     	    1.0d0/9.0d0)
	    feq(i, j, 3)= rho(i, j)*(-1.0d0/6.0d0*u(1,i, j)**2 + (1.0d0/3.0d0)*u(2,i, j)**2 + (1.0d0/3.0d0)*u(2,i, j) + &
            1.0d0/9.0d0)
	    feq(i, j, 4)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 + (1.0d0/12.0d0)*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 + ( &
            1.0d0/12.0d0)*u(2,i, j) + (1.0d0/8.0d0)*(u(1,i, j) + u(2,i, j))**2 + 1.0d0/36.0d0)
	    feq(i, j, 5)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 + (1.0d0/12.0d0)*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 - &
            1.0d0/12.0d0*u(2,i, j) + (1.0d0/8.0d0)*(u(1,i, j) - u(2,i, j))**2 + 1.0d0/36.0d0)
	    feq(i, j, 6)= rho(i, j)*((1.0d0/3.0d0)*u(1,i, j)**2 - 1.0d0/3.0d0*u(1,i, j) - 1.0d0/6.0d0*u(2,i, j)**2 + 1.0d0 &
            /9.0d0)
	    feq(i, j, 7)= rho(i, j)*(-1.0d0/6.0d0*u(1,i, j)**2 + (1.0d0/3.0d0)*u(2,i, j)**2 - 1.0d0/3.0d0*u(2,i, j) + &
            1.0d0/9.0d0)
	    feq(i, j, 8)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 - 1.0d0/12.0d0*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 - &
            1.0d0/12.0d0*u(2,i, j) + (1.0d0/8.0d0)*(-u(1,i, j) - u(2,i, j))**2 + 1.0d0/36.0d0)
	    feq(i, j, 9)= rho(i, j)*(-1.0d0/24.0d0*u(1,i, j)**2 - 1.0d0/12.0d0*u(1,i, j) - 1.0d0/24.0d0*u(2,i, j)**2 + ( &
            1.0d0/12.0d0)*u(2,i, j) + (1.0d0/8.0d0)*(-u(1,i, j) + u(2,i, j))**2 + 1.0d0/36.0d0)
          End If
	 End Do
       End Do
       !$OMP END DO
       !$OMP END PARALLEL
End Subroutine get_feq
!
Subroutine  collide_srt(f, rho, u, Fv, nodetype, tau, ly, lx)
!subroutine to compute collision step for single relaxation time scheme
!for D2Q9 lattice with order 2 terms of u considered in feq
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(In) :: rho(ly,lx)	!zeroth order moment of f
	Real(kind = 8), Intent(InOut) :: u(2,ly,lx)	!velocity in y direction
	Real(kind = 8), Intent(In) :: Fv(2,ly,lx)	!forcing
	Real(kind = 8), Intent(In) :: nodetype(ly,lx) !nodetype
	Real(kind = 8), Intent(In) :: tau(ly,lx)	!tau
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: i, j, n
	Real(kind = 8):: omega, omega1, feq(9), SS(9), fswap
     SS=0.0d0
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, n, feq, fswap, SS,omega,omega1)
	!$OMP DO SCHEDULE(static)
	Do j= 1, lx
	  Do i= 1, ly
	    If (nodetype(i,j) <= 0) Then
		  omega = 1.0d0/tau(i,j)
		  omega1 = 1.0d0 - omega
	      feq(1)= rho(i,j)*(-2.0d0/3.0d0*u(1,i,j)**2 - 2.0d0/3.0d0*u(2,i,j)**2 + 4.0d0/9.0d0)
	      feq(2)= rho(i,j)*((1.0d0/3.0d0)*u(1,i,j)**2 + (1.0d0/3.0d0)*u(1,i,j) - 1.0d0/6.0d0*u(2,i,j)**2 + &
              1.0d0/9.0d0)
	      feq(3)= rho(i,j)*(-1.0d0/6.0d0*u(1,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j) + &
              1.0d0/9.0d0)
	      feq(4)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 + (1.0d0/12.0d0)*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 + ( &
              1.0d0/12.0d0)*u(2,i,j) + (1.0d0/8.0d0)*(u(1,i,j) + u(2,i,j))**2 + 1.0d0/36.0d0)
	      feq(5)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 + (1.0d0/12.0d0)*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 - &
              1.0d0/12.0d0*u(2,i,j) + (1.0d0/8.0d0)*(u(1,i,j) - u(2,i,j))**2 + 1.0d0/36.0d0)
	      feq(6)= rho(i,j)*((1.0d0/3.0d0)*u(1,i,j)**2 - 1.0d0/3.0d0*u(1,i,j) - 1.0d0/6.0d0*u(2,i,j)**2 + 1.0d0 &
              /9.0d0)
	      feq(7)= rho(i,j)*(-1.0d0/6.0d0*u(1,i,j)**2 + (1.0d0/3.0d0)*u(2,i,j)**2 - 1.0d0/3.0d0*u(2,i,j) + &
              1.0d0/9.0d0)
	      feq(8)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 - 1.0d0/12.0d0*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 - &
              1.0d0/12.0d0*u(2,i,j) + (1.0d0/8.0d0)*(-u(1,i,j) - u(2,i,j))**2 + 1.0d0/36.0d0)
	      feq(9)= rho(i,j)*(-1.0d0/24.0d0*u(1,i,j)**2 - 1.0d0/12.0d0*u(1,i,j) - 1.0d0/24.0d0*u(2,i,j)**2 + ( &
              1.0d0/12.0d0)*u(2,i,j) + (1.0d0/8.0d0)*(-u(1,i,j) + u(2,i,j))**2 + 1.0d0/36.0d0)
           Call get_forcing(SS,Fv(1,i,j),Fv(2,i,j))
	      f(i,j,:)= omega1*f(i,j,:) + omega*feq + SS
           u(1,i,j) = u(1,i,j) + Fv(1,i,j)/(2.0d0*rho(i,j))
           u(2,i,j) = u(2,i,j) + Fv(2,i,j)/(2.0d0*rho(i,j))
	      Do n=2,5
	        fswap = f(i,j,n)
	        f(i,j,n)=f(i,j,n+4)
	        f(i,j,n+4)=fswap
	      End Do
	    End If
	  End Do
	End Do
	!$OMP END DO
	!$OMP END PARALLEL
     Contains
       Subroutine get_forcing(SS,Fx,Fy)
         Implicit None
         Real(kind = 8), Intent(InOut):: SS(9)
         Real(kind = 8), Intent(In):: Fx, Fy
         SS(1) = 0
         SS(2) = (1.0d0/3.0d0)*Fx
         SS(3) = (1.0d0/3.0d0)*Fy
         SS(4) = (1.0d0/12.0d0)*Fx + (1.0d0/12.0d0)*Fy
         SS(5) = (1.0d0/12.0d0)*Fx - 1.0d0/12.0d0*Fy
         SS(6) = -1.0d0/3.0d0*Fx
         SS(7) = -1.0d0/3.0d0*Fy
         SS(8) = -1.0d0/12.0d0*Fx - 1.0d0/12.0d0*Fy
         SS(9) = -1.0d0/12.0d0*Fx + (1.0d0/12.0d0)*Fy
       End Subroutine
End Subroutine collide_srt
!
Subroutine  stream_and_bounce(f, nodetype, ly, lx)
	Implicit None
	Real(kind = 8), Intent(Inout) :: f(ly,lx,9)	!Distribution function
	Real(kind = 8), Intent(In) :: nodetype(ly,lx)	!nodetype
	Integer, Intent(In) :: ly,lx	!length in y and x direction respectively
	Integer :: ex(9)	!lattice velocity in x direction
	Integer :: ey(9)	!lattice velocity in y direction
	Integer :: i,j, n, half, nextX, nextY
	Real(kind = 8) :: ftemp
	ex = (/0, 1, 0, 1, 1, -1, 0, -1, -1/)
	ey = (/0, 0, 1, 1, -1, 0, -1, -1, 1/)
	half= 4
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, n, nextX, nextY, ftemp)
	!$OMP DO SCHEDULE(static)
	Do j= 1, lx
	  Do i= 1, ly
	    If (nodetype(i,j) <= 0) Then
              Do  n= 2,5
                nextX = j + ex (n)
                nextY = i - ey (n)
                !Apply periodicity
                If (nextX < 1) nextX = lx
                If (nextY < 1) nextY = ly
                If (nextX > lx) nextX = 1
                If (nextY > ly) nextY = 1
                If (nextX > 0 .And. nextX < lx+1 &
&                   .And. nextY > 0  .And. nextY < ly+1) Then
                  If (nodetype(nextY, nextX) <= 0 .And. &
&                     nodetype(nextY,j)<=0 .And. nodetype(i,nextX)<=0) Then
                    ftemp = f (nextY, nextX, n)
	            f(nextY, nextX, n) = f (i, j, n+half)
                    f(i, j, n+half) = ftemp
                  End If
                End If
              End Do
            End If
          End Do
        End Do
	!$OMP END DO
	!$OMP END PARALLEL
End Subroutine stream_and_bounce
!
Subroutine apply_bc(f,rho,rho0,u,nodetype,interp,left_type,left_val, &
&     right_type,right_val,top_type,top_val,bottom_type,bottom_val,ly,lx)
    !Subroutine to apply boundary condition for navier stokes equation
    !for d2q9 lattice. Pressure boundary implemented as extrapolated scheme and 
    !velocity boundary implemented as anti-bounceback scheme
    Use omp_lib
    Implicit none
    Real(kind = 8), Intent(InOut):: f(ly,lx,9) ! distribution function
    Real(kind = 8), Intent(In):: rho(ly,lx) ! density matrix
    Real(kind = 8), Intent(In):: rho0 ! density of fluid
    Real(kind = 8), Intent(In):: u(2,ly,lx) ! velocity
    Real(kind = 8), Intent(In):: nodetype(ly,lx) ! nodetype indicator
    Character(Len=*), Intent(In):: left_type ! left boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: right_type ! right boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Character(Len=*), Intent(In):: top_type ! top boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    character(Len=*), Intent(In):: bottom_type ! bottom boundary type e.g., 'pressure','velocity','normal_vel_grdient','copy_neighbour'
    Real(kind = 8), Intent(In):: left_val(2) ! value for left boundary
    Real(kind = 8), Intent(In):: right_val(2) ! value for right boundary
    Real(kind = 8), Intent(In):: top_val(2) ! value for top boundary
    Real(kind = 8), Intent(In):: bottom_val(2) ! value for bottom boundary
    Integer, Intent(In):: interp ! 1 to allow interpolation to set boundary condtion 0 to not interpolate at the boundary node
    Integer, Intent(In):: ly ! number of nodes in y direction
    Integer, Intent(In):: lx ! number of nodes in x direction
    !other variables...
    Integer:: i 
    Real(kind = 8):: bc_val1,bc_val2,es2,ubx,uby,l
    es2= 0.333333333333
    !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(i,bc_val1,bc_val2,ubx,uby,l)    
    !left bc
    select case (left_type)
     case ('p')
       !$OMP DO SCHEDULE (static)
        Do i= 1,ly
          If (nodetype(i,1)<=0) Then
            !convert pressure to density
            bc_val1 = left_val(1)/es2 + rho0
            If (interp==1) Then
              bc_val1 = 2.0 *bc_val1 - sum(f(i,2,:))
            End If 
            call left_pressure_bc(f(i,1,:),f(i,2,:),bc_val1, &
&             rho(i,2),u(1,i,2),u(2,i,2))
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO SCHEDULE (static)
        Do i= 1, ly
          If (nodetype(i,1)<=0) Then
            bc_val1 = left_val(1)
            bc_val2 = left_val(2)
            !If (interp==1) Then
            !  bc_val1 = 2.0*bc_val1 - u(1,i,2)
            !  bc_val2 = 2.0*bc_val2 - u(2,i,2)
            !End If
            call left_velocity_bc(f(i,1,:),rho0,bc_val1,bc_val2)
          End If
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = left_val(1)
        bc_val2 = left_val(2)
        if ((nint(ly*bc_val2)-ly*bc_val2)<0) Then
          l =nint(ly*bc_val2) + 1
        else
          l=nint(ly*bc_val2)
        endif
        !$OMP DO SCHEDULE (static)
        Do i= 1, ly
          If (nodetype(i,1)<=0) Then
            If (interp==1) Then
              ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
              uby= 0.0d0
            Else
              ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
              uby= 0.0d0
            End If
            call left_velocity_bc(f(i,1,:),rho0,ubx,uby)
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do i= 1, ly
          If (nodetype(i,1)<=0) Then
            f(i,1,:) = f(i,2,:)
          End If
        End Do
        !$OMP END DO
    End select
    select case (right_type)
     case ('p')
       !$OMP DO SCHEDULE (static)
        Do i= 1,ly
          If (nodetype(i,lx)<=0) Then
            !convert pressure to density
            bc_val1 = right_val(1)/es2 + rho0
            If (interp==1) Then
              bc_val1 = 2.0 *bc_val1 - sum(f(i,lx-1,:))
            End If 
            call right_pressure_bc(f(i,lx,:),f(i,lx-1,:),bc_val1, &
&             rho(i,lx-1),u(1,i,lx-1),u(2,i,lx-1))
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO SCHEDULE (static)
        Do i= 1, ly
          If (nodetype(i,lx)<=0) Then
            bc_val1 = right_val(1)
            bc_val2 = right_val(2)
            !If (interp==1) Then
            !  bc_val1 = 2.0*bc_val1 - u(1,i,lx-1)
            !  bc_val2 = 2.0*bc_val2 - u(2,i,lx-1)
            !End If
            call right_velocity_bc(f(i,lx,:),rho0,bc_val1,bc_val2)
          End If
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = right_val(1)
        bc_val2 = right_val(2)
        if ((nint(ly*bc_val2)-ly*bc_val2)<0) Then
          l =nint(ly*bc_val2) + 1
        else
          l=nint(ly*bc_val2)
        endif
        !$OMP DO SCHEDULE (static)
        Do i= 1, ly
          If (nodetype(i,lx)<=0) Then
            If (interp==1) Then
              ubx= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
              uby= 0.0d0
            Else
              ubx= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
              uby= 0.0d0
            End If
            call right_velocity_bc(f(i,lx,:),rho0,ubx,uby)
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do i= 1, ly
          If (nodetype(i,lx)<=0) Then
            f(i,lx,:) = f(i,lx-1,:)
          End If
        End Do
        !$OMP END DO
    End select
    select case (top_type)
     case ('p')
       !$OMP DO SCHEDULE (static)
        Do i= 1,lx
          If (nodetype(1,i)<=0) Then
            !convert pressure to density
            bc_val1 = top_val(1)/es2 + rho0
            If (interp==1) Then
              bc_val1 = 2.0 *bc_val1 - sum(f(2,i,:))
            End If 
            call top_pressure_bc(f(1,i,:),f(2,i,:),bc_val1, &
&             rho(2,i),u(1,2,i),u(2,2,i))
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO SCHEDULE (static)
        Do i= 1, lx
          If (nodetype(1,i)<=0) Then
            bc_val1 = top_val(1)
            bc_val2 = top_val(2)
            !If (interp==1) Then
            !  bc_val1 = 2.0*bc_val1 - u(1,2,i)
            !  bc_val2 = 2.0*bc_val2 - u(2,2,i)
            !End If
            call top_velocity_bc(f(1,i,:),rho0,bc_val1,bc_val2)
          End If
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = top_val(1)
        bc_val2 = top_val(2)
        if ((nint(lx*bc_val2)-lx*bc_val2)<0) Then
          l =nint(lx*bc_val2) + 1
        else
          l=nint(lx*bc_val2)
        endif
        !$OMP DO SCHEDULE (static)
        Do i= 1, lx
          If (nodetype(1,i)<=0) Then
            If (interp==1) Then
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
            Else
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
            End If
            call top_velocity_bc(f(1,i,:),rho0,ubx,uby)
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do i= 1, lx
          If (nodetype(1,i)<=0) Then
            f(1,i,:) = f(2,i,:)
          End If
        End Do
        !$OMP END DO
    End select
    select case (bottom_type)
     case ('p')
       !$OMP DO SCHEDULE (static)
        Do i= 1,lx
          If (nodetype(ly,i)<=0) Then
            !convert pressure to density
            bc_val1 = bottom_val(1)/es2 + rho0
            If (interp==1) Then
              bc_val1 = 2.0 *bc_val1 - sum(f(ly-1,i,:))
            End If 
            call bottom_pressure_bc(f(ly,i,:),f(ly-1,i,:),bc_val1, &
&             rho(ly-1,i),u(1,ly-1,i),u(2,ly-1,i))
          End If
        End Do
        !$OMP END DO
      case ('u')
        !$OMP DO SCHEDULE (static)
        Do i= 1, lx
          If (nodetype(ly,i)<=0) Then
            bc_val1 = bottom_val(1)
            bc_val2 = bottom_val(2)
            !If (interp==1) Then
            !  bc_val1 = 2.0*bc_val1 - u(1,ly-1,i)
            !  bc_val2 = 2.0*bc_val2 - u(2,ly-1,i)
            !End If
            call bottom_velocity_bc(f(ly,i,:),rho0,bc_val1,bc_val2)
          End If
        End Do
        !$OMP END DO
      case ('laminar_flow')
        bc_val1 = bottom_val(1)
        bc_val2 = bottom_val(2)
        if ((nint(lx*bc_val2)-lx*bc_val2)<0) Then
          l =nint(lx*bc_val2) + 1
        else
          l=nint(lx*bc_val2)
        endif
        !$OMP DO SCHEDULE (static)
        Do i= 1, lx
          If (nodetype(ly,i)<=0) Then
            If (interp==1) Then
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-2.0d0)**2*(i-1.5d0)*(l-i-0.5d0)
            Else
              ubx= 0.0d0
              uby= (4.0d0*bc_val1)/(l-1.0d0)**2*(i-1.0d0)*(l-i)
            End If
            call bottom_velocity_bc(f(ly,i,:),rho0,ubx,uby)
          End If
        End Do
        !$OMP END DO
      case ('copy_neighbour')
        !$OMP DO SCHEDULE (static)
        Do i= 1, lx
          If (nodetype(ly,i)<=0) Then
            f(ly,i,:) = f(ly-1,i,:)
          End If
        End Do
        !$OMP END DO
    End select
    !$OMP END PARALLEL
    Contains
    Subroutine left_pressure_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh)
      Real(kind = 8), Intent(In):: fnbh(9),m0,m0nbh,uxnbh,uynbh
      Real(kind = 8), Intent(InOut):: f(9)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*((1.0d0/3.0d0)*uxnbh**2 + (1.0d0/3.0d0)*uxnbh - 1.0d0/6.0d0*uynbh**2 + 1.0d0/ &
      9.0d0)
      feqnbh= m0nbh*((1.0d0/3.0d0)*uxnbh**2 + (1.0d0/3.0d0)*uxnbh - 1.0d0/6.0d0*uynbh**2 + 1.0d0/ &
      9.0d0)
      f(2)= feq + (fnbh(2) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh + (1.0d0/8.0d0)*(uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh + (1.0d0/8.0d0)*(uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      f(4)= feq + (fnbh(4) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      f(5)= feq + (fnbh(5) - feqnbh)
    End Subroutine left_pressure_bc
    Subroutine right_pressure_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh)
      Real(kind = 8), Intent(In):: fnbh(9),m0,m0nbh,uxnbh,uynbh
      Real(kind = 8), Intent(InOut):: f(9)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*((1.0d0/3.0d0)*uxnbh**2 - 1.0d0/3.0d0*uxnbh - 1.0d0/6.0d0*uynbh**2 + 1.0d0/ &
      9.0d0)
      feqnbh= m0nbh*((1.0d0/3.0d0)*uxnbh**2 - 1.0d0/3.0d0*uxnbh - 1.0d0/6.0d0*uynbh**2 + 1.0d0/ &
      9.0d0)
      f(6)= feq + (fnbh(6) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      f(8)= feq + (fnbh(8) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      f(9)= feq + (fnbh(9) - feqnbh)
    End Subroutine right_pressure_bc
    Subroutine top_pressure_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh)
      Real(kind = 8), Intent(In):: fnbh(9),m0,m0nbh,uxnbh,uynbh
      Real(kind = 8), Intent(InOut):: f(9)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0 &
      /12.0d0*uynbh + (1.0d0/8.0d0)*(uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      f(5)= feq + (fnbh(5) - feqnbh)
      feq= m0*(-1.0d0/6.0d0*uxnbh**2 + (1.0d0/3.0d0)*uynbh**2 - 1.0d0/3.0d0*uynbh + 1.0d0/ &
      9.0d0)
      feqnbh= m0nbh*(-1.0d0/6.0d0*uxnbh**2 + (1.0d0/3.0d0)*uynbh**2 - 1.0d0/3.0d0*uynbh + 1.0d0/ &
      9.0d0)
      f(7)= feq + (fnbh(7) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 - 1.0d0/ &
      12.0d0*uynbh + (1.0d0/8.0d0)*(-uxnbh - uynbh)**2 + 1.0d0/36.0d0)
      f(8)= feq + (fnbh(8) - feqnbh)
    End Subroutine top_pressure_bc
    Subroutine bottom_pressure_bc(f,fnbh,m0,m0nbh,uxnbh,uynbh)
      Real(kind = 8), Intent(In):: fnbh(9),m0,m0nbh,uxnbh,uynbh
      Real(kind = 8), Intent(InOut):: f(9)
      Real(kind = 8) :: feq,feqnbh
      feq= m0*(-1.0d0/6.0d0*uxnbh**2 + (1.0d0/3.0d0)*uynbh**2 + (1.0d0/3.0d0)*uynbh + 1.0d0/ &
      9.0d0)
      feqnbh= m0nbh*(-1.0d0/6.0d0*uxnbh**2 + (1.0d0/3.0d0)*uynbh**2 + (1.0d0/3.0d0)*uynbh + 1.0d0/ &
      9.0d0)
      f(3)= feq + (fnbh(3) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh + (1.0d0/8.0d0)*(uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 + (1.0d0/12.0d0)*uxnbh - 1.0d0/24.0d0*uynbh**2 + ( &
      1.0d0/12.0d0)*uynbh + (1.0d0/8.0d0)*(uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      f(4)= feq + (fnbh(4) - feqnbh)
      feq= m0*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      feqnbh= m0nbh*(-1.0d0/24.0d0*uxnbh**2 - 1.0d0/12.0d0*uxnbh - 1.0d0/24.0d0*uynbh**2 + (1.0d0/ &
      12.0d0)*uynbh + (1.0d0/8.0d0)*(-uxnbh + uynbh)**2 + 1.0d0/36.0d0)
      f(9)= feq + (fnbh(9) - feqnbh)
    End Subroutine bottom_pressure_bc
    Subroutine left_velocity_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(2)= f(6)-(-2.0d0/3.0d0*m0*ux)
      f(4)= f(8)-((1.0d0/6.0d0)*m0*(-ux - uy))
      f(5)= f(9)-((1.0d0/6.0d0)*m0*(-ux + uy))
    End Subroutine left_velocity_bc
    Subroutine right_velocity_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(6)= f(2)-((2.0d0/3.0d0)*m0*ux)
      f(8)= f(4)-((1.0d0/6.0d0)*m0*(ux + uy))
      f(9)= f(5)-((1.0d0/6.0d0)*m0*(ux - uy))
    End Subroutine right_velocity_bc
    Subroutine top_velocity_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(5)= f(9)-((1.0d0/6.0d0)*m0*(-ux + uy))
      f(7)= f(3)-((2.0d0/3.0d0)*m0*uy)
      f(8)= f(4)-((1.0d0/6.0d0)*m0*(ux + uy))
    End Subroutine top_velocity_bc
    Subroutine bottom_velocity_bc(f,m0,ux,uy)
      Real(kind = 8), Intent(In):: m0,ux,uy
      Real(kind = 8), Intent(InOut):: f(9)
      f(3)= f(7)-(-2.0d0/3.0d0*m0*uy)
      f(4)= f(8)-((1.0d0/6.0d0)*m0*(-ux - uy))
      f(9)= f(5)-((1.0d0/6.0d0)*m0*(ux - uy))
    End Subroutine bottom_velocity_bc
End Subroutine apply_bc   

