!=======================================================================================
!This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
!multiphyics simulations
!=======================================================================================
!
!Copyright (C) 2016-2017  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
!
!This program is free software: you can redistribute it and/or modify it under the
!terms of the GNU General Public License as published by the Free Software 
!Foundation, either version 3 of the License, or any later version.
!This program is distributed in the hope that it will be useful, but WITHOUT ANY 
!WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
!PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!Alogrithims to help with geometry update for 2D model
!
!=======================================================================================
!
Subroutine reassign_nodetype(nodetype, ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
      Implicit none
      Real (8), Intent (InOut) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: ex(4),ey(4),i,j,k,next_i,next_j,flag
      ex=(/1,0,-1,0/)
      ey=(/0,1,0,-1/)
!      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
!      !$OMP DO
      Do i = 1, ly
        Do j = 1, lx
          flag = 0
          check: Do k = 1, 4
            next_i = i+ey(k)
            next_j = j+ex(k)
            If (next_i < 1) next_i = 1
            If (next_j < 1) next_j =1
            If (next_i > ly) next_i = ly
            If (next_j > lx) next_j = lx
            If ((nodetype(next_i,next_j) > 0) .AND. (nodetype(i,j) <=0)) Then
              flag = 1
              exit check
            End If
            If ((nodetype(next_i,next_j) <=0) .AND. (nodetype(i,j)  >0)) Then
              flag = 1
              exit check
            End if 
          End Do check
          If (nodetype(i,j)>0) Then
            If (flag == 1) Then
              nodetype(i,j)=1
            else
              nodetype(i,j)=2
            End If
          Else
            If (flag == 1 .AND. nodetype(i,j)/=-5) Then
              nodetype(i,j)=0
            else If (nodetype(i,j)/=-5) Then
              nodetype(i,j)=-1
            End If
          End If 
        End Do
      End Do
!      !$OMP END DO
!      !$OMP END PARALLEL
End Subroutine reassign_nodetype
!

Subroutine reassign_nodetype_8direc(nodetype, ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
      Implicit none
      Real (8), Intent (InOut) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx
      Integer :: ex(8),ey(8),i,j,k,next_i,next_j,flag
      ex=(/1,1,0,-1,-1,-1,0,1/)
      ey=(/0,1,1,1,0,-1,-1,-1/)
!      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
!      !$OMP DO
      Do i = 1, ly
        Do j = 1, lx
          flag = 0
          check: Do k = 1, 8
            next_i = i+ey(k)
            next_j = j+ex(k)
            If (next_i < 1) next_i = 1
            If (next_j < 1) next_j =1
            If (next_i > ly) next_i = ly
            If (next_j > lx) next_j = lx
            If ((nodetype(next_i,next_j) > 0) .AND. (nodetype(i,j) <=0)) Then
              flag = 1
              exit check
            End If
            If ((nodetype(next_i,next_j) <=0) .AND. (nodetype(i,j)  >0)) Then
              flag = 1
              exit check
            End if 
          End Do check
          If (nodetype(i,j)>0) Then
            If (flag == 1) Then
              nodetype(i,j)=1
            else
              nodetype(i,j)=2
            End If
          Else
            If (flag == 1 .AND. nodetype(i,j)/=-5) Then
              nodetype(i,j)=0
            else If (nodetype(i,j)/=-5) Then
              nodetype(i,j)=-1
            End If
          End If 
        End Do
      End Do
!      !$OMP END DO
!      !$OMP END PARALLEL
End Subroutine reassign_nodetype_8direc
!




!
Subroutine redistribute(nodetype,c,f,cphi, ly, lx)
!ditribute the c in newly formed precipited solid node to the neibouring nodes 

      Implicit none
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (Inout) :: c(ly,lx)
      Real (8), Intent (Inout) :: f(ly,lx,5)
      Integer, Intent (In) :: ly, lx
      Integer :: ex(4),ey(4),i,j,k,nextI,nextJ
      REAL (8):: wi,cphi
      ex=(/1,0,-1,0/)
      ey=(/0,1,0,-1/)
!      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k)
!      !$OMP DO
      Do i = 1, ly
        Do j = 1, lx
           If (nodetype(i, j) ==1) Then

              !get weights for neighbours
              wi = 0.d0
              Do k = 1, 4
                nextI = ey (k) + i
                nextJ = ex (k) + j

                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1.and. nodetype(nextI, nextJ) <= 0.and. nodetype(nextI, nextJ) /=-10) Then
                    wi = wi + 1.d0
                End If
              End Do
              If (wi > 0) Then
                wi = 1.d0 / wi
              Else 
                wi = 0.d0
              End If


              Do k = 1, 4
                nextI = ey (k) + i
                nextJ = ex (k) + j
                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1.and. nodetype(nextI, nextJ) <= 0 .and. nodetype(nextI, nextJ) /=-10) Then
                    
                    f(nextI,nextJ,2)=f(nextI,nextJ,2)+(f(I,J,1)+f(I,J,2)+f(I,J,3)+f(I,J,4)+f(I,J,5))*cphi*0.5*wi
                    f(nextI,nextJ,3)=f(nextI,nextJ,3)+(f(I,J,1)+f(I,J,2)+f(I,J,3)+f(I,J,4)+f(I,J,5))*cphi*0.5*wi
                    f(nextI,nextJ,4)=f(nextI,nextJ,4)+(f(I,J,1)+f(I,J,2)+f(I,J,3)+f(I,J,4)+f(I,J,5))*cphi*0.5*wi
                    f(nextI,nextJ,5)=f(nextI,nextJ,5)+(f(I,J,1)+f(I,J,2)+f(I,J,3)+f(I,J,4)+f(I,J,5))*cphi*0.5*wi
                    f(nextI,nextJ,1)=f(nextI,nextJ,1)+(f(I,J,1)+f(I,J,2)+f(I,J,3)+f(I,J,4)+f(I,J,5))*(1-4*cphi*0.5)*wi


                endif
              enddo 

              f(I,J,1)=0.0
              f(I,J,2)=0.0
              f(I,J,3)=0.0
              f(I,J,4)=0.0
              f(I,J,5)=0.0
              c(I,J)=0.0

          End If 
        End Do
      End Do
!      !$OMP END DO
!      !$OMP END PARALLEL
End Subroutine redistribute
!




!
Subroutine reaction_surface(phaseqty,dphase,nodetype,RA,np,ly, lx)
!lessen precipitation at extrdos corner 

      Implicit none
      Real (8), Intent (In) :: phaseqty (np, ly, lx)
      Real (8), Intent (In) :: dphase (np, ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (Inout) :: RA(ly,lx)
      Integer, Intent (In) :: ly, lx,np

      Integer :: ex(4),ey(4),i,j,k,p,nextI,nextJ
      Integer :: n1(4),n10,n20
      REAL (8):: wi
      ex=(/1,0,-1,0/)
      ey=(/0,1,0,-1/)
!      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k)
!      !$OMP DO

      Do i = 1, ly
        Do j = 1, lx
          RA(i, j) =0.0
        Enddo
      Enddo


      Do j = 1, lx
        if(nodetype(1,j)==0)RA(1, j) =1.0
        if(nodetype(ly,j)==0)RA(ly, j) =1.0
      Enddo

      Do i = 1, ly
        if(nodetype(i,1)==0)RA(i, 1) =1.0
        if(nodetype(i,lx)==0)RA(i, lx) =1.0

      Enddo

      Do i = 1, ly
        Do j = 1, lx
           If (nodetype(i, j) ==1.or.nodetype(i, j) ==-10) Then

              !get weights for neighbours
              wi = 0.d0

              n1(k)=0 
              n10=0
              n20=0

              Do k = 1, 4
                nextI = ey (k) + i
                nextJ = ex (k) + j

                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1.and. nodetype(nextI, nextJ) == 0) Then
                    wi = wi + 1.d0
                    n10=n10+1
                    n1(n10)=k

                End If

                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1.and. (nodetype(nextI, nextJ) >= 1.or.nodetype(nextI, nextJ) ==-10)) Then

                    n20=n20+1

                End If
              End Do

              If (n10==1 .And. n20==3)Then
	          RA(ey (n1(1)) + i,ex (n1(1)) + j)=1.0

              End If

              If ((n10== 2 .And. n20==2) .And. ((abs(n1(1)-n1(2))==1) .or. (abs(n1(1)-n1(2))==3)))Then
                  If (RA(ey (n1(1)) + i,ex (n1(1)) + j)<0.5)Then
	              RA(ey (n1(1)) + i,ex (n1(1)) + j)=0.5
                  endif
                  If (RA(ey (n1(2)) + i,ex (n1(2)) + j)<0.5)Then 
                      RA(ey (n1(2)) + i,ex (n1(2)) + j)=0.5
                  endif	
              End If

              If (n10== 3 .And. n20==1)Then
                  If (RA(ey (n1(1)) + i,ex (n1(1)) + j)<1./3.)Then
	              RA(ey (n1(1)) + i,ex (n1(1)) + j)=1./3.
                  endif

                  If (RA(ey (n1(2)) + i,ex (n1(2)) + j)<1./3.)Then 
                      RA(ey (n1(2)) + i,ex (n1(2)) + j)=1./3.
                  endif

                  If (RA(ey (n1(3)) + i,ex (n1(3)) + j)<1./3.)Then 
                      RA(ey (n1(3)) + i,ex (n1(3)) + j)=1./3.
                  endif	
              End If

              If (n10== 4)Then

                  If (RA(ey (n1(1)) + i,ex (n1(1)) + j)<0.25)Then
	              RA(ey (n1(1)) + i,ex (n1(1)) + j)=0.25
                  endif

                  If (RA(ey (n1(2)) + i,ex (n1(2)) + j)<0.25)Then 
                      RA(ey (n1(2)) + i,ex (n1(2)) + j)=0.25
                  endif

                  If (RA(ey (n1(3)) + i,ex (n1(3)) + j)<0.25)Then 
                      RA(ey (n1(3)) + i,ex (n1(3)) + j)=0.25
                  endif
                  If (RA(ey (n1(4)) + i,ex (n1(4)) + j)<0.25)Then 
                      RA(ey (n1(4)) + i,ex (n1(4)) + j)=0.25
                  endif
	
              End If


          End If 
        End Do
      End Do



! some internal solid mass must keeping dissoving
    Do p = 1, np
      Do i = 1, ly
        Do j = 1, lx
           If (nodetype(i, j) ==0 .and. phaseqty(p,i,j)>0.and. dphase(p,i,j)<=0) Then
              if (RA(i,j)==0)RA(i,j)=1.0
           Endif
        End Do
      End Do
    Enddo



!      !$OMP END DO
!      !$OMP END PARALLEL
End Subroutine reaction_surface
!




!
Subroutine reaction_surface2(phaseqty,dphase,nodetype,RA,np,ly, lx)
!lessen precipitation at extrdos corner 

      Implicit none
      Real (8), Intent (In) :: phaseqty (np, ly, lx)
      Real (8), Intent (In) :: dphase (np, ly, lx)
      Real (8), Intent (In) :: nodetype (ly, lx)
      Real (8), Intent (Inout) :: RA(ly,lx)
      Integer, Intent (In) :: ly, lx,np

      Integer :: ex(4),ey(4),i,j,k,p,nextI,nextJ
      Integer :: n1(4),n10,n20
      REAL (8):: wi
      ex=(/1,0,-1,0/)
      ey=(/0,1,0,-1/)
!      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k)
!      !$OMP DO

      Do i = 1, ly
        Do j = 1, lx
          RA(i, j) =0.0
        Enddo
      Enddo


! wall nodes at boundary
      Do j = 1, lx
        if(nodetype(1,j)==0)RA(1, j) =1.0
        if(nodetype(ly,j)==0)RA(ly, j) =1.0
      Enddo

      Do i = 1, ly
        if(nodetype(i,1)==0)RA(i, 1) =1.0
        if(nodetype(i,lx)==0)RA(i, lx) =1.0

      Enddo

      Do i = 1, ly
        Do j = 1, lx
           If (nodetype(i, j) ==0) Then

              !get weights for solid neighbours
              wi = 0.d0

              Do k = 1, 4
                nextI = ey (k) + i
                nextJ = ex (k) + j

                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1.and. (nodetype(nextI, nextJ)> 0 .or. nodetype(nextI, nextJ)==-10)) Then
                    wi = wi + 1.d0

                End If
              End Do
              RA(i,j)=1.0*wi           

          End If 
        End Do
      End Do



! some internal solid mass must keeping dissoving
    Do p = 1, np
      Do i = 1, ly
        Do j = 1, lx
           If (nodetype(i, j) ==0 .and. phaseqty(p,i,j)>0.and. dphase(p,i,j)<=0) Then
              if (RA(i,j)==0)RA(i,j)=1.0
           Endif
        End Do
      End Do
    Enddo



!      !$OMP END DO
!      !$OMP END PARALLEL
End Subroutine reaction_surface2
!




Subroutine update_nodetype(nodetype, vol, totvol, frac, ly, lx)
!update node type to fluid if  vol/totvol < frac and 
!vol/totvol > frac update nodetype to solid
    Implicit None
    Real (8), Intent(Inout) :: nodetype(ly,lx)
    Real (8), Intent(In) :: vol(ly,lx)
    Real (8), Intent(In) :: totvol
    Real (8), Intent(In) :: frac
    Integer, Intent(In) :: ly, lx
    Integer:: i, j        
    Real(8) :: phi
!    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,phi)
!    !$OMP DO
    Do i = 1,ly
      Do j = 1,lx
        phi = vol(i, j)/totvol
        If (nodetype(i, j) > 0) Then
          !dissolve
          If (phi < frac) nodetype(i, j) = 0
        Else If (nodetype(i, j) <= 0) Then
          !precipitate
          If (phi > frac) nodetype(i, j) = 1
        End If
      End Do
    End Do
!    !$OMP END DO
!    !$OMP END PARALLEL
End Subroutine update_nodetype


!
Subroutine reassign_nodetype_poros(phaseqty,RA,nodetype, np,ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
      Implicit none
      Real (8), Intent (In) :: phaseqty (np, ly, lx)
      Real (8), Intent (InOut) :: RA (ly, lx)
      Real (8), Intent (InOut) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx,np
      Integer :: ex(4),ey(4),i,j,p,k,next_i,next_j,flag, flag2
      Real (8) :: wi
      Integer :: ni,nextI,nextJ
      ex=(/1,0,-1,0/)
      ey=(/0,1,0,-1/)
!      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
!      !$OMP DO
  Do p = 1,np
      Do i = 1, ly
        Do j = 1, lx
          flag = 0
          flag2 =0
          check: Do k = 1, 4
            next_i = i+ey(k)
            next_j = j+ex(k)
            If (next_i < 1) next_i = 1
            If (next_j < 1) next_j =1
            If (next_i > ly) next_i = ly
            If (next_j > lx) next_j = lx
            If ((nodetype(next_i,next_j) > 0) .AND. (nodetype(i,j) <=0)) Then
              flag = 1
            End If

            If ((nodetype(next_i,next_j) <=0) .AND. (nodetype(i,j)  >0)) Then
              flag = 1
            End if 

            If ((nodetype(next_i,next_j)==-10) .AND. (nodetype(i,j)<=0)) Then

              flag2 = 1

            End if

            If ((nodetype(i,j)==-10) .AND. (nodetype(next_i,next_j)<=0)) Then

              flag2 = 1

            End if

          End Do check
          If (nodetype(i,j)>0) Then
            If (flag == 1) Then
              nodetype(i,j)=1
            else
              nodetype(i,j)=2
            End If
          Else
            If ((flag == 1 .or. flag2==1) .AND. (nodetype(i,j)/=-5 .AND. nodetype(i,j)/=-10)) Then
              nodetype(i,j)=0
            else If (nodetype(i,j)/=-5 .AND. nodetype(i,j)/=-10) Then
                
              If(phaseqty(p,i,j)<=0)nodetype(i,j)=-1
               
            End If
            
          End If

        End Do
      End Do

  Enddo
!      !$OMP END DO
!      !$OMP END PARALLEL
End Subroutine reassign_nodetype_poros


!
Subroutine reassign_nodetype_poros_8direc(phaseqty,RA,nodetype, np,ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
      Implicit none
      Real (8), Intent (In) :: phaseqty (np, ly, lx)
      Real (8), Intent (InOut) :: RA (ly, lx)
      Real (8), Intent (InOut) :: nodetype (ly, lx)
      Integer, Intent (In) :: ly, lx,np
      Integer :: ex(8),ey(8),i,j,p,k,next_i,next_j,flag, flag2
      Real (8) :: wi
      Integer :: ni,nextI,nextJ
      ex=(/1,1,0,-1,-1,-1,0,1/)
      ey=(/0,1,1,1,0,-1,-1,-1/)
!      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
!      !$OMP DO
  Do p = 1,np
      Do i = 1, ly
        Do j = 1, lx
          flag = 0
          flag2 =0
          check: Do k = 1, 8
            next_i = i+ey(k)
            next_j = j+ex(k)
            If (next_i < 1) next_i = 1
            If (next_j < 1) next_j =1
            If (next_i > ly) next_i = ly
            If (next_j > lx) next_j = lx
            If ((nodetype(next_i,next_j) > 0) .AND. (nodetype(i,j) <=0)) Then
              flag = 1
            End If

            If ((nodetype(next_i,next_j) <=0) .AND. (nodetype(i,j)  >0)) Then
              flag = 1
            End if 

            If ((nodetype(next_i,next_j)==-10) .AND. (nodetype(i,j)<=0)) Then

              flag2 = 1

            End if

            If ((nodetype(i,j)==-10) .AND. (nodetype(next_i,next_j)<=0)) Then

              flag2 = 1

            End if

          End Do check
          If (nodetype(i,j)>0) Then
            If (flag == 1) Then
              nodetype(i,j)=1
            else
              nodetype(i,j)=2
            End If
          Else
            If ((flag == 1 .or. flag2==1) .AND. (nodetype(i,j)/=-5 .AND. nodetype(i,j)/=-10)) Then
              nodetype(i,j)=0
            else If (nodetype(i,j)/=-5 .AND. nodetype(i,j)/=-10) Then
                
              If(phaseqty(p,i,j)<=0)nodetype(i,j)=-1
               
            End If
            
          End If

        End Do
      End Do

  Enddo
!      !$OMP END DO
!      !$OMP END PARALLEL
End Subroutine reassign_nodetype_poros_8direc


Subroutine update_nodetype_poros(phaseqty,dphase,nodetype, vol, totvol, frac, np,ly, lx)
!update node type to fluid if  vol/totvol < frac and 
!vol/totvol > frac update nodetype to solid
    Implicit None
    Real (8), Intent (In) :: phaseqty (np, ly, lx)
    Real (8), Intent(Inout) :: nodetype(ly,lx)
    Real (8), Intent(In) :: dphase(np,ly,lx)
    Real (8), Intent(In) :: vol(ly,lx)
    Real (8), Intent(In) :: totvol
    Real (8), Intent(In) :: frac
    Integer, Intent(In) :: ly, lx,np
    Integer:: i, j ,p 
    Real(8) :: flag(ly, lx)      
    Real(8) :: phi
!    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,phi)
!    !$OMP DO
    flag = 0.0d0
    Do p = 1,np
        Do i = 1,ly
            Do j = 1,lx
                phi = vol(i, j)/totvol
                If (nodetype(i, j) > 0) Then
                    !dissolve
                    If (phi < frac) nodetype(i, j) = 0
                Else If (nodetype(i, j) == -10) Then
                    !dissolve
                    If (phi< frac) nodetype(i, j) = 0
                Else If (nodetype(i, j) == 0) Then
                    !precipitate
                    If (dphase(p,i,j)>0 .And. (phi > frac)) nodetype(i, j) = -10
                    !         the last liquid interal node
                    If (phaseqty(p,i,j)==0 .and. phaseqty(p,i-1,j)==0 .and. phaseqty(p,i+1,j)==0 &
                        & .and. phaseqty(p,i,j-1)==0 .and. phaseqty(p,i,j+1)==0) Then
                            flag(i,j) = flag(i,j) + 1
                            If (flag(i,j) == np) nodetype(i, j) = -1
                    End If
                        !        Else If (nodetype(i, j) == -1) Then
                        !         !precipitate
                        !          pass
                End If
            End Do
        End Do
    End Do
!    !$OMP END DO
!    !$OMP END PARALLEL
End Subroutine update_nodetype_poros
!
Subroutine update_non_diffusive_phases(dphase, phaseqty, nodetype, vol, mvol, &
& totvol, np, ly, lx)
  implicit none
  Real (8), Intent (Inout) :: phaseqty (np, ly, lx)
  Real (8), Intent (Inout) :: vol (ly, lx)
  Real (8), Intent (In) :: dphase (np, ly, lx)
  Real (8), Intent (In) :: nodetype (ly, lx)
  Real (8), Intent (In) :: mvol(np) 
  Real(8), Intent (In) :: totvol
  Integer, Intent(In) :: ly, lx, np
  Real (8) :: wi, rdphase, vreq, vavail, vreq_t
  Integer :: ex (4), ey (4), i, j, p, nextI, nextJ, ni, itr
  Integer :: Is,Js
  ex = (/ 1, 0, -1, 0 /)
  ey = (/ 0, 1, 0, -1 /)
  Do p = 1,np
    Do j = 1, lx
      Do  i = 1, ly
        If (nodetype(i, j)<=0) then
          If (dphase(p, i, j) < 0) Then
            !dissolve
            !get weights for neighbours
            wi = 0.d0
            Do ni = 1, 4
              nextI = ey (ni) + i
              nextJ = ex (ni) + j
              If (nextI > 0 .And. nextI < ly + 1 .And. nextJ > 0 .And. nextJ &
&                 < lx + 1) Then
                If (nodetype(nextI, nextJ) > 0 .And. phaseqty(p, nextI, &
&                   nextJ) > 0) Then
                  wi = wi + 1.d0
                End If
              End If
            End Do
            If (wi > 0) Then
              wi = 1.d0 / wi
            Else 
              wi = 0.d0
            End If
            !update phase and volumes
            If (phaseqty(p, i, j) > 0) Then
              If (phaseqty(p, i, j) >-1.*dphase(p, i, j)) Then
                vol (i, j) = vol (i, j) + dphase (p, i, j) * mvol(p)
                phaseqty (p, i, j) = phaseqty (p, i, j) + dphase (p, i, j)
              Else
                rdphase = dphase (p, i, j) + phaseqty (p, i, j)

!yuli modified from '-' to '+'
                vol (i, j) = vol (i, j) - phaseqty (p, i, j) * mvol(p)
                phaseqty (p, i, j) = phaseqty(p, i, j)  + (-1.d0 * phaseqty(p,i,j)) 
                Do ni = 1, 4
                  nextI = ey (ni) + i
                  nextJ = ex (ni) + j
                  If (nextI > 0 .And. nextI < ly +1 .And. nextJ > 0 .And. &
&                     nextJ < lx+1) Then
                    If (nodetype(nextI, nextJ) > 0 .And. &
&                       phaseqty(p,nextI, nextJ) > 0) Then
                      vol (nextI, nextJ) = vol (nextI, nextJ) + &
&                       rdphase * mvol(p) * wi
                      phaseqty (p, nextI, nextJ) = phaseqty &
&                       (p, nextI, nextJ) + rdphase * wi 
                    End If
                  End If
                End Do
              End If
            Else
              Do ni = 1, 4
                nextI = ey (ni) + i
                nextJ = ex (ni) + j
                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1) Then
                  If (nodetype(nextI, nextJ) > 0 .And. &
&                     phaseqty(p, nextI, nextJ) > 0) Then
                    vol (nextI, nextJ) = vol (nextI, nextJ) + &
&                     dphase (p, i, j) * mvol(p) * wi
                    phaseqty (p, nextI, nextJ) = phaseqty &
&                     (p, nextI, nextJ) + dphase (p, i, j) * wi
                  End If
                End If
              End Do
            End If
          Else If (dphase(p, i, j) > 0) Then
            !precipitate
            vreq = dphase (p, i, j) * mvol(p) !required volume to be filled
            Do itr = 1, 10
              !get weights for neighbours
              wi = 0.d0
              Do ni = 1, 4
                nextI = ey (ni) + i
                nextJ = ex (ni) + j
                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1) Then
                  vavail = totvol - vol (nextI, nextJ)

                  If (nodetype(nextI, nextJ) > 0 .And. vavail > 0) Then
                    wi = wi + 1.d0
                  End If
                End If
              End Do
              If (wi > 0) Then
                wi = 1.d0 / wi
              Else 
                wi = 0.d0
              End If 
    
              vreq_t = wi * vreq                
              !get neighbours filled
              Do ni = 1, 4
                nextI = ey (ni) + i
                nextJ = ex (ni) + j
                If (nextI > 0 .And. nextI < ly+1 .And. nextJ > 0 .And. &
&                   nextJ < lx+1) Then
                  vavail = totvol - vol (nextI, nextJ)
                  If (nodetype(nextI, nextJ) > 0 .And. vavail >  0) Then
                    If (vavail > vreq_t) Then
                      vreq = vreq - vreq_t
                      vol (nextI, nextJ) = vol (nextI, &
&                       nextJ) + vreq_t
                      phaseqty (p, nextI, nextJ) = phaseqty &
&                       (p, nextI, nextJ) + (vreq_t/mvol(p)) 
                    Else
                      vreq = vreq - vavail
                      vol (nextI, nextJ) = vol(nextI, nextJ) + vavail
                      phaseqty (p, nextI, nextJ) = phaseqty &
&                       (p, nextI, nextJ) + (vavail/mvol(p))

 
                    End If
                  End If
                End If
              End Do

              If (vreq <= 1e-30 .Or. wi <= 0) Exit
            End Do

            If (vreq > 0) Then
              vol (i, j) = vol (i, j) + vreq
              phaseqty (p, i, j) = phaseqty (p, i, j) + (vreq/mvol(p)) 
            End If
          End If
        End If
      End Do
    End Do
  End Do
End Subroutine
!
subroutine update_diffusive_phases(dphase,phaseqty,nodetype,vol,mvol,np,ly,lx)
   Implicit None
   real(8),intent(inout):: nodetype(ly,lx)
   real(8),intent(in):: dphase(np,ly,lx)
   real(8),intent(in):: mvol(np)
   real(8),intent(inout):: phaseqty(np,ly,lx)
   real(8),intent(inout):: vol(ly,lx)
   integer,intent(in):: np,ly,lx
   integer:: i, j, p

   Real (8) :: wi, rdphase, vreq, vavail, vreq_t
   Integer :: ex (4), ey (4), nextI, nextJ, ni, itr
   Integer :: Is,Js
   ex = (/ 1, 0, -1, 0 /)
   ey = (/ 0, 1, 0, -1 /)


!   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,p)
!   !$OMP DO REDUCTION(+:vol)
   Do p = 1,np
     Do i = 1,ly
       Do j = 1,lx

         !------- precipitation 
         If (nodetype(i,j) <=0 .AND. dphase(p,i,j) >= 0) then
           phaseqty(p,i,j) = phaseqty(p,i,j) + dphase(p,i,j)
           vol(i,j) = vol(i,j) + dphase(p,i,j) * mvol(p)  
         End If  

         !------- dissolution 
         If (nodetype(i,j) <=0 .AND. dphase(p,i,j) < 0) then

           if (-1.0*dphase(p,i,j)<=phaseqty(p,i,j)) then   !solid mass is enough for dissolution
             phaseqty(p,i,j)=phaseqty(p,i,j)+dphase(p,i,j)
             vol(i,j) = vol(i,j) + dphase(p,i,j) * mvol(p) 

           else !dissove dephase of the neighbouring nodes
             rdphase = dphase (p, i, j) + phaseqty (p, i, j)  !mass needs to be dissoved from neighbouring nodes
             vol (i, j) = vol (i, j) - phaseqty (p, i, j) * mvol(p)
             phaseqty (p,i, j) = phaseqty(p, i, j)  + (-1.d0 * phaseqty(p,i,j)) 




            !get weights for neighbours
             wi = 0.d0
             Do ni = 1, 4
               nextI = ey (ni) + i
               nextJ = ex (ni) + j
               If (nextI > 0 .And. nextI < ly + 1 .And. nextJ > 0 .And. nextJ < lx + 1) Then
                 If (phaseqty(p, nextI, nextJ) > 0) Then
                   wi = wi + 1.d0

                 End If
               End If
             End Do
             If (wi > 0) Then
               wi = 1.d0 / wi
             Else 
               wi = 0.d0
             End If


             Do ni = 1, 4
               nextI = ey (ni) + i
               nextJ = ex (ni) + j
               If (nextI > 0 .And. nextI < ly +1 .And. nextJ > 0 .And. nextJ < lx+1) Then
                 If (phaseqty(p,nextI, nextJ) > 0) Then
                   vol (nextI, nextJ) = vol (nextI, nextJ) + rdphase * mvol(p) * wi
                   phaseqty (p, nextI, nextJ) = phaseqty(p, nextI, nextJ) + rdphase * wi 
                 End If
               End If
             End Do
           End If
         End If

     
       End Do
     End Do
   End Do
!   !$OMP END DO
!   !$OMP END PARALLEL
end subroutine update_diffusive_phases
!
subroutine porosity(poros,vol,totvol,nodetype,ly,lx)
   Implicit None
   Real(8), Intent(Out):: poros(ly,lx)
   Real(8), Intent(In):: vol(ly,lx)
   Real(8), Intent(In):: nodetype(ly,lx)
   Real(8), Intent(In):: totvol
   Integer, Intent(In):: ly,lx
   Integer:: i,j
   poros = 0.d0
!   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j)
!   !$OMP DO
   Do  i = 1,ly
     Do j = 1,lx 
        If (nodetype(i,j) <=0) Then
            poros(i, j) = (1.-vol(i,j))/totvol
        End If
     End Do
   End Do
!   !$OMP END DO
!   !$OMP END PARALLEL
end subroutine porosity
!
subroutine phaseqty_for_phrqc(phaseqty_phrqc, nodetype, phaseqty,np, ly, lx)
         Implicit None
         Real (8), Intent (Out) :: phaseqty_phrqc (np, ly, lx)
         Real (8), Intent (In) :: nodetype (ly, lx)
         Real (8), Intent (In) :: phaseqty (np, ly, lx)
         Integer, Intent (In) :: np,ly, lx
         Real (8):: wi 
         Integer ::ex(4),ey(4), ni,nextI,nextJ,i,j,p
         ex = (/ 1, 0, -1, 0 /)
         ey = (/ 0, 1, 0, -1 /)
         phaseqty_phrqc  = 0.d0
!         !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(wi,ni,nextI,nextJ,i,j,p)
!         !$OMP DO REDUCTION(+:phaseqty_phrqc)
         Do p = 1,np
           Do i = 1, ly
             Do j = 1, lx
               wi = 0.d0
               If (nodetype(i, j) > 0 .AND. phaseqty(p,i,j) > 0) Then
                 !Get weights
                 Do ni = 1, 4
                   nextI = ey (ni) + i
                   nextJ = ex (ni) + j
                   If (nextI > 0 .And. nextI < ly + 1 .And. nextJ > 0 .And. nextJ &
&                      < lx + 1) Then
                     If (nodetype(nextI, nextJ) <= 0) Then
                       wi = wi + 1.d0
                     End If
                   End If
                 End Do
                 If (wi > 0) Then
                   wi = 1.d0/wi 
                 Else 
                   wi = 0.d0
                 End If
                 !Distribute values
                 Do ni = 1, 4
                   nextI = ey (ni) + i
                   nextJ = ex (ni) + j
                   If (nextI > 0 .And. nextI < ly + 1 .And. nextJ > 0 .And. nextJ &
&                      < lx + 1) Then
                     If (nodetype(nextI, nextJ) <= 0) Then
                       phaseqty_phrqc(p,nextI,nextJ) = phaseqty_phrqc(p,nextI,nextJ) + wi * phaseqty(p, i, j)
                     End If
                   End If
                 End Do
               Else If (nodetype(i,j) <=0) Then
                 phaseqty_phrqc(p,i,j)=phaseqty_phrqc(p,i,j)+phaseqty(p,i,j)
               End If
             End Do
           End Do
         End Do
!         !$OMP END DO
!         !$OMP END PARALLEL
End subroutine


! by Anna
!==================================================================
!
Subroutine reassign_mlvl(nodetype, ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
    Implicit none
    Real (8), Intent (InOut) :: nodetype (ly, lx)
    Integer, Intent (In) :: ly, lx    
    Integer :: ex(4), ey(4)
    Integer :: i,j,k
    Integer :: next_i, next_j, flag
   
    ex=(/1,0,-1,0/)
    ey=(/0,1,0,-1/)
    !      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
    !      !$OMP DO
    Do i = 1, ly
        Do j = 1, lx
            flag = 0
            check: Do k = 1, 4
                next_i = i+ey(k)
                next_j = j+ex(k)
                If (next_i < 1) next_i = 1
                If (next_j < 1) next_j =1
                If (next_i > ly) next_i = ly
                If (next_j > lx) next_j = lx
                If ((nodetype(next_i,next_j) == -5) .AND. (nodetype(i,j) == -1 )) Then
                    flag = 1
                    exit check
                End If
            End Do check
            If (nodetype(i,j) == -1 ) Then
                If (flag == 1) Then
                    nodetype(i,j)= - 2
                else
                    nodetype(i,j)= -1
                End If
            End If 
                
        End Do
    End Do
    !      !$OMP END DO
    !      !$OMP END PARALLEL
End Subroutine reassign_mlvl
!

!
Subroutine reassign_mlvl_solid(nodetype, ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
    Implicit none
    Real (8), Intent (InOut) :: nodetype (ly, lx)
    Integer, Intent (In) :: ly, lx    
    Integer :: ex(4), ey(4)
    Integer :: i,j,k
    Integer :: next_i, next_j, flag
   
    ex=(/1,0,-1,0/)
    ey=(/0,1,0,-1/)
    !      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
    !      !$OMP DO
    Do i = 1, ly
        Do j = 1, lx
            flag = 0
            check: Do k = 1, 4
                next_i = i+ey(k)
                next_j = j+ex(k)
                If (next_i < 1) next_i = 1
                If (next_j < 1) next_j =1
                If (next_i > ly) next_i = ly
                If (next_j > lx) next_j = lx
                If ((nodetype(next_i,next_j) == 1) .AND. (nodetype(i,j) == -1 )) Then
                    flag = 1
                    exit check
                End If
            End Do check
            If (nodetype(i,j) == -1 ) Then
                If (flag == 1) Then
                    nodetype(i,j)= - 2
                else
                    nodetype(i,j)= -1
                End If
            End If 
                
        End Do
    End Do
    !      !$OMP END DO
    !      !$OMP END PARALLEL
End Subroutine reassign_mlvl_solid
!

!
Subroutine reassign_mlvl_8dir(nodetype, vol, totvol, frac, ly, lx)
!marks fluid interface node as 0 and solid interface node
!as 1
    Implicit none
    Real (8), Intent (InOut) :: nodetype (ly, lx)
    Real(8),intent(inout):: vol(ly,lx)
    Real (8), Intent(In) :: totvol
    Real (8), Intent(In) :: frac
    Integer, Intent (In) :: ly, lx
    
    Real (8) ::phi, phi_next
    Integer :: ex(8), ey(8)
    Integer :: i,j,k
    Integer :: next_i, next_j, flag
   
    ex=(/1,1,0,-1,-1,-1,0,1/)
    ey=(/0,1,1,1,0,-1,-1,-1/)
    !      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,next_i,next_j,k,flag)
    !      !$OMP DO
    Do i = 1, ly
        Do j = 1, lx
            phi = vol(i,j)/ totvol
            flag = 0
            check: Do k = 1, 8
                next_i = i+ey(k)
                next_j = j+ex(k)
                If (next_i < 1) next_i = 1
                If (next_j < 1) next_j =1
                If (next_i > ly) next_i = ly
                If (next_j > lx) next_j = lx
                phi_next = vol(next_i,next_j)/ totvol
                If ((nodetype(next_i,next_j) == -5) .And. (phi_next >= frac) .AND. (nodetype(i,j) == -1 )) Then
                    flag = 1
                    exit check
                End If
            End Do check
            If ((nodetype(i,j) == -5) .And.  (phi <= 0) ) Then
                nodetype(i,j)= -1
            End If
            If (nodetype(i,j) == -1 ) Then
                If (flag == 1) Then
                    nodetype(i,j)= - 5
                else
                    nodetype(i,j)= -1
                End If
            End If 
                
        End Do
    End Do
    !      !$OMP END DO
    !      !$OMP END PARALLEL
End Subroutine reassign_mlvl_8dir
!

subroutine update_mlvl_phases(dphase,phaseqty,nodetype,vol,mvol,totvol,frac,np,ly,lx)
   Implicit None
   real(8),intent(inout):: nodetype(ly,lx)
   real(8),intent(in):: dphase(np,ly,lx)
   real(8),intent(in):: mvol(np)
   real(8),intent(inout):: phaseqty(np,ly,lx)
   real(8),intent(inout):: vol(ly,lx)
   Real (8), Intent(In) :: totvol
   Real (8), Intent(In) :: frac
   integer,intent(in):: np,ly,lx
   integer:: i, j, p

   Real (8) ::phi
   Integer :: ex (4), ey (4), nextI, nextJ, ni, itr
   Integer :: Is,Js
   ex = (/ 1, 0, -1, 0 /)
   ey = (/ 0, 1, 0, -1 /)


!   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,p)
!   !$OMP DO REDUCTION(+:vol)
   Do p = 1,np
     Do i = 1,ly
       Do j = 1,lx    
          phi = vol(i,j)/ totvol     
         !------- precipitation 
         !If (nodetype(i,j) == -5 ) then
         If (p == 2) then
             If(dphase(p,i,j) >= 0) then
                 If (phi < frac) then
                     phaseqty(p,i,j) = phaseqty(p,i,j) + dphase(p,i,j)
                     vol(i,j) = vol(i,j) + dphase(p,i,j) * mvol(p)
                 End If
             End If
         End If
         
         !------- dissolution 
         If (nodetype(i,j) == -5) then
             If(dphase(p,i,j) <= 0) then
                 phaseqty(p,i,j) = phaseqty(p,i,j) + dphase(p,i,j)
                 vol(i,j) = vol(i,j) + dphase(p,i,j) * mvol(p)
             End If
         End If
         
         
       End Do
     End Do
   End Do
!   !$OMP END DO
!   !$OMP END PARALLEL
end subroutine update_mlvl_phases