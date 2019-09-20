!|-------------------------------------------------------------------------------|
!|single mineral reaction contains subroutines for hetrogenous and               |
!|homogenous mineral reaction in 2D(D2Q5) and 3D(D3Q7)                           |
!|                                                                               |
!|*******************************************************************************|
!|          DEVELOPER: Ravi Patel||EMAIL:rpatel@sckcen.be,rapatel@ugent.be       |
!|*******************************************************************************|
REAL(8) FUNCTION hetrxn(fin,ceq,k1,k2,tau,cs2,D)
   !implements D*dc/dx = -k1(k2c-ceq)
   IMPLICIT NONE
   REAL(8),Intent(IN)::fin,ceq,k1,k2,tau,cs2,D
   REAL(8):: a1,a2,a3,an,bn,wi
   a1 =k1*k2
   a2 = -D
   a3 = k1*ceq
   wi = (1.d0/cs2)
   an = wi * a1 - wi * a2 / tau
   bn = wi * a1 + wi * a2 / tau
   hetrxn = (a3-fin*bn) / an
END FUNCTION hetrxn

SUBROUTINE hom_reaction(f,conc,nodetype,phaseqty,ceq,k1,k2,ly,lx)
 !single mineral reaction is implemented 
    USE omp_lib 		
    IMPLICIT NONE
    REAL(8),INTENT(INOUT)::f(ly,lx,5),phaseqty(ly,lx)
    REAL(8),INTENT(IN):: conc(ly,lx),ceq,nodetype(ly,lx),k1,k2
    INTEGER,INTENT(IN):: ly,lx
    INTEGER:: i,j
    REAL(8)::w1sk,w2sk,creq, w1, w2
    w1 = 2.d0/6.d0
    w2 = 1.d0/6.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,w1sk,w2sk,creq) 
    !$OMP DO SCHEDULE(DYNAMIC,1)
    DO j=1,lx
      DO i=1,ly
        IF (nodetype(i,j)<=0) THEN	
              creq=-1.d0* k1 * (k2* conc(i,j)-ceq)
              !creq= ceq- conc(i,j)
	  IF (phaseqty(i,j)>0) THEN
	    IF(creq <phaseqty(i,j)) THEN
	      w1sk=creq *w1
              w2sk=creq *w2
	      phaseqty(i,j)=phaseqty(i,j)-creq
	    ELSE
	      w1sk=phaseqty(i,j) *w1
              w2sk=phaseqty(i,j) *w2 
	      phaseqty(i,j)=0.d0
            END IF
          ELSE
            w1sk=0.d0
	    w2sk=0.d0
          END IF
          f(i,j,1)= f(i,j,1) + w1sk
          f(i,j,2)= f(i,j,2) + w2sk
          f(i,j,3)= f(i,j,3) + w2sk
          f(i,j,4)= f(i,j,4) + w2sk
          f(i,j,5)= f(i,j,5) + w2sk 
        END IF
      END DO         
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END SUBROUTINE hom_reaction

SUBROUTINE het_reaction(f,conc,nodetype,phaseqty,tau,ceq,k1,k2,D0,ly,lx)
    !this imposes reactions at the boundary of solid-fluid interface
    !to be applied after propogation step with bounce back
    USE omp_lib 		
    IMPLICIT NONE
    REAL(8),INTENT(INOUT)::f(ly,lx,5),phaseqty(ly,lx)
    REAL(8),INTENT(IN):: conc(ly,lx),ceq,tau,k1,k2,D0,nodetype(ly,lx)
    INTEGER,INTENT(IN):: ly,lx
    INTEGER:: i,j,k,nextX,nextY,ex(5),ey(5)
    REAL(8):: cs2,hetrxn,ftemp
    ex= (/ 0, 1, 0,-1, 0/)
    ey= (/ 0, 0, 1, 0,-1/)
    cs2 = 1.d0/3.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,nextX,nextY,ftemp) 
    !$OMP DO SCHEDULE(static)
    DO j=1,lx
      DO i=1,ly
	IF (nodetype(i,j)==1) THEN
	  DO k=2,5
	    nextX= j+ex(k)
	    nextY= i-ey(k)
	    IF (nextX>0 .AND. nextX<lx+1  .AND. &
	      & nextY>0 .AND. nextY<ly+1) THEN
	      IF (nodetype(nextY,nextX)<=0) THEN
	        ftemp = f(nextY,nextX,k)
		f(nextY,nextX,k) =  hetrxn(ftemp,ceq,k1,k2,tau,cs2,D0)
	        phaseqty(i,j)= phaseqty(i,j)-(f(nextY,nextX,k)-ftemp)
              END IF
            END IF
          END DO 	
	END IF
      END DO
    END DO  
    !$OMP END DO
    !$OMP END PARALLEL       
END SUBROUTINE  het_reaction

SUBROUTINE update_nodetype(nodetype,phaseqty,mv,veff,thres,ly,lx)
    !this imposes reactions at the boundary of solid-fluid interface
    !to be applied after propogation step with bounce back
    USE omp_lib 		
    IMPLICIT NONE
    REAL(8),INTENT(INOUT)::nodetype(ly,lx)
    REAL(8),INTENT(IN)::phaseqty(ly,lx),thres,mv,veff
    REAL(8):: frac
    INTEGER,INTENT(IN):: ly,lx
    INTEGER:: i,j
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j) 
    !$OMP DO SCHEDULE(static)
    DO j=1,lx
      DO i=1,ly
	  IF (nodetype(i,j)==1) THEN
             frac = (phaseqty(i,j) * mv )/veff
             IF(frac<=thres) nodetype(i,j) = -1.0
	  END IF
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL      
END SUBROUTINE update_nodetype


