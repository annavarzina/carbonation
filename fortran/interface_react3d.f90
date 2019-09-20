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

SUBROUTINE hom_reaction(f,conc,nodetype,phaseqty,ceq,k1,k2,lz,ly,lx)
 !single mineral reaction is implemented 
    IMPLICIT NONE
    REAL(8),INTENT(INOUT):: f(lz,ly,lx,7),phaseqty(lz,ly,lx)
    REAL(8),INTENT(IN):: conc(lz,ly,lx),ceq,nodetype(lz,ly,lx),k1,k2
    INTEGER,INTENT(IN):: lz,ly,lx
    INTEGER:: i,j, k
    REAL(8)::w1sk,w2sk,creq, w1, w2
    w1 = 1.d0/7.d0
    w2 = 1.d0/7.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,w1sk,w2sk,creq) 
    !$OMP DO SCHEDULE(static)
    DO k=1,lx
        DO j=1,ly
            DO i=1,lz
                IF (nodetype(i,j,k)<=0) THEN
                    IF (phaseqty(i,j,k)>0) THEN
                        creq=-1.d0* k1 * (k2* conc(i,j,k)-ceq)
                        IF(creq <phaseqty(i,j,k)) THEN
                            w1sk=creq *w1
                            w2sk=creq *w2
                    	       phaseqty(i,j,k)=phaseqty(i,j,k)-creq
                        ELSE
                            w1sk=phaseqty(i,j,k) *w1
                            w2sk=phaseqty(i,j,k) *w2 
                    	       phaseqty(i,j,k)=0.d0
                        END IF
                    ELSE
                        w1sk=0.d0
                	       w2sk=0.d0
                    END IF
                    f(i,j,k,1)= f(i,j,k,1) + w1sk
                    f(i,j,k,2)= f(i,j,k,2) + w2sk
                    f(i,j,k,3)= f(i,j,k,3) + w2sk
                    f(i,j,k,4)= f(i,j,k,4) + w2sk
                    f(i,j,k,5)= f(i,j,k,5) + w2sk 
                    f(i,j,k,6)= f(i,j,k,6) + w2sk
                    f(i,j,k,7)= f(i,j,k,7) + w2sk 
                END IF
            END DO         
        END DO     
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END SUBROUTINE hom_reaction

SUBROUTINE het_reaction(f,conc,nodetype,phaseqty,tau,ceq,k1,k2,D0,lz,ly,lx)
    !this imposes reactions at the boundary of solid-fluid interface
    !to be applied after propogation step with bounce back
    IMPLICIT NONE
    REAL(8),INTENT(INOUT)::f(lz,ly,lx,7),phaseqty(lz,ly,lx)
    REAL(8),INTENT(IN):: conc(lz,ly,lx),ceq,tau,k1,k2,D0,nodetype(lz,ly,lx)
    INTEGER,INTENT(IN):: lz,ly,lx
    INTEGER:: i,j,k,l,nextX,nextY,nextZ,ex(7),ey(7), ez(7)
    REAL(8):: cs2,hetrxn,ftemp
    ex=(/ 0, 1, 0, 0,-1, 0, 0/)
    ey=(/ 0, 0, 1, 0, 0,-1, 0/)	
    ez=(/ 0, 0, 0, 1, 0, 0,-1/)	
    cs2 = 1.d0/3.d0
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,nextX,nextY,nextZ,ftemp) 
    !$OMP DO SCHEDULE(static)
    DO k=1,lx
        DO j=1,ly
            DO i=1,lz
                IF (nodetype(i,j,k)==1) THEN
                    DO l=2,7
                        nextX= k+ex(l)
                        nextY= j-ey(l)
                        nextZ= i+ez(l)
                        IF (nextX>0 .AND. nextX<lx+1  .AND. &
                            & nextY>0 .AND. nextY<ly+1 .AND. &
                            & nextZ>0 .AND. nextZ<lz+1) THEN
                            IF (nodetype(nextZ,nextY,nextX)<=0) THEN
                                ftemp = f(nextZ,nextY,nextX,l)
                                f(nextZ,nextY,nextX,l) =  hetrxn(ftemp,ceq,k1,k2,tau,cs2,D0)
                                phaseqty(i,j,k)= phaseqty(i,j,k)-(f(nextZ,nextY,nextX,l)-ftemp)
                            END IF
                        END IF
                    END DO 
                END IF
            END DO
        END DO  
    END DO  
    !$OMP END DO
    !$OMP END PARALLEL       
END SUBROUTINE  het_reaction


SUBROUTINE update_nodetype(nodetype,phaseqty,mv,veff,thres,lz,ly,lx)
    !this imposes reactions at the boundary of solid-fluid interface
    !to be applied after propogation step with bounce back
    IMPLICIT NONE
    REAL(8),INTENT(INOUT)::nodetype(lz,ly,lx)
    REAL(8),INTENT(IN)::phaseqty(lz,ly,lx),thres,mv,veff
    REAL(8):: frac
    INTEGER,INTENT(IN):: lz,ly,lx
    INTEGER:: i,j,k
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,frac) 
    !$OMP DO SCHEDULE(static)
    DO k=1,lx
        DO j=1,ly
            DO i=1,lz
                IF (nodetype(i,j,k)==1) THEN
                    frac = (phaseqty(i,j,k) * mv )/veff
                    IF(frac<=thres) nodetype(i,j,k) = -1.0
                END IF
            END DO
        END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL      
END SUBROUTINE update_nodetype


