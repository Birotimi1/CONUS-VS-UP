SUBROUTINE BWR(SUMV,SUMPSV,SUMSH,Nks,W)
! VTI adaptation of Garcia-Jerez HV-DFA damped body-wave contour segment.
! Separate PSV/SH meshes end at their own half-space shear branch point.
! Preserve upstream endpoint/rectangle convention, including k=0.1*dk.
USE hv_vti
USE GLOBALES,ONLY:I,twopi,SHDAMP,PSVDAMP,EmatrixNorma
USE MARC,ONLY:H
IMPLICIT NONE
INTEGER,INTENT(IN)::Nks
REAL,INTENT(IN)::W
REAL(LONG_FLOAT),INTENT(OUT)::SUMV,SUMPSV,SUMSH
REAL(LONG_FLOAT)::dkp,dks,kp,ks,t,expo,norma,hpsv,v,sh
INTEGER::ik,j
COMPLEX(LONG_CMPLX)::yp(4,2),ys(2,1),dp(4,2),ds(2,1),lp(4,4),ip(4,4),ls(2,2),is(2,2)
COMPLEX(LONG_CMPLX)::q(2,2),ep(4,4),es(2,2),gam,nu,nsh,wp,ws,denom
! Upstream conjugate sqrt is equivalent to outgoing roots at w+i*damp*w.
! Its extra -w column factor cancels in the surface Green function.
wp=cmplx(W,W*PSVDAMP,kind=long_cmplx);ws=cmplx(W,W*SHDAMP,kind=long_cmplx)
dkp=DBLE(W)*sqrt(rho(ncapas)/hl(ncapas))/Nks
dks=DBLE(W)*sqrt(rho(ncapas)/hn(ncapas))/Nks
SUMV=0;SUMPSV=0;SUMSH=0
do ik=0,Nks
  t=DBLE(ik);if(ik==0)t=0.1d0
  kp=t*dkp;ks=t*dks;j=ncapas
  call vti_psv(j,wp,kp,gam,nu,lp)
  nsh=vti_sh(j,ws,ks);call LmatrixSH(ls,hl(j),nsh)
  yp=lp(:,1:2);ys(:,1)=ls(:,1)
  do j=ncapas-1,1,-1
    call vti_psv(j,wp,kp,gam,nu,lp,ip)
    nsh=vti_sh(j,ws,ks)
    call LmatrixSH(ls,hl(j),nsh);call LmatrixInv2(is,hl(j),nsh)
    dp=matmul(ip,yp);ds=matmul(is,ys)
    norma=1/sqrt(real(dot_product(dp(:,1),dp(:,1)))*real(dot_product(dp(:,2),dp(:,2))))
    q=reshape((/dp(2,2),-dp(2,1),-dp(1,2),dp(1,1)/)*norma,(/2,2/))
    dp=matmul(dp,q);dp(1,2)=0;dp(2,1)=0
    ds=ds/sqrt(real(dot_product(ds(:,1),ds(:,1))))
    call EmatrixNorma(ep,expo,-H(j),gam,nu,3)
    call EmatrixNorma2(es,expo,-H(j),nsh)
    yp=matmul(matmul(lp,ep),dp);ys=matmul(matmul(ls,es),ds)
  enddo
  sh=real(I*ks*ys(1,1)/ys(2,1))
  yp(:,1)=yp(:,1)/maxval(abs(yp(:,1)));yp(:,2)=yp(:,2)/maxval(abs(yp(:,2)))
  denom=-yp(3,2)*yp(4,1)+yp(3,1)*yp(4,2)
  v=real(I*kp*(yp(2,2)*yp(3,1)-yp(2,1)*yp(3,2))/denom)
  hpsv=real(I*kp*(-yp(1,2)*yp(4,1)+yp(1,1)*yp(4,2))/denom)
  SUMV=SUMV-v*dkp/twopi;SUMPSV=SUMPSV-hpsv*dkp/(2*twopi)
  SUMSH=SUMSH-sh*dks/(2*twopi)
enddo
END SUBROUTINE
