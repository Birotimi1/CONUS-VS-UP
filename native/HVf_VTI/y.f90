
      REAL FUNCTION Y(SLOW)
      USE Marc,ONLY:ISRAYLEIGH
      USE Globales,ONLY:long_float,long_cmplx
      REAL,INTENT(IN)::SLOW
      REAL(long_float)YL,YR
      IF(ISRAYLEIGH)THEN
        Y=real(YR(SLOW),kind=kind(Y))
      ELSE
        Y=real(YL(SLOW),kind=kind(Y))
      ENDIF
      END FUNCTION

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!                R A Y L E I G H                 !!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function yR(slow)
      use hv_vti
      use Marc,ONLY:G_OMEGA,ALFA,BTA,H,MU,NCAPAS
      use Globales,ONLY:i,pi,twopi,det2x2,EmatrixNorma
      Use Types
      IMPLICIT NONE
      real(long_float)yR
      REAL,INTENT(IN)::SLOW
      complex(long_cmplx)gam_val,nu_val
      ! Local variables
      real(long_float)w,c,exposYb,expo_prov1,expo_prov2,norma
      complex(long_cmplx)Yb(4,2),L(4,4),Linv(4,4),Einv1(4,4),Einv2(4,4),Q(2,2),gam,nu,D(4,2)
      complex(long_cmplx)phaseQ,detQ
      real(long_float) hsum
      integer index
      ! Initializing
      expo_prov1=0;expo_prov2=0
      ! Body
      w=DBLE(G_OMEGA)
      c = 1.D0 / DBLE(slow)
      exposYb=0
      phaseQ=cmplx(1.d0,0.d0,kind=long_cmplx)
      ! Go!
      index=ncapas;
  call vti_psv(index,cmplx(w,0.d0,kind=long_cmplx),w/c,gam,nu,L)
      D = reshape((/1.,0.,0.,0., 0.,1.,0.,0./),(/4,2/))
      Yb=matmul(L,D)
      ! Merge identical numerical sublayers ONLY in the root function.
      ! Keep the original thin sublayers in Wang's energy integration to
      ! avoid overflow; the physical Earth is identical in both calculations.
      hsum=0.d0
      do index=ncapas-1,1,-1
        hsum=hsum+H(index)
        if(index>1)then
          if(ha(index)==ha(index-1).and.hc(index)==hc(index-1).and.hf(index)==hf(index-1).and.hl(index)==hl(index-1).and.hn(index)==hn(index-1).and.rho(index)==rho(index-1))cycle
        endif
  call vti_psv(index,cmplx(w,0.d0,kind=long_cmplx),w/c,gam,nu,L)
  call inverse4(L,Linv)
        D=matmul(Linv,Yb)
        norma=1/sqrt(real(dot_product(D(:,1),D(:,1)))*real(dot_product(D(:,2),D(:,2))))
        Q=reshape((/D(2,2),-D(2,1),-D(1,2),D(1,1)/)*norma,(/2,2/))
        D=matmul(D,Q);D(1,2)=0;D(2,1)=0
        ! det(Y*conjg(Q1*...*Qn)) = det(Y)*conjg(prod(det(Qj))).
        ! Only its phase affects zeros/sign changes. Accumulating the full
        ! Q product loses rank in deep evanescent layers and creates false
        ! roots near the surface shear branch. Remove positive magnitudes
        ! at every step without changing the secular equation's sign.
        detQ=det2x2(Q/maxval(abs(Q)))
        if(abs(detQ)==0) error stop 'Singular modal basis'
        phaseQ=phaseQ*detQ/abs(detQ)
        call EmatrixNorma(Einv1,expo_prov1,-hsum,gam,nu,1)
        Yb(:,1)=matmul(matmul(L,Einv1),D(:,1))
        call EmatrixNorma(Einv2,expo_prov2,-hsum,gam,nu,2)
        Einv2(1,1)=0
        Yb(:,2)=matmul(matmul(L,Einv2),D(:,2))
        hsum=0.d0
      enddo
      yR=real(det2x2(Yb(3:4,:))*conjg(phaseQ))
      return
      end function


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!                   L O V E                      !!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION YL(SLOW)
      use hv_vti
      USE Marc,ONLY:G_OMEGA,BTA,H,MU,NCAPAS
      USE Globales,ONLY:long_float,long_cmplx,i,pi,twopi
      implicit none
      real(long_float)yl
      real,intent(in)::slow
      complex(long_cmplx) nu_val
      ! Local variables
      integer index
      real(long_float)::hsum
      real(long_float)::w,c,expo_prov,Q!,Qacum,exposQacum,exposYb
      complex(long_cmplx)::Yb(2,1),L(2,2),Linv(2,2),nu,D(2,1),Einv1(2,2)
      ! Body
      w=G_OMEGA
      c=1/slow
      index=ncapas
      nu=vti_sh(index,cmplx(w,0.d0,kind=long_cmplx),w/c)
      CALL LmatrixSH(L,mu(index),nu)
      D(:,1)=(/1.,0./)
      Yb=matmul(L,D)
      hsum=0.d0
      do index=ncapas-1,1,-1
        hsum=hsum+H(index)
        if(index>1)then
          if(ha(index)==ha(index-1).and.hc(index)==hc(index-1).and.hf(index)==hf(index-1).and.hl(index)==hl(index-1).and.hn(index)==hn(index-1).and.rho(index)==rho(index-1))cycle
        endif
        nu=vti_sh(index,cmplx(w,0.d0,kind=long_cmplx),w/c)
        CALL LmatrixSH(L,mu(index),nu)
        CALL LmatrixInv2(Linv,mu(index),nu)
        D=matmul(Linv,Yb)
        Q=1/sqrt(real(D(1,1)*conjg(D(1,1))+D(2,1)*conjg(D(2,1))))
        D=D*Q
        CALL EmatrixNorma2(Einv1,expo_prov,-hsum,nu)
        Yb=matmul(matmul(L,Einv1),D)
        hsum=0.d0
      enddo
      yL=real(Yb(2,1))
      RETURN
      END FUNCTION
