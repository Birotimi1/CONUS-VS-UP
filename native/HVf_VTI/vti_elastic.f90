! CONUS-VS-UP VTI extension of the GPL HV-DFA propagator, November 2026.
! Same transformed state as upstream: [u,v,L(u'-k*v),F*k*u+C*v'].
! Physical vertical displacement differs by a factor i; preserve the existing
! Green-function/ellipticity convention. See README_VTI.md for derivation.
module hv_vti
  use types
  use marc, only: ncapas,alfa,bta,rho,mu
  implicit none
  logical :: vti_input=.false.
  logical :: require_complete_modes=.false.
  real(long_float) :: max_phase_step=3.141592653589793d0/32
  real(long_float), allocatable :: xi(:),phi(:),eta(:),ha(:),hc(:),hf(:),hl(:),hn(:)
contains
  function phase_change(s1,s2) result(change)
    ! Limit accumulated vertical P/S or SH phase change between search points.
    ! This resolves dense deep-Earth overtones without tiny steps throughout
    ! the evanescent range. Tightening this limit is a convergence check.
    use marc, only: h,g_omega,israyleigh
    real(long_float),intent(in)::s1,s2
    real(long_float)::change,a,c,l,f,bb,cc,sp(2),qq(2,2)
    complex(long_cmplx)::dd,g2,n2
    integer j,k
    sp=[s1,s2]; change=0
    do j=1,ncapas-1
      a=ha(j)/rho(j); c=hc(j)/rho(j); l=hl(j)/rho(j); f=hf(j)/rho(j)
      do k=1,2
        if(israyleigh)then
          bb=(f+l)**2*sp(k)**2-c*(a*sp(k)**2-1)-l*(l*sp(k)**2-1)
          cc=(a*sp(k)**2-1)*(l*sp(k)**2-1)
          dd=sqrt(cmplx(bb*bb-4*l*c*cc,0.d0,kind=long_cmplx))
          g2=(-bb+dd)/(2*l*c); n2=(-bb-dd)/(2*l*c)
          qq(:,k)=[abs(aimag(sqrt(g2))),abs(aimag(sqrt(n2)))]
        else
          qq(1,k)=sqrt(max((rho(j)-hn(j)*sp(k)**2)/hl(j),0.d0));qq(2,k)=0
        endif
      enddo
      change=change+h(j)*sum(abs(qq(:,1)-qq(:,2)))
    enddo
    change=abs(g_omega)*change
  end function

  subroutine initialize_vti()
    if(precision(1.0)<15) error stop 'Build HVf with double-precision default REAL; use build.sh'
    if (.not.allocated(xi)) then
      allocate(xi(ncapas),phi(ncapas),eta(ncapas)); xi=1;phi=1;eta=1
    endif
    allocate(ha(ncapas),hc(ncapas),hf(ncapas),hl(ncapas),hn(ncapas))
    if(any(xi<=0).or.any(phi<=0).or.any(rho<=0)) error stop 'Invalid VTI inputs'
    ha=rho*alfa**2*5/(phi+4); hc=ha*phi
    hl=rho*bta**2*3/(xi+2); hn=hl*xi; hf=eta*(ha-2*hl)
    if(any(ha<=hn).or.any((ha-hn)*hc<=hf**2)) error stop 'Unstable VTI tensor'
    mu=hl
  end subroutine

  function outgoing_root(z) result(q)
    complex(long_cmplx),intent(in)::z
    complex(long_cmplx)::q
    q=sqrt(z)
    if(real(q)<0.or.(real(q)==0.and.aimag(q)>0)) q=-q
    ! Avoid singular eigenvector matrices exactly on a branch point.
    if(abs(q)<tiny(1.d0)**0.25d0) q=cmplx(tiny(1.d0)**0.25d0,0.d0,kind=long_cmplx)
  end function

  subroutine vti_psv(j,w,k,gam,nu,lmat,linv)
    integer,intent(in)::j
    complex(long_cmplx),intent(in)::w
    real(long_float),intent(in)::k
    complex(long_cmplx),intent(out)::gam,nu,lmat(4,4)
    complex(long_cmplx),optional,intent(out)::linv(4,4)
    real(long_float)::a,c,l,f
    complex(long_cmplx)::bb,cc,dd,g2,n2,u,v,scale
    integer::branch,up,col
    a=ha(j)/rho(j);c=hc(j)/rho(j);l=hl(j)/rho(j);f=hf(j)/rho(j)
    bb=(f+l)**2*k*k-c*(a*k*k-w*w)-l*(l*k*k-w*w)
    cc=(a*k*k-w*w)*(l*k*k-w*w)
    dd=sqrt(bb*bb-4*l*c*cc)
    g2=(-bb+dd)/(2*l*c); n2=(-bb-dd)/(2*l*c)
    gam=outgoing_root(g2); nu=outgoing_root(n2)
    do branch=1,2
      if(branch==1)then
        u=k; v=(a*k*k-l*g2-w*w)/((f+l)*gam); scale=alfa(j)/w
      else
        v=k; u=(c*n2-l*k*k+w*w)/((f+l)*nu); scale=bta(j)/w
      endif
      do up=0,1
        col=branch+2*up
        lmat(1,col)=u*scale; lmat(2,col)=v*scale*(1-2*up)
        if(branch==1)then
          lmat(3,col)=hl(j)*(-gam*u-k*v)*scale*(1-2*up)
          lmat(4,col)=(hf(j)*k*u-hc(j)*gam*v)*scale
        else
          lmat(3,col)=hl(j)*(-nu*u-k*v)*scale*(1-2*up)
          lmat(4,col)=(hf(j)*k*u-hc(j)*nu*v)*scale
        endif
      enddo
    enddo
    if(present(linv)) call inverse4(lmat,linv)
  end subroutine

  function vti_sh(j,w,k) result(nu)
    integer,intent(in)::j
    complex(long_cmplx),intent(in)::w
    real(long_float),intent(in)::k
    complex(long_cmplx)::nu
    nu=outgoing_root((hn(j)*k*k-rho(j)*w*w)/hl(j))
  end function

  subroutine inverse4(a,b)
    complex(long_cmplx),intent(in)::a(4,4)
    complex(long_cmplx),intent(out)::b(4,4)
    complex(long_cmplx)::aug(4,8),row(8),pivot
    real(long_float)::s
    integer::j,k,p
    aug=0;aug(:,1:4)=a
    do j=1,4
      aug(j,j+4)=1;s=maxval(abs(a(j,:)))
      if(s==0) error stop 'Singular VTI eigenvectors'
      aug(j,:)=aug(j,:)/s
    enddo
    do j=1,4
      p=j-1+maxloc(abs(aug(j:4,j)),1)
      row=aug(j,:);aug(j,:)=aug(p,:);aug(p,:)=row
      pivot=aug(j,j)
      if(abs(pivot)<tiny(1.d0)) error stop 'Singular VTI eigenvectors'
      aug(j,:)=aug(j,:)/pivot
      do k=1,4
        if(k/=j) aug(k,:)=aug(k,:)-aug(k,j)*aug(j,:)
      enddo
    enddo
    b=aug(:,5:8)
  end subroutine

  function halfspace_secular(j,s) result(d)
    integer,intent(in)::j
    real(long_float),intent(in)::s
    real(long_float)::d
    complex(long_cmplx)::g,n,lmat(4,4)
    call vti_psv(j,cmplx(1.d0,0.d0,kind=long_cmplx),s,g,n,lmat)
    d=real(lmat(3,1)*lmat(4,2)-lmat(3,2)*lmat(4,1))/maxval(abs(lmat(3:4,1:2)))**2
  end function

  function vti_rayleigh_slow(j) result(s)
    integer,intent(in)::j
    real(long_float)::s,lo,hi,mid,fl,fm,sv
    integer::iteration
    sv=sqrt(rho(j)/hl(j));lo=sv*(1+1.d-7);hi=sv*4
    fl=halfspace_secular(j,lo)
    if(fl*halfspace_secular(j,hi)>=0) error stop 'VTI Rayleigh bound not bracketed'
    do iteration=1,60
      mid=(lo+hi)/2;fm=halfspace_secular(j,mid)
      if(fl*fm<=0)then
        hi=mid
      else
        lo=mid;fl=fm
      endif
    enddo
    s=(lo+hi)/2
  end function
end module
