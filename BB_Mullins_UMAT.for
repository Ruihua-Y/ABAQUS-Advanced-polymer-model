!DIR$ FREEFORM

subroutine umat(stress,statev,ddsdde,sse,spd,scd, &
     rpl,ddsddt,drplde,drpldt,stran,dstran,time,dtime,temp,dtemp, &
     predef,dpred,cmname,ndi,nshr,ntens,nstatv,props,nprops,coords, &
     drot,pnewdt,celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

  include 'aba_param.inc'

  character*80 cmname
  integer ndi,nshr,ntens,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc
  integer i,j

  double precision stress(ntens),statev(nstatv),ddsdde(ntens,ntens)
  double precision sse,spd,scd,rpl,ddsddt(ntens),drplde(ntens),drpldt
  double precision stran(ntens),dstran(ntens),time(2),dtime,temp,dtemp
  double precision predef(*),dpred(*),props(nprops),coords(3),drot(3,3)
  double precision pnewdt,celent,dfgrd0(3,3),dfgrd1(3,3)

  double precision kappa,mu_eq,mu_v,n_eq,n_v,gamma_ref,c_exp,m_exp
  double precision reg_eps,newton_tol
  integer newton_maxit,ierr

  double precision r_dmg,uhat_dmg,beta_dmg,amin_dmg,b_dmg
  logical use_damage

  double precision ciinv_n(3,3),ciinv_np1(3,3)
  double precision sigma(3,3),sigma_p(3,3),I3(3,3),deps(3,3),fpert(3,3)
  double precision sv_base(6),sv_pert(6),h
  double precision udev,udev_p,eta,eta_p,alpha_d,alpha_p
  double precision umax_n,umin_n,umax_new,umin_new,umax_p,umin_p
  logical initialized

  call ident3(I3)

  do i=1,ntens
     stress(i)=0.d0
     ddsddt(i)=0.d0
     drplde(i)=0.d0
     do j=1,ntens
        ddsdde(i,j)=0.d0
     end do
  end do
  sse=0.d0
  spd=0.d0
  scd=0.d0
  rpl=0.d0
  drpldt=0.d0

  if (ndi.ne.3 .or. nshr.ne.3 .or. ntens.ne.6) then
     return
  end if

  if (nprops.lt.8) then
     return
  end if

  kappa     = props(1)
  mu_eq     = props(2)
  mu_v      = props(3)
  n_eq      = props(4)
  n_v       = props(5)
  gamma_ref = props(6)
  c_exp     = props(7)
  m_exp     = props(8)

  reg_eps = 1.d-6
  if (nprops.ge.9) reg_eps = max(props(9),1.d-12)

  newton_tol = 1.d-10
  if (nprops.ge.10) newton_tol = max(props(10),1.d-14)

  newton_maxit = 30
  if (nprops.ge.11) newton_maxit = max(5,int(props(11)))

  use_damage = (nprops.ge.16 .and. nstatv.ge.10)
  r_dmg     = 1.d30
  uhat_dmg  = 1.d30
  beta_dmg  = 0.d0
  amin_dmg  = 1.d0
  b_dmg     = 0.d0
  if (nprops.ge.16) then
     r_dmg    = max(props(12),1.d0+1.d-12)
     uhat_dmg = max(props(13),1.d-30)
     beta_dmg = max(props(14),0.d0)
     amin_dmg = max(min(props(15),1.d0),0.d0)
     b_dmg    = max(props(16),0.d0)
  end if

  call state_to_ciinv(statev,nstatv,ciinv_n,initialized)
  if (.not.initialized) call ident3(ciinv_n)

  umax_n = 0.d0
  umin_n = 0.d0
  if (nstatv.ge.8) then
     umax_n = statev(7)
     umin_n = statev(8)
  end if

  call bb_update_point(dfgrd1,ciinv_n,dtime,kappa,mu_eq,mu_v,n_eq,n_v, &
       gamma_ref,c_exp,m_exp,reg_eps,newton_tol,newton_maxit, &
       sigma,udev,ciinv_np1,ierr)

  if (ierr.ne.0) then
     pnewdt = 0.5d0
  end if

  if (use_damage) then
     call or_damage_update(udev,umax_n,umin_n,r_dmg,uhat_dmg, &
          beta_dmg,amin_dmg,b_dmg,eta,alpha_d,umax_new,umin_new)
  else
     eta = 1.d0
     alpha_d = 1.d0
     umax_new = umax_n
     umin_new = umin_n
  end if

  call scale_dev_stress(sigma,eta)

  call sym_to_voigt6(sigma,stress)
  call ciinv_to_state(ciinv_np1,statev,nstatv)
  if (nstatv.ge.8) then
     statev(7)=umax_new
     statev(8)=umin_new
  end if
  if (nstatv.ge.10) then
     statev(9)=eta
     statev(10)=alpha_d
  end if

  h = 1.d-7
  do j=1,6
     call make_deps_voigt(j,h,deps)
     fpert = matmul(I3 + deps, dfgrd1)

     call bb_update_point(fpert,ciinv_n,dtime,kappa,mu_eq,mu_v,n_eq, &
          n_v,gamma_ref,c_exp,m_exp,reg_eps,newton_tol,newton_maxit, &
          sigma_p,udev_p,ciinv_np1,ierr)

     if (use_damage) then
        call or_damage_update(udev_p,umax_n,umin_n,r_dmg,uhat_dmg, &
             beta_dmg,amin_dmg,b_dmg,eta_p,alpha_p,umax_p,umin_p)
     else
        eta_p = 1.d0
     end if
     call scale_dev_stress(sigma_p,eta_p)

     call sym_to_voigt6(sigma_p,sv_pert)
     call sym_to_voigt6(sigma,sv_base)
     do i=1,6
        ddsdde(i,j) = (sv_pert(i)-sv_base(i))/h
     end do
  end do

  do i=1,6
     do j=i+1,6
        ddsdde(i,j)=0.5d0*(ddsdde(i,j)+ddsdde(j,i))
        ddsdde(j,i)=ddsdde(i,j)
     end do
  end do

  return
end subroutine umat

subroutine bb_update_point(F,ciinv_n,dtime,kappa,mu_eq,mu_v,n_eq,n_v, &
     gamma_ref,c_exp,m_exp,reg_eps,newton_tol,newton_maxit, &
     sigma,udev,ciinv_np1,ierr)

  implicit none

  integer i,j,a,it,newton_maxit,ierr
  double precision F(3,3),ciinv_n(3,3),sigma(3,3),ciinv_np1(3,3)
  double precision dtime,kappa,mu_eq,mu_v,n_eq,n_v,gamma_ref,c_exp,m_exp
  double precision reg_eps,newton_tol,udev

  double precision Jdet,detF,detCi
  double precision I3(3,3),Fbar(3,3),Fbar_inv(3,3),Ft(3,3),Fbar_t(3,3)
  double precision b(3,3),bbar(3,3),be_tr(3,3),be(3,3)
  double precision taue_bar(3,3),tauv_bar(3,3),tau_iso(3,3),tau(3,3)
  double precision eval_tr(3),evec(3,3),lam2(3),eps(3),eps_tr(3),epsp(3)
  double precision r(3),rp(3),dr(3),jac(3,3),de(3),rhs(3)
  double precision taup(3),devtau(3),devtau_p(3)
  double precision tau_eff,tau_eff_p,gamma_dot,gamma_dot_p
  double precision hnum,trbbar,lam_r_eq,qeq,trbe,lam_r_v,qv
  double precision Ci_n(3,3),tmp(3,3),mean_eps,lam_chain_i,trCi
  double precision det3,norm3,inv_langevin_ratio
  double precision psiA,psiA0,lamr0,ne_guard,nv_guard
  logical ok

  ierr = 0
  call ident3(I3)

  detF = det3(F)
  if (detF.le.1.d-16) then
     ierr = 1
     sigma = 0.d0
     udev = 0.d0
     ciinv_np1 = ciinv_n
     return
  end if
  Jdet = detF

  Fbar = Jdet**(-1.d0/3.d0) * F
  call transpose33(Fbar,Fbar_t)
  be_tr = matmul(matmul(Fbar,ciinv_n),Fbar_t)
  be_tr = 0.5d0*(be_tr + transpose(be_tr))

  call inv3(ciinv_n,Ci_n,ok)
  if (.not.ok) then
     ierr=1
     sigma=0.d0
     udev=0.d0
     ciinv_np1=ciinv_n
     return
  end if
  trCi = Ci_n(1,1)+Ci_n(2,2)+Ci_n(3,3)
  lam_chain_i = sqrt(max(trCi/3.d0,1.d-16))

  call spectral_sym3(be_tr,eval_tr,evec)
  do a=1,3
     eval_tr(a) = max(eval_tr(a),1.d-16)
     eps_tr(a) = 0.5d0*log(eval_tr(a))
     eps(a) = eps_tr(a)
  end do

  hnum = 1.d-8
  do it=1,newton_maxit
     call bb_residual_principal(eps,eps_tr,lam_chain_i,dtime,mu_v,n_v, &
          gamma_ref,c_exp,m_exp,reg_eps,r,taup,devtau,tau_eff,gamma_dot)
     if (norm3(r).lt.newton_tol) exit

     do j=1,3
        epsp = eps
        epsp(j) = epsp(j) + hnum
        call bb_residual_principal(epsp,eps_tr,lam_chain_i,dtime,mu_v, &
             n_v,gamma_ref,c_exp,m_exp,reg_eps,rp,taup,devtau_p, &
             tau_eff_p,gamma_dot_p)
        do i=1,3
           jac(i,j) = (rp(i)-r(i))/hnum
        end do
     end do

     do i=1,3
        rhs(i) = -r(i)
     end do
     call solve3x3(jac,rhs,de,ok)
     if (.not.ok) then
        ierr=1
        ciinv_np1=ciinv_n
        sigma=0.d0
        udev=0.d0
        return
     end if

     eps = eps + de
     mean_eps = (eps(1)+eps(2)+eps(3))/3.d0
     eps = eps - mean_eps
     if (norm3(de).lt.newton_tol) exit
  end do
  if (it.gt.newton_maxit) ierr = 2

  do a=1,3
     lam2(a) = exp(2.d0*eps(a))
  end do
  call rebuild_from_spectral(evec,lam2,be)

  nv_guard = max(n_v,1.d-12)
  trbe = be(1,1)+be(2,2)+be(3,3)
  lam_r_v = sqrt(max(trbe/(3.d0*nv_guard),1.d-16))
  qv = inv_langevin_ratio(lam_r_v)
  tauv_bar = (mu_v/3.d0) * qv * be

  call transpose33(F,Ft)
  b = matmul(F,Ft)
  bbar = Jdet**(-2.d0/3.d0) * b
  trbbar = bbar(1,1)+bbar(2,2)+bbar(3,3)
  ne_guard = max(n_eq,1.d-12)
  lam_r_eq = sqrt(max(trbbar/(3.d0*ne_guard),1.d-16))
  qeq = inv_langevin_ratio(lam_r_eq)
  taue_bar = (mu_eq/3.d0) * qeq * bbar

  call chain8_energy(lam_r_eq,mu_eq,ne_guard,psiA)
  lamr0 = 1.d0/sqrt(max(ne_guard,1.d0))
  call chain8_energy(lamr0,mu_eq,ne_guard,psiA0)
  udev = psiA - psiA0

  call dev3(taue_bar+tauv_bar,tau_iso)
  tau = tau_iso
  do i=1,3
     tau(i,i) = tau(i,i) + kappa*(Jdet-1.d0)
  end do
  sigma = tau / Jdet
  sigma = 0.5d0*(sigma + transpose(sigma))

  call inv3(Fbar,Fbar_inv,ok)
  if (.not.ok) then
     ierr=3
     ciinv_np1=ciinv_n
     return
  end if
  ciinv_np1 = matmul(matmul(Fbar_inv,be),transpose(Fbar_inv))
  ciinv_np1 = 0.5d0*(ciinv_np1 + transpose(ciinv_np1))
  detCi = det3(ciinv_np1)
  if (detCi.gt.1.d-16) ciinv_np1 = detCi**(-1.d0/3.d0) * ciinv_np1

  return
end subroutine bb_update_point

subroutine bb_residual_principal(eps,eps_tr,lam_chain_i,dtime,mu_v,n_v, &
     gamma_ref,c_exp,m_exp,reg_eps,r,tau,devtau,tau_eff,gamma_dot)

  implicit none

  integer a
  double precision eps(3),eps_tr(3),lam_chain_i,dtime,mu_v,n_v
  double precision gamma_ref,c_exp,m_exp,reg_eps,r(3),tau(3),devtau(3)
  double precision tau_eff,gamma_dot
  double precision lam(3),lam2(3),i1e,lam_r,qv,mtau,fac,nv_guard
  double precision inv_langevin_ratio

  do a=1,3
     lam(a) = exp(eps(a))
     lam2(a) = lam(a)*lam(a)
  end do
  i1e = lam2(1)+lam2(2)+lam2(3)
  nv_guard = max(n_v,1.d-12)
  lam_r = sqrt(max(i1e/(3.d0*nv_guard),1.d-16))
  qv = inv_langevin_ratio(lam_r)

  do a=1,3
     tau(a) = (mu_v/3.d0) * qv * lam2(a)
  end do
  mtau = (tau(1)+tau(2)+tau(3))/3.d0
  do a=1,3
     devtau(a) = tau(a)-mtau
  end do
  tau_eff = sqrt(max(devtau(1)*devtau(1)+devtau(2)*devtau(2)+ &
       devtau(3)*devtau(3),1.d-30))

  fac = max(lam_chain_i - 1.d0 + reg_eps, reg_eps)
  gamma_dot = gamma_ref * fac**c_exp * tau_eff**m_exp

  do a=1,3
     r(a) = eps(a) - eps_tr(a) + dtime*gamma_dot*devtau(a)/tau_eff
  end do

  return
end subroutine bb_residual_principal

subroutine or_damage_update(udev,umax_n,umin_n,r,uhat,beta,amin,b, &
     eta,alpha_d,umax_new,umin_new)

  implicit none

  double precision udev,umax_n,umin_n,r,uhat,beta,amin,b
  double precision eta,alpha_d,umax_new,umin_new
  double precision tiny,den,arg

  tiny = 1.d-30

  if (umax_n.lt.tiny) then
     alpha_d = 1.d0
  else
     alpha_d = max(amin,(min(umin_n,umax_n)/umax_n)**b)
  end if

  umax_new = max(umax_n, alpha_d*udev)

  den = uhat + beta*umax_new
  arg = (umax_new - alpha_d*udev)/max(den,tiny)
  eta = 1.d0 - (1.d0/r)*erf(arg)

  if (alpha_d*udev.ge.umax_n) then
     umin_new = umax_new
  else
     umin_new = min(alpha_d*udev, umin_n)
  end if

  return
end subroutine or_damage_update

subroutine scale_dev_stress(sig,eta)

  implicit none

  integer i
  double precision sig(3,3),eta,mean

  mean = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
  do i=1,3
     sig(i,i) = eta*sig(i,i) + (1.d0-eta)*mean
  end do
  sig(1,2)=eta*sig(1,2)
  sig(2,1)=eta*sig(2,1)
  sig(1,3)=eta*sig(1,3)
  sig(3,1)=eta*sig(3,1)
  sig(2,3)=eta*sig(2,3)
  sig(3,2)=eta*sig(3,2)

  return
end subroutine scale_dev_stress

subroutine chain8_energy(lamr,mu,n,psi)

  implicit none

  double precision lamr,mu,n,psi,beta,x2,lns

  x2 = lamr*lamr
  if (x2.gt.0.999999d0) x2 = 0.999999d0
  beta = lamr*(3.d0-x2)/(1.d0-x2)

  if (beta.lt.0.05d0) then

     lns = -(beta*beta/6.d0)*(1.d0+beta*beta/30.d0+(beta**4)/472.5d0)
  else

     lns = log(2.d0*beta) - beta - log(1.d0-exp(-2.d0*beta))
  end if

  psi = mu*n*(lamr*beta+lns)

  return
end subroutine chain8_energy

subroutine make_deps_voigt(j,h,deps)
  implicit none
  integer j
  double precision h,deps(3,3)
  deps = 0.d0
  if (j.eq.1) deps(1,1)=h
  if (j.eq.2) deps(2,2)=h
  if (j.eq.3) deps(3,3)=h
  if (j.eq.4) then
     deps(1,2)=0.5d0*h
     deps(2,1)=0.5d0*h
  end if
  if (j.eq.5) then
     deps(1,3)=0.5d0*h
     deps(3,1)=0.5d0*h
  end if
  if (j.eq.6) then
     deps(2,3)=0.5d0*h
     deps(3,2)=0.5d0*h
  end if
  return
end subroutine make_deps_voigt

subroutine state_to_ciinv(statev,nstatv,ciinv,initialized)
  implicit none
  integer nstatv
  double precision statev(nstatv),ciinv(3,3),snorm
  logical initialized

  call ident3(ciinv)
  if (nstatv.lt.6) then
     initialized = .false.
     return
  end if

  snorm = abs(statev(1))+abs(statev(2))+abs(statev(3))+ &
       abs(statev(4))+abs(statev(5))+abs(statev(6))
  if (snorm.lt.1.d-16) then
     initialized = .false.
     return
  end if

  initialized = .true.
  ciinv(1,1)=statev(1)
  ciinv(2,2)=statev(2)
  ciinv(3,3)=statev(3)
  ciinv(1,2)=statev(4)
  ciinv(2,1)=statev(4)
  ciinv(1,3)=statev(5)
  ciinv(3,1)=statev(5)
  ciinv(2,3)=statev(6)
  ciinv(3,2)=statev(6)
  ciinv = 0.5d0*(ciinv + transpose(ciinv))
  return
end subroutine state_to_ciinv

subroutine ciinv_to_state(ciinv,statev,nstatv)
  implicit none
  integer nstatv
  double precision ciinv(3,3),statev(nstatv)
  if (nstatv.lt.6) return
  statev(1)=ciinv(1,1)
  statev(2)=ciinv(2,2)
  statev(3)=ciinv(3,3)
  statev(4)=ciinv(1,2)
  statev(5)=ciinv(1,3)
  statev(6)=ciinv(2,3)
  return
end subroutine ciinv_to_state

subroutine sym_to_voigt6(a,v)
  implicit none
  double precision a(3,3),v(6)
  v(1)=a(1,1)
  v(2)=a(2,2)
  v(3)=a(3,3)
  v(4)=a(1,2)
  v(5)=a(1,3)
  v(6)=a(2,3)
  return
end subroutine sym_to_voigt6

subroutine ident3(I3)
  implicit none
  double precision I3(3,3)
  I3 = 0.d0
  I3(1,1)=1.d0
  I3(2,2)=1.d0
  I3(3,3)=1.d0
  return
end subroutine ident3

subroutine transpose33(a,at)
  implicit none
  double precision a(3,3),at(3,3)
  integer i,j
  do i=1,3
     do j=1,3
        at(i,j)=a(j,i)
     end do
  end do
  return
end subroutine transpose33

subroutine dev3(a,adev)
  implicit none
  double precision a(3,3),adev(3,3),m
  integer i
  adev = a
  m = (a(1,1)+a(2,2)+a(3,3))/3.d0
  do i=1,3
     adev(i,i)=adev(i,i)-m
  end do
  return
end subroutine dev3

double precision function det3(a)
  implicit none
  double precision a(3,3)
  det3 = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
       - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
       + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  return
end function det3

subroutine inv3(a,ainv,ok)
  implicit none
  double precision a(3,3),ainv(3,3),d
  double precision det3
  logical ok
  d = det3(a)
  if (abs(d).lt.1.d-16) then
     ainv=0.d0
     ok=.false.
     return
  end if
  ok=.true.
  ainv(1,1)= (a(2,2)*a(3,3)-a(2,3)*a(3,2))/d
  ainv(1,2)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))/d
  ainv(1,3)= (a(1,2)*a(2,3)-a(1,3)*a(2,2))/d
  ainv(2,1)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))/d
  ainv(2,2)= (a(1,1)*a(3,3)-a(1,3)*a(3,1))/d
  ainv(2,3)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))/d
  ainv(3,1)= (a(2,1)*a(3,2)-a(2,2)*a(3,1))/d
  ainv(3,2)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))/d
  ainv(3,3)= (a(1,1)*a(2,2)-a(1,2)*a(2,1))/d
  return
end subroutine inv3

subroutine spectral_sym3(a,eval,evec)
  implicit none
  integer i,j,p,q,iter
  double precision a(3,3),eval(3),evec(3,3),b(3,3)
  double precision maxoff,theta,c,s,bpp,bqq,bpq,bip,biq,vip,viq,tol

  tol = 1.d-12
  b = 0.5d0*(a + transpose(a))
  call ident3(evec)

  do iter=1,60
     maxoff = 0.d0
     p=1
     q=2
     do i=1,3
        do j=i+1,3
           if (abs(b(i,j)).gt.maxoff) then
              maxoff = abs(b(i,j))
              p=i
              q=j
           end if
        end do
     end do
     if (maxoff.lt.tol) exit

     theta = 0.5d0*atan2(2.d0*b(p,q),b(q,q)-b(p,p))
     c = cos(theta)
     s = sin(theta)

     bpp = b(p,p)
     bqq = b(q,q)
     bpq = b(p,q)
     b(p,p)=c*c*bpp - 2.d0*s*c*bpq + s*s*bqq
     b(q,q)=s*s*bpp + 2.d0*s*c*bpq + c*c*bqq
     b(p,q)=0.d0
     b(q,p)=0.d0

     do i=1,3
        if (i.ne.p .and. i.ne.q) then
           bip = b(i,p)
           biq = b(i,q)
           b(i,p)=c*bip - s*biq
           b(p,i)=b(i,p)
           b(i,q)=s*bip + c*biq
           b(q,i)=b(i,q)
        end if
     end do

     do i=1,3
        vip = evec(i,p)
        viq = evec(i,q)
        evec(i,p)=c*vip - s*viq
        evec(i,q)=s*vip + c*viq
     end do
  end do

  eval(1)=b(1,1)
  eval(2)=b(2,2)
  eval(3)=b(3,3)
  call sort_eigs(eval,evec)
  return
end subroutine spectral_sym3

subroutine sort_eigs(eval,evec)
  implicit none
  integer i,j,k
  double precision eval(3),evec(3,3),tmp,vec(3)
  do i=1,2
     k=i
     do j=i+1,3
        if (eval(j).gt.eval(k)) k=j
     end do
     if (k.ne.i) then
        tmp=eval(i)
        eval(i)=eval(k)
        eval(k)=tmp
        vec=evec(:,i)
        evec(:,i)=evec(:,k)
        evec(:,k)=vec
     end if
  end do
  return
end subroutine sort_eigs

subroutine rebuild_from_spectral(evec,lam2,a)
  implicit none
  integer i,j,k
  double precision evec(3,3),lam2(3),a(3,3)
  a=0.d0
  do k=1,3
     do i=1,3
        do j=1,3
           a(i,j)=a(i,j)+lam2(k)*evec(i,k)*evec(j,k)
        end do
     end do
  end do
  a=0.5d0*(a+transpose(a))
  return
end subroutine rebuild_from_spectral

subroutine solve3x3(A,b,x,ok)
  implicit none
  integer i,j,k,piv
  double precision A(3,3),b(3),x(3),M(3,4),tmp,f,maxv
  logical ok

  M(:,1:3)=A
  M(:,4)=b

  do k=1,3
     piv=k
     maxv=abs(M(k,k))
     do i=k+1,3
        if (abs(M(i,k)).gt.maxv) then
           maxv=abs(M(i,k))
           piv=i
        end if
     end do
     if (maxv.lt.1.d-16) then
        ok=.false.
        x=0.d0
        return
     end if
     if (piv.ne.k) then
        do j=k,4
           tmp=M(k,j)
           M(k,j)=M(piv,j)
           M(piv,j)=tmp
        end do
     end if

     tmp=M(k,k)
     do j=k,4
        M(k,j)=M(k,j)/tmp
     end do
     do i=1,3
        if (i.ne.k) then
           f=M(i,k)
           do j=k,4
              M(i,j)=M(i,j)-f*M(k,j)
           end do
        end if
     end do
  end do

  x(1)=M(1,4)
  x(2)=M(2,4)
  x(3)=M(3,4)
  ok=.true.
  return
end subroutine solve3x3

double precision function norm3(v)
  implicit none
  double precision v(3)
  norm3 = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
  return
end function norm3

double precision function inv_langevin_ratio(lam_r)
  implicit none
  double precision lam_r,x2,den
  x2 = lam_r*lam_r
  if (x2.lt.1.d-14) then
     inv_langevin_ratio = 3.d0
     return
  end if
  if (x2.gt.1.d0-1.d-8) x2 = 1.d0-1.d-8
  den = 1.d0 - x2
  inv_langevin_ratio = (3.d0 - x2)/den
  return
end function inv_langevin_ratio

