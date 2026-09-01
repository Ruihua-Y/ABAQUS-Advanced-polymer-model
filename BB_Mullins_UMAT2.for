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

        double precision mu,lamL,kappa,s,xi,cexp,tauBase,mexp,tauCut
        double precision reg_eps,newton_tol
        integer newton_maxit,ierr,ierrp

        double precision fv_n(3,3),fv_np1(3,3),fv_dum(3,3)
        double precision sigma(3,3),sigma_p(3,3),I3(3,3),deps(3,3),fpert(3,3)
        double precision sv_base(6),sv_pert(6),h
        double precision del_gam,del_dum,trb,chain_strain

        double precision r_dmg,uhat_dmg,beta_dmg,amin_dmg,b_dmg
        double precision udev,udev_p,eta,eta_p,alpha_d,alpha_p
        double precision umax_n,umin_n,umax_new,umin_new,umax_p,umin_p
        logical use_damage,initialized

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

        if (ndi.ne.3 .or. nshr.ne.3 .or. ntens.ne.6) return
        if (nprops.lt.9) return
        if (nstatv.lt.13) return

        mu      = props(1)
        lamL    = props(2)
        kappa   = props(3)
        s       = props(4)
        xi      = props(5)
        cexp    = props(6)
        tauBase = props(7)
        mexp    = props(8)
        tauCut  = props(9)

        reg_eps = 1.d-6
        if (nprops.ge.10) reg_eps = max(props(10),1.d-12)

        newton_tol = 1.d-10
        if (nprops.ge.11) newton_tol = max(props(11),1.d-14)

        newton_maxit = 30
        if (nprops.ge.12) newton_maxit = max(5,int(props(12)))

        use_damage = (nprops.ge.17 .and. nstatv.ge.17)
        r_dmg     = 1.d30
        uhat_dmg  = 1.d30
        beta_dmg  = 0.d0
        amin_dmg  = 1.d0
        b_dmg     = 0.d0
        if (nprops.ge.17) then
           r_dmg    = max(props(13),1.d0+1.d-12)
           uhat_dmg = max(props(14),1.d-30)
           beta_dmg = max(props(15),0.d0)
           amin_dmg = max(min(props(16),1.d0),0.d0)
           b_dmg    = max(props(17),0.d0)
        end if

        call state_to_fv(statev,nstatv,fv_n,initialized)
        if (.not.initialized) call ident3(fv_n)
        call normalize_fv(fv_n)

        umax_n = 0.d0
        umin_n = 0.d0
        if (nstatv.ge.15) then
           umax_n = statev(14)
           umin_n = statev(15)
        end if

        call bb_update_point(dfgrd1,fv_n,dtime,mu,lamL,kappa,s,xi, &
             cexp,tauBase,mexp,tauCut,reg_eps,newton_tol,newton_maxit, &
             sigma,fv_np1,del_gam,udev,ierr)

        if (ierr.ne.0) pnewdt = 0.5d0

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
        call sym_to_voigt6(sigma,sv_base)
        call fv_to_state(fv_np1,statev,nstatv)

        statev(1) = time(2)
        statev(2) = statev(2) + del_gam
        trb = dfgrd1(1,1)**2 + dfgrd1(1,2)**2 + dfgrd1(1,3)**2 + &
              dfgrd1(2,1)**2 + dfgrd1(2,2)**2 + dfgrd1(2,3)**2 + &
              dfgrd1(3,1)**2 + dfgrd1(3,2)**2 + dfgrd1(3,3)**2
        chain_strain = 0.5d0*log(max(trb/3.d0,1.d-16))
        statev(3) = chain_strain
        statev(4) = 0.d0
        if (nstatv.ge.15) then
           statev(14) = umax_new
           statev(15) = umin_new
        end if
        if (nstatv.ge.17) then
           statev(16) = eta
           statev(17) = alpha_d
        end if

        h = 1.d-7
        do j=1,6
           call make_deps_voigt(j,h,deps)
           do i=1,3
              fpert(i,1) = (I3(i,1)+deps(i,1))*dfgrd1(1,1) + &
                           (I3(i,2)+deps(i,2))*dfgrd1(2,1) + &
                           (I3(i,3)+deps(i,3))*dfgrd1(3,1)
              fpert(i,2) = (I3(i,1)+deps(i,1))*dfgrd1(1,2) + &
                           (I3(i,2)+deps(i,2))*dfgrd1(2,2) + &
                           (I3(i,3)+deps(i,3))*dfgrd1(3,2)
              fpert(i,3) = (I3(i,1)+deps(i,1))*dfgrd1(1,3) + &
                           (I3(i,2)+deps(i,2))*dfgrd1(2,3) + &
                           (I3(i,3)+deps(i,3))*dfgrd1(3,3)
           end do

           call bb_update_point(fpert,fv_n,dtime,mu,lamL,kappa,s,xi, &
                cexp,tauBase,mexp,tauCut,reg_eps,newton_tol,newton_maxit, &
                sigma_p,fv_dum,del_dum,udev_p,ierrp)

           if (use_damage) then
              call or_damage_update(udev_p,umax_n,umin_n,r_dmg,uhat_dmg, &
                   beta_dmg,amin_dmg,b_dmg,eta_p,alpha_p,umax_p,umin_p)
           else
              eta_p = 1.d0
           end if
           call scale_dev_stress(sigma_p,eta_p)

           call sym_to_voigt6(sigma_p,sv_pert)
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

      subroutine bb_update_point(F,fv_n,dtime,mu,lamL,kappa,s,xi, &
           cexp,tauBase,mexp,tauCut,reg_eps,newton_tol,newton_maxit, &
           sigma,fv_np1,del_gam,udev,ierr)

        implicit none

        integer i,j,a,it,newton_maxit,ierr,ls
        double precision F(3,3),fv_n(3,3),sigma(3,3),fv_np1(3,3)
        double precision dtime,mu,lamL,kappa,s,xi,cexp,tauBase,mexp,tauCut
        double precision reg_eps,newton_tol
        double precision del_gam,udev
        double precision lamA,psiA,psiA0,lamLg

        double precision Jdet,detF,facJ
        double precision b(3,3),bbar(3,3)
        double precision fv_inv(3,3),fe_tr(3,3),be_tr(3,3),bbar_e_tr(3,3)
        double precision evec(3,3),eval_tr(3),eps_tr(3),eps(3)
        double precision Cbar(3,3),tmp1(3,3)
        double precision r(3),rp(3),jac(3,3),de(3),rhs(3)
        double precision Nvec(3),Np(3),tau_eff,tau_effp,gamdot,gamdotp
        double precision lamV,lamVp
        double precision epst(3),rt(3),Nt(3),tau_efft,gamdott,lamVt
        double precision rnorm,rnormt,alpha,hnum,mean_eps,de_norm
        double precision lam2(3),bbar_e(3,3),expN(3),eN(3,3)
        double precision finv(3,3),m1(3,3),m2(3,3)
        double precision det3,norm3
        logical ok,converged

        ierr = 0
        del_gam = 0.d0
        udev = 0.d0

        detF = det3(F)
        if (detF.le.1.d-16) then
           ierr = 1
           sigma = 0.d0
           fv_np1 = fv_n
           return
        end if
        Jdet = detF
        facJ = Jdet**(-2.d0/3.d0)

        b = matmul(F,transpose(F))
        b = 0.5d0*(b+transpose(b))
        bbar = facJ*b

        lamLg = max(lamL,1.d-12)
        lamA  = sqrt(max((bbar(1,1)+bbar(2,2)+bbar(3,3))/3.d0,1.d-16))
        call chain8_energy(lamA/lamLg,mu,lamLg*lamLg,psiA)
        call chain8_energy(1.d0/lamLg,mu,lamLg*lamLg,psiA0)
        udev = psiA - psiA0

        call inv3(fv_n,fv_inv,ok)
        if (.not.ok) then
           ierr = 1
           sigma = 0.d0
           fv_np1 = fv_n
           return
        end if
        fe_tr = matmul(F,fv_inv)
        be_tr = matmul(fe_tr,transpose(fe_tr))
        be_tr = 0.5d0*(be_tr+transpose(be_tr))
        bbar_e_tr = facJ*be_tr
        bbar_e_tr = 0.5d0*(bbar_e_tr+transpose(bbar_e_tr))

        call spectral_sym3(bbar_e_tr,eval_tr,evec)
        do a=1,3
           eval_tr(a) = max(eval_tr(a),1.d-16)
           eps_tr(a) = 0.5d0*log(eval_tr(a))
        end do
        mean_eps = (eps_tr(1)+eps_tr(2)+eps_tr(3))/3.d0
        do a=1,3
           eps_tr(a) = eps_tr(a) - mean_eps
           eps(a) = eps_tr(a)
        end do

        tmp1 = matmul(bbar,evec)
        Cbar = matmul(transpose(evec),tmp1)

        hnum = 1.d-8
        converged = .false.
        do it=1,newton_maxit
           call bb_residual_principal(eps,eps_tr,Cbar,dtime,Jdet,s*mu,lamL, &
                xi,cexp,tauBase,mexp,tauCut,reg_eps, &
                r,Nvec,tau_eff,gamdot,lamV)
           rnorm = norm3(r)
           if (rnorm.lt.newton_tol) then
              converged = .true.
              exit
           end if

           do j=1,3
              do i=1,3
                 epst(i) = eps(i)
              end do
              epst(j) = epst(j) + hnum
              call bb_residual_principal(epst,eps_tr,Cbar,dtime,Jdet,s*mu, &
                   lamL,xi,cexp,tauBase,mexp,tauCut,reg_eps, &
                   rp,Np,tau_effp,gamdotp,lamVp)
              do i=1,3
                 jac(i,j) = (rp(i)-r(i))/hnum
              end do
           end do

           do i=1,3
              rhs(i) = -r(i)
           end do
           call solve3x3(jac,rhs,de,ok)
           if (.not.ok) then
              ierr = 1
              exit
           end if

           if (de(1).ne.de(1) .or. de(2).ne.de(2) .or. de(3).ne.de(3)) then
              ierr = 2
              exit
           end if
           if (abs(de(1)).gt.1.d50 .or. abs(de(2)).gt.1.d50 .or. &
                abs(de(3)).gt.1.d50) then
              ierr = 2
              exit
           end if
           de_norm = norm3(de)

           alpha = 1.d0
           ok = .false.
           do ls=1,10
              do i=1,3
                 epst(i) = eps(i) + alpha*de(i)
              end do
              mean_eps = (epst(1)+epst(2)+epst(3))/3.d0
              do i=1,3
                 epst(i) = epst(i) - mean_eps
              end do
              call bb_residual_principal(epst,eps_tr,Cbar,dtime,Jdet,s*mu, &
                   lamL,xi,cexp,tauBase,mexp,tauCut,reg_eps, &
                   rt,Nt,tau_efft,gamdott,lamVt)
              if (rt(1).ne.rt(1) .or. rt(2).ne.rt(2) .or. &
                   rt(3).ne.rt(3)) then
                 alpha = 0.5d0*alpha
                 cycle
              end if
              if (abs(rt(1)).gt.1.d100 .or. abs(rt(2)).gt.1.d100 .or. &
                   abs(rt(3)).gt.1.d100) then
                 alpha = 0.5d0*alpha
                 cycle
              end if
              rnormt = norm3(rt)
              if (rnormt.le.(1.d0-0.1d0*alpha)*rnorm) then
                 ok = .true.
                 exit
              end if
              alpha = 0.5d0*alpha
           end do
           if (.not.ok) then
              ierr = 2
              exit
           end if

           do i=1,3
              eps(i) = epst(i)
              r(i) = rt(i)
              Nvec(i) = Nt(i)
           end do
           tau_eff = tau_efft
           gamdot = gamdott
           lamV = lamVt

           if (alpha*de_norm.lt.100.d0*newton_tol .and. &
                rnormt.lt.100.d0*newton_tol) then
              converged = .true.
              exit
           end if
        end do
        if (it.gt.newton_maxit .and. .not.converged) ierr = 2

        if (ierr.ne.0) then

           call bb_stress(Jdet,bbar,bbar_e_tr,mu,lamL,kappa,s,sigma)
           fv_np1 = fv_n
           del_gam = 0.d0
           return
        end if

        do a=1,3
           lam2(a) = exp(2.d0*eps(a))
        end do
        call rebuild_from_spectral(evec,lam2,bbar_e)

        del_gam = dtime*gamdot
        do a=1,3
           expN(a) = exp(del_gam*Nvec(a))
        end do
        call rebuild_from_spectral(evec,expN,eN)

        call inv3(F,finv,ok)
        if (.not.ok) then
           ierr = 3
           call bb_stress(Jdet,bbar,bbar_e,mu,lamL,kappa,s,sigma)
           fv_np1 = fv_n
           del_gam = 0.d0
           return
        end if
        m1 = matmul(fv_n,finv)
        m2 = matmul(eN,F)
        fv_np1 = matmul(m1,m2)
        call normalize_fv(fv_np1)

        call bb_stress(Jdet,bbar,bbar_e,mu,lamL,kappa,s,sigma)

        return
      end subroutine bb_update_point

      subroutine bb_residual_principal(eps,eps_tr,Cbar,dtime,J,muB,lamL, &
           xi,cexp,tauBase,mexp,tauCut,reg_eps, &
           r,Nvec,tau_eff,gamdot,lamV)

        implicit none

        integer a
        double precision eps(3),eps_tr(3),Cbar(3,3),dtime,J,muB,lamL
        double precision xi,cexp,tauBase,mexp,tauCut,reg_eps
        double precision r(3),Nvec(3),tau_eff,gamdot,lamV
        double precision lam2(3),i1e,lamBe,qre,ql,prefac
        double precision devsig(3),tau2,lamV2,base,flowarg
        double precision inv_langevin_ratio

        do a=1,3
           lam2(a) = exp(2.d0*eps(a))
        end do
        i1e = lam2(1)+lam2(2)+lam2(3)
        lamBe = sqrt(max(i1e/3.d0,1.d-16))
        qre = inv_langevin_ratio(lamBe/lamL)
        ql  = inv_langevin_ratio(1.d0/lamL)
        prefac = muB*qre/(J*ql)

        do a=1,3
           devsig(a) = prefac*(lam2(a)-i1e/3.d0)
        end do
        tau2 = devsig(1)*devsig(1)+devsig(2)*devsig(2)+devsig(3)*devsig(3)
        tau_eff = sqrt(max(tau2,1.d-30))

        lamV2 = (Cbar(1,1)*exp(-2.d0*eps(1)) + &
                 Cbar(2,2)*exp(-2.d0*eps(2)) + &
                 Cbar(3,3)*exp(-2.d0*eps(3)))/3.d0
        lamV = sqrt(max(lamV2,1.d-16))

        base = max(lamV-1.d0+xi,reg_eps)
        flowarg = max(tau_eff/tauBase-tauCut,0.d0)
        gamdot = base**cexp * flowarg**mexp
        gamdot = min(gamdot,1.d15)

        do a=1,3
           Nvec(a) = devsig(a)/tau_eff
           r(a) = eps(a)-eps_tr(a)+dtime*gamdot*Nvec(a)
        end do

        return
      end subroutine bb_residual_principal

      subroutine bb_stress(J,bbar,bbar_e,mu,lamL,kappa,s,sigma)

        implicit none

        integer i
        double precision J,bbar(3,3),bbar_e(3,3),mu,lamL,kappa,s
        double precision sigma(3,3)
        double precision trbbar,trbe,lamA,lamBe,qA,qB,ql,facA,facB,p
        double precision sigA_dev(3,3),sigB_dev(3,3)
        double precision inv_langevin_ratio

        trbbar = bbar(1,1)+bbar(2,2)+bbar(3,3)
        lamA = sqrt(max(trbbar/3.d0,1.d-16))
        qA = inv_langevin_ratio(lamA/lamL)
        ql = inv_langevin_ratio(1.d0/lamL)

        trbe = bbar_e(1,1)+bbar_e(2,2)+bbar_e(3,3)
        lamBe = sqrt(max(trbe/3.d0,1.d-16))
        qB = inv_langevin_ratio(lamBe/lamL)

        call dev3(bbar,sigA_dev)
        call dev3(bbar_e,sigB_dev)

        facA = mu*qA/(J*ql)
        facB = s*mu*qB/(J*ql)
        sigma = facA*sigA_dev + facB*sigB_dev

        p = 2.d0*kappa*(J-1.d0)
        do i=1,3
           sigma(i,i) = sigma(i,i) + p
        end do
        sigma = 0.5d0*(sigma+transpose(sigma))

        return
      end subroutine bb_stress

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

      subroutine state_to_fv(statev,nstatv,fv,initialized)
        implicit none
        integer nstatv,i
        double precision statev(nstatv),fv(3,3),snorm
        logical initialized

        call ident3(fv)
        if (nstatv.lt.13) then
           initialized = .false.
           return
        end if

        snorm = 0.d0
        do i=5,13
           snorm = snorm + abs(statev(i))
        end do
        if (snorm.lt.1.d-16) then
           initialized = .false.
           return
        end if

        initialized = .true.
        fv(1,1)=statev(5)
        fv(1,2)=statev(6)
        fv(1,3)=statev(7)
        fv(2,1)=statev(8)
        fv(2,2)=statev(9)
        fv(2,3)=statev(10)
        fv(3,1)=statev(11)
        fv(3,2)=statev(12)
        fv(3,3)=statev(13)
        return
      end subroutine state_to_fv

      subroutine fv_to_state(fv,statev,nstatv)
        implicit none
        integer nstatv
        double precision fv(3,3),statev(nstatv)
        if (nstatv.lt.13) return
        statev(5)=fv(1,1)
        statev(6)=fv(1,2)
        statev(7)=fv(1,3)
        statev(8)=fv(2,1)
        statev(9)=fv(2,2)
        statev(10)=fv(2,3)
        statev(11)=fv(3,1)
        statev(12)=fv(3,2)
        statev(13)=fv(3,3)
        return
      end subroutine fv_to_state

      subroutine normalize_fv(fv)
        implicit none
        double precision fv(3,3),d,fac
        double precision det3
        d = det3(fv)
        if (d.le.1.d-16) then
           call ident3(fv)
           return
        end if
        fac = d**(-1.d0/3.d0)
        fv = fac*fv
        return
      end subroutine normalize_fv

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

