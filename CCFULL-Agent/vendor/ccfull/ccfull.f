C
C     Program for coupled-channels equations for heavy-ion 
C     fusion reactions
C
C     The couplings are always treated to all orders.
C
C     Anharmonic couplings and generalised rotational couplings (but for E2 
C     only: E4 couplings are yet to be worked out) are included in this 
C     version of CCFULL. 
C
C     last modification: February 16th, 2015
c                        June 29th, 2010: closed channels
c                        December 14th, 2010: delete the option "P=1-R"
c                        February 16th, 2015: symmetrization
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 ANS
C
C  Parameters:
C    NLEVELMAX - maximum number of intrinsic levels to be included
C
      PARAMETER (NLEVELMAX=30,NRMAX=5000)
C
C Principal global variables
C
C    P         - penetrability 
C    SIGMA     - fusion cross section, unit mb
C    SPIN      - mean angular momentum
C
C    AP,ZP,RP  - atomic #, proton # and radius of the projectile
C    AT,ZT,RT  - those of the target
C    RMASS     - reduced mass, unit MeV/c**2
C    ANMASS    - nucleon mass, unit MeV/c**2
C    E         - bombarding energy in the center of mass frame, unit MeV
C    V0,R0,A0  - depth, range, and dissuseness parameters of uncoupled 
C                nuclear potential, which is assumed to be a Woods-Saxon form
C    IVIBROT   - option for intrinsic degree of freedom 
C                = -1; no excitation (inert)/ =0 ; vibrational coupling/ 
C                =1; rotational coupling
C    BETA      - defeormation parameter
C    OMEGA     - excitation energy of the oscillator
C    LAMBDA    - multipolarity
C    NPHONON   - the number of phonon to be included
C    BETA2     - static quadrupole deformation parameter
C    BETA4     - static hexadecapole deformation parameter
C    E2        - excitation energy of 2+ state in a rotational band
C    NROT      - the number of levels in the rotational band to be included 
C                (up to I^pi=2*NROT+ states are included)
C    L         - angular momentum of the relative motion
C    CPOT      - coupling matrix
C    

        dimension cpot0(nlevelmax,nlevelmax)

	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/oscp/betap,betapn,omegap,lambdap,nphononp
	common/rotp/beta2p,beta4p,e2p,nrotp
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/dimen/nlevel
        common/pot/v0,r0,a0,pow
	common/dyn/e,rmass
        common/angular/l
	common/shape/rb,vb,curv
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,0:nrmax)
        common/couph/cpoth(nlevelmax,nlevelmax)                  
        common/trans/ftr,qtrans,iq,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)
        common/eigen/ev(nlevelmax),evec(nlevelmax,nlevelmax)

C
C Output files: 
C   OUTPUT: detailed information of the calculation
C   cross.dat : fusion cross sections
C   spin.dat  : averaged angular momenta
C   partial.dat: partial cross sections (being calculated only when a single 
C                          value of the energy is entered in the input file)
C
	OPEN(7,FILE='cross.dat',STATUS='UNKNOWN')
	OPEN(8,FILE='spin.dat',STATUS='UNKNOWN')
	OPEN(9,FILE='OUTPUT',STATUS='UNKNOWN')
	OPEN(10,FILE='ccfull.inp',STATUS='UNKNOWN')
	open(15,file='s-wave.dat',status='unknown')

C
C Define two constants used in various places
C
      hbar=197.329d0
      pi=3.141592653d0
C
C Input phase
C
      anmass=938.d0

      ap=16.
      zp=8.
      at=148.
      zt=62.
	read(10,*)ap,zp,at,zt
C  r_coup, unit fm
      rp=1.2d0*ap**(1.d0/3.d0)
      rt=1.06d0*at**(1.d0/3.d0)
	read(10,*)rp,ivibrotp,rt,ivibrott
	rp=rp*ap**(1.d0/3.d0)
	rt=rt*at**(1.d0/3.d0)

      rmass=ap*at/(ap+at)*anmass

C-----
c      IVIBROTT=0
C             = -1 for inert, =0 for vibration, =1 for rotation
	if(ivibrott.eq.-1) then
	      ntarg=0
              read(10,*)omega,beta,alambda,nphonon
	elseif(ivibrott.eq.0) then
c (input for target phonon)
              omegat=0.55d0
              betat=0.182d0
              lambdat=2
              nphonont=3
              read(10,*)omegat,betat,lambdat,nphonont
	      ntarg=nphonont
	else
c (input for target rotation)
	      e2t=0.d0
	      beta2t=0.d0
	      beta4t=0.d0
              nrott=0
              lambdat=2
              read(10,*)e2t,beta2t,beta4t,nrott
	      ntarg=nrott
	endif
C-----
C (input for target phonon (the second mode of excitation))
              nphonont2=3
c                     =0; no second mode in the target
              omegat2=1.16d0
              betat2=0.236d0
              lambdat2=3
              read(10,*)omegat2,betat2,lambdat2,nphonont2
	      if(ivibrott.eq.-1) nphonont2=0
c-----
c      ivibrotp=-1
c             = -1 for inert, =0 for vibration, =1 for rotation
	if(ivibrotp.eq.-1) then
	      nproj=0
              read(10,*)omegap,betap,alambdap,nphononp
	else if(ivibrotp.eq.0) then
c (input for projectile phonon)
              omegap=3.74d0
              betap=0.339d0
              lambdap=3
              nphononp=1
              read(10,*)omegap,betap,lambdap,nphononp
              nproj=nphononp
	else
c (input for projectile rotation)
              e2p=0.d0
	      beta2p=0.d0
	      beta4p=0.d0
	      nrotp=0
              read(10,*)e2p,beta2p,beta4p,nrotp
	      nproj=nrotp
	endif

      ntrans=0
c ntrans= 0 ; no transfer channel/ =1 ; with transfer channel
c ftrans(r)=ftr*dvn(r)/dr
c
      ftr=7.d0
      qtrans=5.21d0
        	read(10,*)ntrans,qtrans,ftr

c--------------- parameters for the nuclear potential
              r1=1.233d0*ap**(1.d0/3.d0)-0.98d0*ap**(-1.d0/3.d0)
              r2=1.233d0*at**(1.d0/3.d0)-0.98d0*at**(-1.d0/3.d0)
              r0=r1+r2+0.29d0
              r12=r1*r2/(r1+r2)
              a0=0.63d0
              gamma=0.95d0*(1.d0-1.8d0*(ap-2.d0*zp)/ap*(at-2.d0*zt)/at)
              v0=16.d0*pi*gamma*r12*a0
              v0=31.67*r1*r2/(r1+r2)

                v0=1551d0
                r0=0.95d0*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
                a0=1.05d0
                        read(10,*)v0,r0,a0
	                r0=r0*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))

        pow=1.d0
	WRITE(6,*)' Modified Woods-Saxon (power of WS) (n/y)?'
	read(5,1)ans
 1	format(a1)
	     IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	     WRITE(6,*)'power=?'
	     read(5,*)pow
                if(pow.eq.0.d0) then
                WRITE(6,*)'Power cannot be zero :^('
                stop
                endif
	     endif

c  the minimum energy, the maximum energy, and the energy step

	emin=53.d0
	emax=72.d0
	de=2.d0*at/(ap+at)/6.d0
        	read(10,*)emin,emax,de

      rmax=50.d0
      dr=0.05d0
        	read(10,*)rmax,dr

c========================================== end of input phase
        lskip=1
        fsym=1.
        if(ap.eq.at.and.zp.eq.zt) then
           lskip=2
           fsym=2.
           endif

	if(ivibrott.ne.-1.and.ntarg.eq.0.and.nphonont2.ne.0) then
	write(6,91)
 91	FORMAT(1H ,'The target excitation should be input'
     &  /' in the first mode :^(')
	stop
	endif

c	if(ivibrott.eq.1.and.ivibrotp.eq.1) then
c	write(6,92)
c 92	format(1h ,'this code cannot handle fusion between'
c     &  / 'two deformed nuclei :^(')
c	stop
c	endif

c	if(ivibrott.eq.1.and.nphonont2.ne.0) then
c 	write(6,*)'!!!!!!!!!!!!!!! warning !!!!!!!!!!!!!!!!!!'
c	write(6,*)'    '
c	write(6,*)'the vibrational phonon in the target is not'
c        write(6,*)'a deformed phonon, but a *spherical* phonon,'
c        write(6,*)'i.e. the rotation-vibration model is not used.'
c	write(6,*)'    '
c 	write(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c	write(6,*)'    '
c	endif

         call nucleus
	 nlevel=nlevel+ntrans
	 if(nlevel.gt.nlevelmax) then
	 write(6,*)'too many channels :^('
	 stop
	 endif
 100	write(6,*) 'nlevel=',nlevel

C Evaluation coupling matrix elements
C------------------------------------------------------
      l=0
C   search a minimum of the Coulomb pocket RMIN
      call potshape
      ngrid=int((rmax-rmin)/dr)+1.d-6
             if(ngrid+1.gt.nrmax) then
	     write(6,*)'too large grid !!  :^('
	     stop
	     endif
      rmax=rmin+dr*ngrid

      call cmat0
      do ir=0,ngrid+1
      r=rmin+dr*ir
      call cmat(r,cpot0)
      do i=1,nlevel
      do j=1,nlevel
         cpot(i,j,ir)=cpot0(i,j)
      enddo
      enddo
      enddo
      
      rh=rmin+dr/2.d0                
      call cmat(rh,cpoth)            

	write(9,99)
	write(9,97)
 99	FORMAT(/,'       Ecm (MeV)    sigma (mb)          <l> ')
 97	format('       -------------------------------------')

        iwarn=0

c  start the energy loop
c-------------------------------------------------------------
	if(de.eq.0.d0) then
	iestep=0
	else
        iestep=(emax-emin)/de+1.d-6
	endif

	if(iestep.eq.0) then
	open(11,file='partial.dat',status='unknown')
	endif

      do 10 i=0,iestep
      e=emin+de*i

      sigma=0.d0
      spin=0.d0
c  angular momentum loop
c-----------------------
      do 20 l=0,200,lskip
      s0=sigma

      rmin0=rmin
      call potshape
      rmin00=rmin
      rmin=rmin0
      if(v(rmin00).gt.e.or.rmin00.lt.0.d0) then
      p=0.d0
        if(l.ne.0.and.iwarn.eq.0) then
        write(6,*)'    '
        WRITE(6,*)'!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!'
        WRITE(6,*)'Shallow potential!!!!'
        WRITE(6,*)'You may get an oscillatioin in BD at E > ', E
        WRITE(6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        WRITE(6,*)'    '
        WRITE(9,*)'    '
        WRITE(9,*)'!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!'
        WRITE(9,*)'Shallow potential!!!!'
        WRITE(9,*)'You may get an oscillatioin in BD at E > ', E
        WRITE(9,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(9,*)'    '
        iwarn=1
        endif
      go to 15
      endif

c
c integration of the c.c. equations
      call mnumerov(p)
      write(15,*)e,p

c	if(l.eq.0) write(9,*)e,p
	if(nlevel.gt.10) write(6,*)e,l*1.,p

      if(iestep.eq.0) then
         write(11,*)e,l*1.
     &          ,fsym*(2.d0*l+1.d0)*p*pi*hbar**2/2.d0/rmass/e*10.d0
      endif

      sigma=sigma+(2.d0*l+1.d0)*p*pi*hbar**2/2.d0/rmass/e*10.d0
      spin=spin+l*(2.d0*l+1.d0)*p*pi*hbar**2/2.d0/rmass/e*10.d0
      if(sigma-s0.lt.s0*1.d-4) goto 15
 20   continue

 15   if(sigma.ne.0.d0) spin=spin/sigma
      sigma=sigma*fsym

      write(7,*)e,sigma
      write(8,*)e,spin
       if(sigma.lt.1.d-2) then
       write(6,81)e,sigma,spin
       write(9,81)e,sigma,spin
       else
       write(6,82)e,sigma,spin
       write(9,82)e,sigma,spin
       endif

 10   continue
 81   format(f15.5,e15.5,f15.5)
 82   format(3f15.5)

      WRITE(6,*)'The program terminated normally :^)'	
      if(iwarn.eq.1) then
      write(6,*)'    '
      WRITE(6,*)'.....but you may get an oscillation in BD'
      WRITE(6,*)' (see OUTPUT file).'
      endif

      stop
      end

c*************************************************************
      subroutine nucleus
c
c subroutine to record several input parameters in the 
c 'output' file
c
c*************************************************************
      implicit real*8(a-h,o-z)
      parameter (nlevelmax=30)
      character*1 ans,sgn
      character*2 text,pro,targ
      dimension text(109)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/oscp/betap,betapn,omegap,lambdap,nphononp
	common/rotp/beta2p,beta4p,e2p,nrotp
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/phonon/nphonon
	common/dimen/nlevel
        common/pot/v0,r0,a0,pow
	common/dyn/e,rmass
	common/shape/rb,vb,curv
	common/grid/rmin,rmax,dr
        common/trans/ftr,qtrans,iq,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)
        common/grott/bett(0:nlevelmax,0:nlevelmax),erott(0:nlevelmax)
        common/grotp/betp(0:nlevelmax,0:nlevelmax),erotp(0:nlevelmax)

      DATA TEXT/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ',
     &'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     &'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     &'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     &'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re',
     &'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     &'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md',
     &'No','Lr','XX','X1','X2','X3','X4','04'/
 
        nlevel=(nproj+1)*(ntarg+1)
	np=zp
	npp=ap
	pro=text(np)
	nt=zt
	ntt=at
	targ=text(nt)
	WRITE(9,*)NPP,PRO,'  +  ',NTT,TARG, '     Fusion reaction'
	WRITE(6,*)NPP,PRO,'  +  ',NTT,TARG, '     Fusion reaction'

	write(9,30)
	write(6,30)

	r0p=rp/ap**(1.d0/3.d0)
	r0t=rt/at**(1.d0/3.d0)

	if(ntarg.ne.0) then
	if(ivibrott.eq.0) then
	write(6,21)betat,omegat,lambdat,nphonont 
 21	FORMAT('Phonon Excitation in the targ.: beta=',F6.3,
     &    ', omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	betatn=betat
	write(6,*)'    '
	WRITE(6,*)' Different beta_N from beta_C for this mode(n/y)?'
	READ(5,1)ANS
	     IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	     WRITE(6,*)'beta_N=?'
	     read(5,*)betatn
	     endif
	write(9,211)betatn,betat,r0t,omegat,lambdat,nphonont 
 211	FORMAT('Phonon Excitation in the targ.: beta_N=',F6.3,
     &    ', beta_C=',F6.3,', r0=',F5.2,'(fm),',
     &     /32X,'omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)
	endif

	if(ivibrott.eq.1) then
	write(9,22)beta2t,beta4t,r0t,e2t,nrott
	write(6,22)beta2t,beta4t,r0t,e2t,nrott
 22     FORMAT('Rotational Excitation in the targ.: beta2=',F6.3,
     &    ', beta4=',F6.3,', r0=',F5.2,'(fm),',
     &      /36X,'E2=',F5.2,'(MeV), Nrot=',I2)
	endif
	endif

	if(nphonont2.ne.0) then
	write(6,23)betat2,omegat2,lambdat2,nphonont2 
 23	FORMAT('Phonon Excitation in the targ.: beta=',F6.3,
     &    ', omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	betat2n=betat2
	write(6,*)'    '
	WRITE(6,*)' Different beta_N from beta_C for this mode(n/y)?'
	READ(5,1)ANS
	     IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	     WRITE(6,*)'beta_N=?'
	     read(5,*)betat2n
	     endif

	write(9,233)betat2n,betat2,r0t,omegat2,lambdat2,nphonont2 
 233	FORMAT('Phonon Excitation in the targ.: beta_N=',F6.3,
     &    ', beta_C=',F6.3,', r0=',F5.2,'(fm),',
     &    /32X, 'omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	call mutual
	endif

	if(nproj.ne.0) then
	if(ivibrotp.eq.0) then
	write(6,24)betap,omegap,lambdap,nphononp 
 24	FORMAT('Phonon Excitation in the proj.: beta=',F6.3,
     &    ', omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)

	betapn=betap
	write(6,*)'    '
	WRITE(6,*)' Different beta_N from beta_C for this mode(n/y)?'
	READ(5,1)ANS
	     IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	     WRITE(6,*)'beta_N=?'
	     read(5,*)betapn
	     endif
	write(9,244)betapn,betap,r0p,omegap,lambdap,nphononp 
 244	FORMAT('Phonon Excitation in the proj.: beta_N=',F6.3,
     &    ', beta_C=',F6.3,', r0=',F5.2,'(fm),',
     &    /32X, 'omega=',F5.2, '(MeV), Lambda=',I2,', Nph=',I2)
	endif

	if(ivibrotp.eq.1) then
	write(9,25)beta2p,beta4p,r0p,e2p,nrotp
	write(6,25)beta2p,beta4p,r0p,e2p,nrotp
 25     FORMAT('Rotational Excitation in the proj.: beta2=',F6.3,
     &    ', beta4=',F6.3,', r0=',F5.2,'(fm),',
     &      /36X,'E2=',F5.2,'(MeV), Nrot=',I2)
	endif
	endif

        call anharmonicity
        call grotation

	if(ntrans.ne.0) then
	write(9,30)
	write(9,26)ftr, qtrans
	write(6,26)ftr, qtrans
 26	FORMAT('Transfer channel: Strength=',F5.2,
     &    ', Q=',F5.2, '(MeV)')
        IQ=0
        WRITE(6,*)'     '
        WRITE(6,*)'Which type (charge change)?'
        WRITE(6,*)'         [0 for n-transfer, 1 for p-pick up,', 
     &              ' -1 for p-stripping, etc.]'
        read(5,*)iq
        if(iq.gt.0) then
        sgn='+'
        write(6,291)sgn,iq
        write(9,291)sgn,iq
        else
        write(6,29)iq
        write(9,29)iq
        endif
 29     FORMAT('   Charge difference of the transfer partition:',I5)
 291    FORMAT('   Charge difference of the transfer partition:  ',
     &          a,i2)
        write(6,*)'     '
	endif

	write(9,30)

	r00=r0/(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
	write(9,27)v0,r00,a0,pow
 27	FORMAT('Potential parameters: V0=', F8.2, '(MeV), r0=', 
     &        F5.2, '(fm), a=', F5.2, '(fm), power=',F5.2)
		CALL POTSHAPE
		WRITE(9,28)RB,VB,CURV
 28 	FORMAT('   Uncoupled barrier: Rb=', F5.2, '(fm), Vb=',F8.2, 
     &        '(MeV), Curv=',F5.2, '(MeV)')
	WRITE(9,30)
 30	FORMAT('-------------------------------------------------')

 1	FORMAT(A1)
      RETURN
      END

c*************************************************************
      subroutine anharmonicity
c
c anharmonic couplings for vibrations
c
c*************************************************************
      implicit real*8(a-h,o-z)
      parameter (nlevelmax=30)
      character*1 ans,sign
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/oscp/betap,betapn,omegap,lambdap,nphononp
	common/rotp/beta2p,beta4p,e2p,nrotp
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/phonon/nphonon
	common/dimen/nlevel
        common/pot/v0,r0,a0,pow
	common/dyn/e,rmass
	common/shape/rb,vb,curv
	common/grid/rmin,rmax,dr
        common/trans/ftr,qtrans,iq,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)

c-----------------------------
c ahv coupling for the first mode in the target
           do it=0,ntarg
           do jt=0,ntarg

           betnahv(it,jt)=0.d0
           betcahv(it,jt)=0.d0

           if(it.eq.jt-1) then
              betnahv(it,jt)=betatn*sqrt(jt*1.d0)
              betcahv(it,jt)=betat*sqrt(jt*1.d0)
           else if(jt.eq.it-1) then
              betnahv(it,jt)=betatn*sqrt(it*1.d0)
              betcahv(it,jt)=betat*sqrt(it*1.d0)
           endif

           enddo
              omeahv(it)=it*omegat
           enddo

        if(ivibrott.eq.0.and.ntarg.gt.1) then

	WRITE(6,90)'AHV couplings for the first mode ',
     &             'in the target phonon (n/y)?'
 90	format(2a)
	read(5,1)ans
        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(9,*)'            '
           WRITE(9,90)'**** AHV Couplings in the target',
     &                ' (the 1st mode) ****'
           write(9,*)'             '
	if((-1.d0)**lambdat.eq.1.d0) sign='+'
	if((-1.d0)**lambdat.eq.-1.d0) sign='-'

           do it=0,ntarg
           do 44 jt=0,ntarg
           if(it.gt.jt) goto 44
           
           if(it.eq.0.and.jt.eq.0) goto 44
           if(it.eq.0.and.jt.eq.1) goto 44

        write(6,*)'             '
	WRITE(6,94)'Transition from the',LAMBDAT,SIGN,'^',IT,
     &                  ' to the',LAMBDAT,SIGN,'^',JT,'state:'
        WRITE(6,95)'   beta_N and beta_C in the HO limit=',
     &             BETNAHV(IT,JT),BETCAHV(IT,JT)
        WRITE(6,90)'    Modify these beta_N a/o beta_C (n/y)?'
	READ(5,1)ANS

        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	WRITE(9,94)'Transition from the',LAMBDAT,SIGN,'^',IT,
     &                  ' to the',LAMBDAT,SIGN,'^',JT,'state:'
        WRITE(9,95)'   beta_N and beta_C in the HO limit=',
     &             BETNAHV(IT,JT),BETCAHV(IT,JT)
        WRITE(6,*)'   beta_N and beta_C =?'
        READ(5,*)BETNAHV(IT,JT),BETCAHV(IT,JT)
        WRITE(9,97)BETNAHV(IT,JT),BETCAHV(IT,JT)
        ENDIF

 94     format(1h ,a,i2,2a,i2,a,i2,2a,i2,a)
 95     format(1h ,a,2g15.5)
 97     format('    beta_n=',g15.5,',   beta_c=',g15.5)

 44        enddo
           enddo

              WRITE(6,*)'           '
              WRITE(6,*)'Excitation energy for the first mode:'
              DO I=2,NTARG
              WRITE(6,98)'    Energy of the',I,
     &              '-phonon state in the HO=',I*OMEGAT
              WRITE(6,*)'   Modify this energy(n/y)?'
              READ(5,1)ANS
              IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
              WRITE(6,*)'   Energy=?'
              READ(5,*)OMEAHV(I)
              WRITE(6,*)'       '
              WRITE(9,61)'Energy of the',I,'-phonon state=',OMEAHV(I),
     &                         '(HO: ',I*OMEGAT,')'
              endif
              enddo
              write(9,*)'          '              
              write(6,*)'          '              
 98           format(a,i2,a,g15.5)
 61           format(1h ,a,i2,a,g15.5,a,g12.5,a)
            endif
            endif
c-----------------------------
c ahv coupling for the second mode in the target
           do it=0,nphonont2
           do jt=0,nphonont2

           betnahv2(it,jt)=0.d0
           betcahv2(it,jt)=0.d0

           if(it.eq.jt-1) then
              betnahv2(it,jt)=betat2n*sqrt(jt*1.d0)
              betcahv2(it,jt)=betat2*sqrt(jt*1.d0)
           else if(jt.eq.it-1) then
              betnahv2(it,jt)=betat2n*sqrt(it*1.d0)
              betcahv2(it,jt)=betat2*sqrt(it*1.d0)
           endif

           enddo
              omeahv2(it)=it*omegat2
           enddo

        if(nphonont2.gt.1) then

	WRITE(6,90)'AHV couplings for the second mode ',
     &             'in the target phonon (n/y)?'
	READ(5,1)ANS
        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(9,*)'            '
           WRITE(9,90)'**** AHV Couplings in the target',
     &                ' (the 2nd mode) ****'
           write(9,*)'             '

	if((-1.d0)**lambdat2.eq.1.d0) sign='+'
	if((-1.d0)**lambdat2.eq.-1.d0) sign='-'

           do it=0,nphonont2
           do 54 jt=0,nphonont2
           if(it.gt.jt) goto 54
           
           if(it.eq.0.and.jt.eq.0) goto 54
           if(it.eq.0.and.jt.eq.1) goto 54

        write(6,*)'             '
	WRITE(6,94)'Transition from the',LAMBDAT2,SIGN,'^',IT,
     &                  ' to the',LAMBDAT2,SIGN,'^',JT,'state:'
        WRITE(6,95)'    beta_N and beta_C in the HO limit=',
     &             BETNAHV2(IT,JT),BETCAHV2(IT,JT)
        WRITE(6,90)'     Modify these beta_N a/o beta_C (n/y)?'
	READ(5,1)ANS

        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	WRITE(9,94)'Transition from the',LAMBDAT2,SIGN,'^',IT,
     &                  ' to the',LAMBDAT2,SIGN,'^',JT,'state:'
        WRITE(9,95)'   beta_N and beta_C in the HO limit=',
     &             BETNAHV(IT,JT),BETCAHV(IT,JT)
        WRITE(6,*)'    beta_N and beta_C =?'
        READ(5,*)BETNAHV2(IT,JT),BETCAHV2(IT,JT)
        WRITE(9,97)BETNAHV2(IT,JT),BETCAHV2(IT,JT)
        ENDIF

 54        ENDDO
           ENDDO

              WRITE(6,*)'           '
              WRITE(6,*)'Excitation energy for the second mode:'
              DO I=2,NPHONONT2
              WRITE(6,98)'    Energy of the',I,
     &              '-phonon state in the HO=',I*OMEGAT2
            WRITE(6,*)'   Modify this energy(n/y)?'
              READ(5,1)ANS
              IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
              WRITE(6,*)'   Energy=?'
              READ(5,*)OMEAHV2(I)
              WRITE(9,61)'Energy of the',I,'-phonon state=',OMEAHV2(I),
     &                         '(HO: ',I*OMEGAT2,')'
              ENDIF
              ENDDO
              WRITE(9,*)'          '              
              WRITE(6,*)'          '              
            ENDIF
            ENDIF

c-----------------------------
c ahv coupling for the projectile phonon 
           do it=0,nproj
           do jt=0,nproj

           betnahvp(it,jt)=0.d0
           betcahvp(it,jt)=0.d0

           if(it.eq.jt-1) then
              betnahvp(it,jt)=betapn*sqrt(jt*1.d0)
              betcahvp(it,jt)=betap*sqrt(jt*1.d0)
           else if(jt.eq.it-1) then
              betnahvp(it,jt)=betapn*sqrt(it*1.d0)
              betcahvp(it,jt)=betap*sqrt(it*1.d0)
           endif

           enddo
              omeahvp(it)=it*omegap
           enddo

        if(ivibrotp.eq.0.and.nproj.gt.1) then

	write(6,72)
 72     FORMAT('AHV couplings for the projectile phonon (n/y)?')
	READ(5,1)ANS
        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(9,*)'            '
           WRITE(9,73)
 73        FORMAT(1H ,'**** AHV Couplings in the projectile ****') 
           write(9,*)'            '
	if((-1.d0)**lambdap.eq.1.d0) sign='+'
	if((-1.d0)**lambdap.eq.-1.d0) sign='-'

           do it=0,nproj
           do 64 jt=0,nproj
           if(it.gt.jt) goto 64
           
           if(it.eq.0.and.jt.eq.0) goto 64
           if(it.eq.0.and.jt.eq.1) goto 64

        write(6,*)'             '
	WRITE(6,94)'Transition from the',LAMBDAP,SIGN,'^',IT,
     &                  ' to the',LAMBDAP,SIGN,'^',JT,'state:'
        WRITE(6,95)'  beta_N and beta_C in the HO limit=',
     &             BETNAHVP(IT,JT),BETCAHVP(IT,JT)
        WRITE(6,90)'   Modify these beta_N a/o beta_C (n/y)?'
	read(5,1)ans

        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	WRITE(9,94)'Transition from the',LAMBDAP,SIGN,'^',IT,
     &                  ' to the',LAMBDAP,SIGN,'^',JT,'state:'
        WRITE(9,95)'   beta_N and beta_C in the HO limit=',
     &             BETNAHVP(IT,JT),BETCAHVP(IT,JT)
        WRITE(6,*)'  beta_N and beta_C =?'
        READ(5,*)BETNAHVP(IT,JT),BETCAHVP(IT,JT)
        WRITE(9,97)BETNAHVP(IT,JT),BETCAHVP(IT,JT)
        ENDIF

 64        ENDDO
           ENDDO

              WRITE(6,*)'           '
           WRITE(6,*)'Excitation energy for the projectile phonon:'
           DO I=2,NPROJ
              WRITE(6,98)'    Energy of the',I,
     &              '-phonon state in the HO=',I*OMEGAP
            WRITE(6,*)'   Modify this energy(n/y)?'
              READ(5,1)ANS
              IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
              WRITE(6,*)'   Energy=?'
              READ(5,*)OMEAHVP(I)
              WRITE(9,61)'Energy of the',I,'-phonon state=',OMEAHVP(I),
     &                         '(HO: ',I*OMEGAP,')'
              ENDIF
              ENDDO
              WRITE(9,*)'          '              
              WRITE(6,*)'          '              
            ENDIF
            ENDIF

 1	FORMAT(A1)

      RETURN
      END

c*************************************************************
      subroutine grotation
c
c generalised e2 rotational couplings
c
c*************************************************************
      implicit real*8(a-h,o-z)
      parameter (nlevelmax=30)
      character*1 ans
	common/hion/ap,zp,rp,at,zt,rt
	common/intr/ivibrotp,ivibrott,nproj,ntarg
	common/rotp/beta2p,beta4p,e2p,nrotp
	common/rott/beta2t,beta4t,e2t,nrott
        common/grott/bett(0:nlevelmax,0:nlevelmax),erott(0:nlevelmax)
        common/grotp/betp(0:nlevelmax,0:nlevelmax),erotp(0:nlevelmax)

c-----------------------------
c the generalised couplings for the target rotor
           do it=0,ntarg
           do jt=0,ntarg
           bett(it,jt)=0.d0
           if(it.eq.jt-1.or.it.eq.jt+1.or.it.eq.jt) bett(it,jt)=beta2t
           enddo
           erott(it)=2.d0*it*(2.d0*it+1.d0)/6.d0*e2t
           enddo

        if(ivibrott.eq.1.and.ntarg.gt.1) then

        WRITE(6,*)'      '
	WRITE(6,*)' Generalised E2 couplings for the target rotor (n/y)?'
	READ(5,1)ANS
        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(9,*)'            '
           WRITE(9,90)'**** Generalised couplings for the ',
     &                ' target rotor ****'
           WRITE(9,*)'             '
           WRITE(6,*)'             '

           WRITE(6,*)'Beta2_targ=',BETA2T

           DO IT=1,NTARG-1
           JT=IT+1

           WRITE(6,95)'Modify beta2 for the transition from',
     &                2*IT,'+ to',2*JT,'+ ? (n/y)'
   	   READ(5,1)ANS

           IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(6,*)'   BETA2=?'           
           READ(5,*)BETT(IT,JT)
               BETT(JT,IT)=BETT(IT,JT)
           WRITE(9,94)'Beta2 for the transition from',2*IT,'+ to',
     &                2*JT,'+ =: ',BETT(IT,JT)
           ENDIF
           ENDDO

              WRITE(6,*)'           '
              WRITE(6,*)'Reorientation terms:'

           DO IT=1,NTARG
           JT=IT

           WRITE(6,95)'Modify beta2 for the transition from',
     &                2*IT,'+ to',2*JT,'+ ? (n/y)'
   	   READ(5,1)ANS

           IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(6,*)'   BETA2=?'           
           READ(5,*)BETT(IT,JT)
               BETT(JT,IT)=BETT(IT,JT)
           WRITE(9,94)'Beta2 for the transition from',2*IT,'+ to',
     &                2*JT,'+ =: ',BETT(IT,JT)
           ENDIF
           ENDDO


              WRITE(6,*)'           '
              WRITE(6,*)'Excitation energies:'
              DO I=2,NTARG
              WRITE(6,98)'    Energy of the',2*I,
     &              '+ state for the pure rotor=',EROTT(I)
              WRITE(6,*)'   Modify this energy(n/y)?'
              READ(5,1)ANS
              IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
                 EROTT0=EROTT(I)
              WRITE(6,*)'   Energy=?'
              READ(5,*)EROTT(I)
              WRITE(6,*)'       '
              WRITE(9,61)'Energy of the',2*I,'+ state=',EROTT(I),
     &                         '(rotor: ',EROTT0,')'
              ENDIF
              ENDDO
              WRITE(9,*)'          '              
              WRITE(6,*)'          '              
 98           FORMAT(A,I2,A,G15.5)
 61           FORMAT(1H ,A,I2,A,G15.5,A,G12.5,A)
            ENDIF
            ENDIF

 94     FORMAT(1H ,A,I2,A,I2,A,G15.5)
 95     FORMAT(1H ,A,I2,A,I2,A)
 90	FORMAT(2A)

C-----------------------------
c The generalised couplings for the projectile rotor
           DO IP=0,NPROJ
           DO JP=0,NPROJ
           BETP(IP,JP)=0.D0
           IF(IP.EQ.JP-1.OR.IP.EQ.JP+1.OR.IP.EQ.JP) BETP(IP,JP)=BETA2P
           ENDDO
           EROTP(IP)=2.D0*IP*(2.D0*IP+1.D0)/6.D0*E2P
           ENDDO

        IF(IVIBROTP.EQ.1.AND.NPROJ.GT.1) THEN

        WRITE(6,*)'      '
	WRITE(6,90)' Generalised E2 couplings for the projectile', 
     &             ' rotor (n/y)?'
	READ(5,1)ANS
        IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(9,*)'            '
           WRITE(9,90)'**** Generalised couplings for the ',
     &                ' projectile rotor ****'
           WRITE(9,*)'             '
           WRITE(6,*)'             '

           WRITE(6,*)'Beta2_proj=',BETA2P

           DO IP=1,NPROJ-1
           JP=IP+1

           WRITE(6,95)'Modify beta2 for the transition from',
     &                2*IP,'+ to',2*JP,'+ ? (n/y)'
   	   READ(5,1)ANS

           IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(6,*)'   BETA2=?'           
           READ(5,*)BETP(IP,JP)
               BETP(JP,IP)=BETP(IP,JP)
           WRITE(9,94)'Beta2 for the transition from',2*IP,'+ to',
     &                2*JP,'+ =: ',BETP(IP,JP)
           ENDIF
           ENDDO

              WRITE(6,*)'           '
              WRITE(6,*)'Reorientation terms:'

           DO IP=1,NPROJ
           JP=IP

           WRITE(6,95)'Modify beta2 for the transition from',
     &                2*IP,'+ to',2*JP,'+ ? (n/y)'
   	   READ(5,1)ANS

           IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
           WRITE(6,*)'   BETA2=?'           
           READ(5,*)BETP(IP,JP)
               BETP(JP,IP)=BETP(IP,JP)
           WRITE(9,94)'Beta2 for the transition from',2*IP,'+ to',
     &                2*JP,'+ =: ',BETP(IP,JP)
           ENDIF
           ENDDO


              WRITE(6,*)'           '
              WRITE(6,*)'Excitation energies:'
              DO I=2,NPROJ
              WRITE(6,98)'    Energy of the',2*I,
     &              '+ state for the pure rotor=',EROTP(I)
              WRITE(6,*)'   Modify this energy(n/y)?'
              READ(5,1)ANS
              IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
                 EROTP0=EROTP(I)
              WRITE(6,*)'   Energy=?'
              READ(5,*)EROTP(I)
              WRITE(6,*)'       '
              WRITE(9,61)'Energy of the',2*I,'+ state=',EROTP(I),
     &                         '(rotor: ',EROTP0,')'
              ENDIF
              ENDDO
              WRITE(9,*)'          '              
              WRITE(6,*)'          '              
            ENDIF
            ENDIF

 1	FORMAT(A1)

        RETURN
        END

c*****************************************************************************
      subroutine mutual
c
c subroutine to select levels to be included in the c.c. calculations
c
c*****************************************************************************
        implicit real*8(a-h,m,o-z)
        parameter (nlevelmax=30)
	character*1 ans,sign1,sign2
	common/intr/ivibrotp,ivibrott,nproj,ntarg
        common/osct/betat,betatn,omegat,lambdat,nphonont
        common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
	common/rott/beta2t,beta4t,e2t,nrott
	common/dimen/nlevel
        common/trans/ftr,qtrans,iq,ntrans
	common/mut/imutual(0:nlevelmax,0:nlevelmax)

	do 10 i=0,nlevelmax
	do 20 j=0,nlevelmax
	imutual(i,j)=0
 20	continue
 10	continue

	write(6,*)'    '
	WRITE(6,*)'Mutual excitations in the *target* nucleus'

	WRITE(6,*)'  Include the mutual excitations (y/n)?'
	READ(5,1)ANS
 1	FORMAT(A1)
	     IF(ANS.EQ.'n' .OR. ANS.EQ.'N') THEN
	     IMUT=0
	     NLEVEL=(NTARG+NPHONONT2+1)*(NPROJ+1)
c	     WRITE(9,99)
c 99	     FORMAT(1H ,'No mutual excitations in the target are included.')
	     DO 80 I=0,NTARG
	     DO 81 J=0,NPHONONT2
	     IF(I.EQ.0.OR.J.EQ.0) IMUTUAL(I,J)=1
 81          CONTINUE
 80	     CONTINUE
	     GOTO 70
	     ENDIF

	WRITE(6,*)'  All the possible mutual excitation channels (n/y)?'
	READ(5,1)ANS
	     IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
	     IMUT=1
	     NLEVEL=(NTARG+1)*(NPHONONT2+1)*(NPROJ+1)
c	     WRITE(9,90)
c 90	     FORMAT(1H ,'All the possible mutual excitation channels'
c     &                  1H ,'in the target are included.')
	     DO 82 I=0,NTARG
	     DO 83 J=0,NPHONONT2
	     IMUTUAL(I,J)=1
 83          CONTINUE
 82	     CONTINUE
	     GOTO 70
	     ENDIF

	IMUT=2
	IF((-1.D0)**LAMBDAT.EQ.1.D0) SIGN1='+'
	IF((-1.D0)**LAMBDAT.EQ.-1.D0) SIGN1='-'
	IF((-1.D0)**LAMBDAT2.EQ.1.D0) SIGN2='+'
	IF((-1.D0)**LAMBDAT2.EQ.-1.D0) SIGN2='-'
	DO 30 I=0,NTARG
	DO 40 J=0,NPHONONT2
		IF(I.EQ.0.OR.J.EQ.0) THEN 
		IMUTUAL(I,J)=1
		GOTO 40
		ENDIF
	WRITE(6,94)'  Include (',LAMBDAT,SIGN1,'^',I,',',LAMBDAT2,SIGN2,
     &            '^',J,') state ? (y/n)'
	READ(5,1)ANS
	     IF(ANS.EQ.'n' .OR. ANS.EQ.'N') THEN
	     IMUTUAL(I,J)=0
	     GOTO 40
	     ENDIF
	IMUTUAL(I,J)=1
 40	CONTINUE
 30	CONTINUE

 70	WRITE(9,*)'    '
	WRITE(9,93)
	WRITE(6,93)
 93	FORMAT(1H ,'  Excited states in the target to be included: ')

	NLEVEL=0
	IF((-1.D0)**LAMBDAT.EQ.1.D0) SIGN1='+'
	IF((-1.D0)**LAMBDAT.EQ.-1.D0) SIGN1='-'
	IF((-1.D0)**LAMBDAT2.EQ.1.D0) SIGN2='+'
	IF((-1.D0)**LAMBDAT2.EQ.-1.D0) SIGN2='-'
	DO 50 I=0,NTARG
	DO 60 J=0,NPHONONT2
		IF(IMUTUAL(I,J).EQ.0) GOTO 60
	WRITE(9,94)'    (',LAMBDAT,SIGN1,'^',I,',',LAMBDAT2,SIGN2,
     &                   '^',J,') state'
	WRITE(6,94)'    (',LAMBDAT,SIGN1,'^',I,',',LAMBDAT2,SIGN2,
     &                   '^',J,') state'
 94	FORMAT(1H ,A,I2,2A,I2,A,I2,2A,I2,A)
	NLEVEL=NLEVEL+1
 60	CONTINUE
 50	CONTINUE
        NLEVEL=NLEVEL*(NPROJ+1)

	RETURN
	END
c*****************************************************************************
      subroutine potshape
c
c shape of the bare coulomb barrier, i.e. the barrier position, the curvature,
c and the position of the coulomb pocket
c
c*****************************************************************************
      implicit real*8(a-h,m,o-z)
      common/shape/rb,vb,curv
      common/grid/rmin,rmax,dr
      common/dyn/e,rmass
      common/const/hbar,pi
      external dv,v

        r=50.5d0
        u0=dv(r)
 10	r=r-1.d0

           if(r.lt.0.d0) then
              rb=-5.d0
	      rmin=-5.d0
              return
              endif

	u1=dv(r)
	if(u0*u1.gt.0.d0) goto 10
		ra=r+1.d0
		rb=r
	        tolk=1.d-6
                n=log10(abs(rb-ra)/tolk)/log10(2.d0)+0.5d0
                do 20 i=1,n
                r=(ra+rb)/2.d0
                u=dv(r)
                if(u0*u.lt.0.d0) then
                rb=r
                else
                ra=r
                endif
 20             continue

	rb=r

        ddv00=(dv(rb+1.d-5)-dv(rb-1.d-5))/2.d-5

	if(ddv00.gt.0.d0) then
	write(6,*)'something is strange :^('
	write(8,*)'something is strange :^('
	stop
	endif

	curv=hbar*sqrt(abs(ddv00)/rmass)
	vb=v(rb)

		ra=rb-0.5d0
		u0=dv(ra)
		rbb=0.5d0
		u1=dv(rbb)
	        tolk=1.d-6
                n=log10(abs(rb-ra)/tolk)/log10(2.d0)+0.5d0
                do 21 i=1,n
                r=(ra+rbb)/2.d0
                u=dv(r)
                if(u0*u.lt.0.d0) then
                rbb=r
                else
                ra=r
                endif
 21             continue

	rmin=r

c	write(6,*)'rmin=',rmin,v(rmin)
c       write(6,*)'rb=',rb,v(rb)
c       write(6,*)'curv=',curv00

	return
	end

c**************************************************************
      subroutine mnumerov00(p)
c
c subroutine for integration of the c.c. eqs. by modified 
c numerov method
c
c**************************************************************
      implicit real*8 (a-h,o-z)
      external v
      parameter (nlevelmax=30,nrmax=5000)

      dimension psi(nlevelmax),psi0(nlevelmax),psi1(nlevelmax)
      dimension xi(nlevelmax),xi0(nlevelmax),xi1(nlevelmax)
      dimension phi0(nlevelmax)
      dimension bb(nlevelmax,nlevelmax),
     &          bin(nlevelmax,nlevelmax)
      dimension bb2(nlevelmax,nlevelmax)
      dimension cc(nlevelmax,nlevelmax),
     &          cin(nlevelmax,nlevelmax)
      dimension dd0(nlevelmax,nlevelmax),dd1(nlevelmax,nlevelmax)
      dimension dd(nlevelmax)
      dimension fcw(0:200),gcw(0:200),fpcw(0:200),gpcw(0:200)
      dimension sigmad(0:200),iexp(0:200)
      dimension ech(nlevelmax),ech2(nlevelmax)
 
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,0:nrmax)

      complex*16 k,kk,k2
      complex*16 psi,psi0,psi1,phi0
      complex*16 xi,xi0,xi1
      complex*16 ai
      complex*16 cwup0,cwdown0,cwup1,cwdown1
      complex*16 bb,bb2,bin,cc,cin,dd
      complex*16 dummy,dummy2

      ai=(0.d0,1.d0)
c define two constants used in this subroutine
      fac=dr**2*(2.d0*rmass/hbar**2) 
      ngrid=int((rmax-rmin)/dr)
      rmax=rmin+ngrid*dr

	do 22 lc=0,200
	fcw(lc)=0.d0
	gcw(lc)=0.d0
	fpcw(lc)=0.d0
	gpcw(lc)=0.d0
	sigmad(lc)=0.d0
	iexp(lc)=0d0
 22	continue

      do 51 i=1,nlevelmax
      do 52 j=1,nlevelmax
      bb(i,j)=0.d0
      bin(i,j)=0.d0
      cc(i,j)=0.d0
      cin(i,j)=0.d0
 52   continue
      dd(i)=0.d0
 51   continue

c integration of the io-th channel wave function from rmin 
c==========================================================
      do 15 io=1,nlevel
c	write(6,*)'io=',io*1.

      do 200 j1=1,nlevel
      psi(j1)=0.d0
      psi0(j1)=0.d0
      psi1(j1)=0.d0
      phi0(j1)=0.d0
      xi0(j1)=0.d0
      xi1(j1)=0.d0
      xi(j1)=0.d0
 200  continue

        if(io.eq.1) then
        do io2=1,nlevel
        ech(io2)=e-v(rmin)-cpot(io2,io2,0)
        ech2(io2)=e-cpot(io2,io2,ngrid)
        enddo
        endif

c  initial conditions
c----------
	if(ech(io).gt.0.d0) then
        k=sqrt(2.d0*rmass/hbar**2*ech(io))
        psi0(io)=exp(-ai*k*rmin)
        phi0(io)=-ai*k*psi0(io)
	else
        k=sqrt(2.d0*rmass/hbar**2*abs(ech(io)))
        psi0(io)=exp(k*rmin)
        phi0(io)=k*psi0(io)
	endif

      call rkutta00(psi0,phi0,psi1)

      do 91 i0=1,nlevel
	xi0(i0)=(1.d0-fac/12.d0*(v(rmin)-e))*psi0(i0)
	xi1(i0)=(1.d0-fac/12.d0*(v(rmin+dr)-e))*psi1(i0)
      do 92 ic=1,nlevel
      xi0(i0)=xi0(i0)-fac/12.d0*cpot(i0,ic,0)*psi0(ic)
      xi1(i0)=xi1(i0)-fac/12.d0*cpot(i0,ic,1)*psi1(ic)
 92   continue
 91   continue

      do 29 ir=2,ngrid+1
c----------------------------------------------  iterations start
      r=rmin+dr*ir
      r0=rmin+dr*(ir-2)
      r1=rmin+dr*(ir-1)

      do 71 i0=1,nlevel
      do 72 ic=1,nlevel
         dd0(i0,ic)=fac/sqrt(12.d0)*cpot(i0,ic,ir-1)
      if(i0.eq.ic) dd0(i0,ic)=dd0(i0,ic)+fac/sqrt(12.d0)*(v(r1)-e)
     &                        +sqrt(3.d0)
 72   continue
 71   continue

      do 73 i0=1,nlevel
      do 74 ic=1,nlevel
      dd1(i0,ic)=0.d0
      if(i0.eq.ic) dd1(i0,ic)=dd1(i0,ic)-1.d0
      do 75 ik=1,nlevel
         dd1(i0,ic)=dd1(i0,ic)+dd0(i0,ik)*dd0(ik,ic)
 75   continue
 74   continue
 73   continue
	
      do 54 i0=1,nlevel
         xi(i0)=-xi0(i0)
      do 53 ic=1,nlevel
         xi(i0)=xi(i0)+dd1(i0,ic)*xi1(ic)
 53   continue
 54   continue

      if(ir.eq.ngrid+1) goto 66
      do 42 i0=1,nlevel
         xi0(i0)=xi1(i0)
         xi1(i0)=xi(i0)
 42   continue
 29   continue
c--------------------------------------------------------------
c  matching to the coulomb wave function at rmax

 66   do 56 i0=1,nlevel
      do 55 ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*cpot(i0,ic,ngrid-1)
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)-fac/12.d0*(v(rmax-dr)-e)+1.d0
 55   continue
 56   continue
      call matinv(nlevel,cc,cin)
      do 94 i0=1,nlevel
         psi0(i0)=0.d0
      do 93 ic=1,nlevel
         psi0(i0)=psi0(i0)+cin(i0,ic)*xi0(ic)
 93   continue
 94   continue

      do 556 i0=1,nlevel
      do 555 ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*cpot(i0,ic,ngrid+1)
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)-fac/12.d0*(v(rmax+dr)-e)+1.d0
 555   continue
 556   continue
      call matinv(nlevel,cc,cin)
      do 944 i0=1,nlevel
         psi(i0)=0.d0
      do 933 ic=1,nlevel
         psi(i0)=psi(i0)+cin(i0,ic)*xi(ic)
 933   continue
 944   continue

      do 60 ii=1,nlevel
c   coulomb wave function
      ec=e-cpot(ii,ii,ngrid-1)
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax-dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup0=(gcw(l)+ai*fcw(l))
      cwdown0=gcw(l)-ai*fcw(l)

      ec=e-cpot(ii,ii,ngrid+1)
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax+dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup1=(gcw(l)+ai*fcw(l))
      cwdown1=gcw(l)-ai*fcw(l)

      bb(io,ii)=(cwup0*psi(ii)-cwup1*psi0(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
      bb2(io,ii)=(cwdown1*psi0(ii)-cwdown0*psi(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
 60   continue
c---------------------------------------------------------------
 15   continue

c===============================================================
c                                        penetration probability
c  calculate the inversion of the matrix bb
      call matinv(nlevel,bb,bin)
c
      p=0.d0
      ref=0.d0

      do 85 io=1,nlevel
      if(ech(io).lt.0.d0) goto 85
      k=sqrt((2.d0*rmass/hbar**2*ech(io)))
      k2=sqrt((2.d0*rmass/hbar**2*ech2(io)))
      kk=sqrt(2.d0*rmass/hbar**2*e)
      
             dummy=0.d0
             dummy2=0.d0
             do io2=1,nlevel
                dummy2=dummy2+bin(1,io2)*bb2(io2,io)
                enddo

      p=p+(abs(bin(1,io)))**2*k/kk
      ref=ref+(abs(dummy2))**2*k2/kk
 85   continue

c      if(p+ref.gt.1.1d0) then
c         write(6,*)'penetrability > 1.1 !!!!!'
c         write(6,*)'something must be strange (accuracy problem) :^('
c         stop
c      endif

      if(p.gt.0.99d0) p=1.d0-ref

      return
      end

c**************************************************************
      subroutine mnumerov0(p)
c
c subroutine for integration of the c.c. eqs. by modified 
c numerov method
c
c a similar stabilization method as that given in 
c   z.h. levine, pra30('84) 1120
c   t.n. rescigno and a.e. orel, pra25('82)2402
c is used in the routine. 
c
c the linear independence condition is imposed at the barrier position.
c
c**************************************************************
      implicit real*8 (a-h,o-z)
      external v
      parameter (nlevelmax=30,istab0=5,nrmax=5000)

      dimension psi(nlevelmax),psi0(nlevelmax),psi1(nlevelmax)
      dimension xi(nlevelmax,nlevelmax),xi0(nlevelmax,nlevelmax),
     &          xi1(nlevelmax,nlevelmax)
      dimension phi0(nlevelmax)
      dimension bb(nlevelmax,nlevelmax),
     &          bin(nlevelmax,nlevelmax)
      dimension bb2(nlevelmax,nlevelmax)
      dimension cc(nlevelmax,nlevelmax),
     &          cin(nlevelmax,nlevelmax)
      dimension dd0(nlevelmax,nlevelmax),dd1(nlevelmax,nlevelmax)
      dimension dd(nlevelmax)
      dimension fcw(0:200),gcw(0:200),fpcw(0:200),gpcw(0:200)
      dimension sigmad(0:200),iexp(0:200)
      dimension ech(nlevelmax),ech2(nlevelmax)
 
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,0:nrmax)
	common/shape/rb,vb,curv
        
      complex*16 k,kk,k2
      complex*16 psi,psi0,psi1,phi0
      complex*16 xi,xi0,xi1
      complex*16 ai
      complex*16 cwup0,cwdown0,cwup1,cwdown1
      complex*16 bb,bb2,bin,cc,cin,dd
      complex*16 dummy,dummy2
      complex*16 xi1d(nlevelmax,nlevelmax)
      complex*16 xi0d(nlevelmax,nlevelmax)
      complex*16 aa(nlevelmax,nlevelmax)
      complex*16 bb0(nlevelmax,nlevelmax)
      complex*16 bb20(nlevelmax,nlevelmax)

      ai=(0.d0,1.d0)
c define two constants used in this subroutine
      fac=dr**2*(2.d0*rmass/hbar**2) 
      ngrid=int((rmax-rmin)/dr)
      rmax=rmin+ngrid*dr

      ibarrier=(rb-rmin)/dr+1.d-6

	do 22 lc=0,200
	fcw(lc)=0.d0
	gcw(lc)=0.d0
	fpcw(lc)=0.d0
	gpcw(lc)=0.d0
	sigmad(lc)=0.d0
	iexp(lc)=0d0
 22	continue

      do 51 i=1,nlevelmax
      do 52 j=1,nlevelmax
      bb(i,j)=0.d0
      bin(i,j)=0.d0
      cc(i,j)=0.d0
      cin(i,j)=0.d0
      aa(i,j)=0.d0
 52   continue
      dd(i)=0.d0
 51   continue

      do i=1,nlevel
         aa(i,i)=1.d0
         enddo

c  initial conditions
c==========================================================
      do 15 io=1,nlevel
c	write(6,*)'io=',io*1.

      do 200 j1=1,nlevel
      psi(j1)=0.d0
      psi0(j1)=0.d0
      psi1(j1)=0.d0
      phi0(j1)=0.d0
 200  continue

        if(io.eq.1) then
        do io2=1,nlevel
        ech(io2)=e-v(rmin)-cpot(io2,io2,0)
        ech2(io2)=e-cpot(io2,io2,ngrid)
        enddo
        endif
        
	if(ech(io).gt.0.d0) then
        k=sqrt(2.d0*rmass/hbar**2*ech(io))
        psi0(io)=exp(-ai*k*rmin)
        phi0(io)=-ai*k*psi0(io)
	else
        k=sqrt(2.d0*rmass/hbar**2*abs(ech(io)))
        psi0(io)=exp(k*rmin)
        phi0(io)=k*psi0(io)
        endif

      call rkutta00(psi0,phi0,psi1)

      do 91 i0=1,nlevel
	xi0(i0,io)=(1.d0-fac/12.d0*(v(rmin)-e))*psi0(i0)
	xi1(i0,io)=(1.d0-fac/12.d0*(v(rmin+dr)-e))*psi1(i0)
      do 92 ic=1,nlevel
      xi0(i0,io)=xi0(i0,io)-fac/12.d0*cpot(i0,ic,0)*psi0(ic)
      xi1(i0,io)=xi1(i0,io)-fac/12.d0*cpot(i0,ic,1)*psi1(ic)
 92   continue
 91   continue
 15   continue

c===============
      istab=0
      do 29 ir=2,ngrid+1
c----------------------------------------------  iterations start
      r=rmin+dr*ir
      r0=rmin+dr*(ir-2)
      r1=rmin+dr*(ir-1)
      istab=istab+1

      do 71 i0=1,nlevel
      do 72 ic=1,nlevel
         dd0(i0,ic)=fac/sqrt(12.d0)*cpot(i0,ic,ir-1)
      if(i0.eq.ic) dd0(i0,ic)=dd0(i0,ic)+fac/sqrt(12.d0)*(v(r1)-e)
     &                        +sqrt(3.d0)
 72   continue
 71   continue

      do 73 i0=1,nlevel
      do 74 ic=1,nlevel
      dd1(i0,ic)=0.d0
      if(i0.eq.ic) dd1(i0,ic)=dd1(i0,ic)-1.d0
      do 75 ik=1,nlevel
         dd1(i0,ic)=dd1(i0,ic)+dd0(i0,ik)*dd0(ik,ic)
 75   continue
 74   continue
 73   continue
	
      do 16 ich=1,nlevel
      do 54 i0=1,nlevel
         xi(i0,ich)=-xi0(i0,ich)
      do 53 ic=1,nlevel
         xi(i0,ich)=xi(i0,ich)+dd1(i0,ic)*xi1(ic,ich)
 53   continue
 54   continue
 16   continue

      if(ir.eq.ngrid+1) goto 66
      if(r.lt.12.d0.and.istab.eq.istab0) then
        call stabilize(xi1,xi,aa,ir)
        istab=0
        endif
c      if(ir.eq.ibarrier) call stabilize(xi1,xi,aa,ir)

      do ich=1,nlevel
      do i0=1,nlevel
         xi0(i0,ich)=xi1(i0,ich)
         xi1(i0,ich)=xi(i0,ich)
      enddo
      enddo

 29   continue
c--------------------------------------------------------------
c  matching to the coulomb wave function at rmax

 66   do 17 io=1,nlevel

      do 56 i0=1,nlevel
      do 55 ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*cpot(i0,ic,ngrid-1)
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)-fac/12.d0*(v(rmax-dr)-e)+1.d0
 55   continue
 56   continue
      call matinv(nlevel,cc,cin)
      do 94 i0=1,nlevel
         psi0(i0)=0.d0
      do 93 ic=1,nlevel
         psi0(i0)=psi0(i0)+cin(i0,ic)*xi0(ic,io)
 93   continue
 94   continue

      do 556 i0=1,nlevel
      do 555 ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*cpot(i0,ic,ngrid+1)
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)-fac/12.d0*(v(rmax+dr)-e)+1.d0
 555   continue
 556   continue
      call matinv(nlevel,cc,cin)
      do 944 i0=1,nlevel
         psi(i0)=0.d0
      do 933 ic=1,nlevel
         psi(i0)=psi(i0)+cin(i0,ic)*xi(ic,io)
 933   continue
 944   continue

      do 60 ii=1,nlevel
c   coulomb wave function
      ec=e-cpot(ii,ii,ngrid-1)
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax-dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup0=(gcw(l)+ai*fcw(l))
      cwdown0=gcw(l)-ai*fcw(l)

      ec=e-cpot(ii,ii,ngrid+1)
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax+dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup1=(gcw(l)+ai*fcw(l))
      cwdown1=gcw(l)-ai*fcw(l)

      bb0(ii,io)=(cwup0*psi(ii)-cwup1*psi0(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
      bb20(ii,io)=(cwdown1*psi0(ii)-cwdown0*psi(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
 60   continue
c---------------------------------------------------------------
 17   continue

c===============================================================
c                                        penetration probability
      do i=1,nlevel
      do j=1,nlevel
         bb(i,j)=0.d0
         bb2(i,j)=0.d0
         do j2=1,nlevel
            bb(i,j)=bb(i,j)+bb0(i,j2)*aa(j2,j)
            bb2(i,j)=bb2(i,j)+bb20(i,j2)*aa(j2,j)
            enddo
            enddo
            enddo

c  calculate the inversion of the matrix bb
      call matinv(nlevel,bb,bin)
c
      p=0.d0
      ref=0.d0

      do 85 io=1,nlevel
      if(ech(io).lt.0.d0) goto 85
      k=sqrt((2.d0*rmass/hbar**2*ech(io)))
      k2=sqrt((2.d0*rmass/hbar**2*ech2(io)))
      kk=sqrt(2.d0*rmass/hbar**2*e)
      
             dummy=0.d0
             dummy2=0.d0
             do io2=1,nlevel
                dummy2=dummy2+bb2(io,io2)*bin(io2,1)
                enddo

      p=p+(abs(bin(io,1)))**2*k/kk
      ref=ref+(abs(dummy2))**2*k2/kk
 85   continue

c      if(p+ref.gt.1.1d0) then
c         write(6,*)'penetrability > 1.1 !!!!!'
c         write(6,*)'something must be strange (accuracy problem) :^('
c         stop
c      endif

      if(p.gt.0.99d0) p=1.d0-ref

      return
      end

c**************************************************************
      subroutine mnumerov(p)
c
c subroutine for integration of the c.c. eqs. by modified 
c numerov method
c
c a similar stabilization method as that given in 
c   z.h. levine, pra30('84) 1120
c   t.n. rescigno and a.e. orel, pra25('82)2402
c is used in the routine. 
c
c the linear independence condition is imposed at the barrier position.
c
c**************************************************************
      implicit real*8 (a-h,o-z)
      external v
      parameter (nlevelmax=30,nrmax=5000)

      dimension psi(nlevelmax),psi0(nlevelmax),psi1(nlevelmax)
      dimension xi(nlevelmax,nlevelmax),xi0(nlevelmax,nlevelmax),
     &          xi1(nlevelmax,nlevelmax)
      dimension phi0(nlevelmax)
      dimension bb(nlevelmax,nlevelmax),
     &          bin(nlevelmax,nlevelmax)
      dimension bb2(nlevelmax,nlevelmax)
      dimension cc(nlevelmax,nlevelmax),
     &          cin(nlevelmax,nlevelmax)
      dimension dd0(nlevelmax,nlevelmax),dd1(nlevelmax,nlevelmax)
      dimension dd(nlevelmax)
      dimension fcw(0:200),gcw(0:200),fpcw(0:200),gpcw(0:200)
      dimension sigmad(0:200),iexp(0:200)
      dimension ech(nlevelmax),ech2(nlevelmax)
 
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,0:nrmax)
	common/shape/rb,vb,curv
        
      complex*16 k,kk,k2
      complex*16 psi,psi0,psi1,phi0
      complex*16 xi,xi0,xi1
      complex*16 ai
      complex*16 cwup0,cwdown0,cwup1,cwdown1
      complex*16 bb,bb2,bin,cc,cin,dd
      complex*16 dummy,dummy2
      complex*16 xi1d(nlevelmax,nlevelmax)
      complex*16 xi0d(nlevelmax,nlevelmax)
      complex*16 aa(nlevelmax,nlevelmax)
      complex*16 bb0(nlevelmax,nlevelmax)
      complex*16 bb20(nlevelmax,nlevelmax)

      ai=(0.d0,1.d0)
c define two constants used in this subroutine
      fac=dr**2*(2.d0*rmass/hbar**2) 
      ngrid=int((rmax-rmin)/dr)
      rmax=rmin+ngrid*dr

      ibarrier=(rb-rmin)/dr+1.d-6

	do lc=0,200
	fcw(lc)=0.d0
	gcw(lc)=0.d0
	fpcw(lc)=0.d0
	gpcw(lc)=0.d0
	sigmad(lc)=0.d0
	iexp(lc)=0d0
        enddo

      do i=1,nlevelmax
      do j=1,nlevelmax
      bb(i,j)=0.d0
      bin(i,j)=0.d0
      cc(i,j)=0.d0
      cin(i,j)=0.d0
      aa(i,j)=0.d0
      enddo
      dd(i)=0.d0
      enddo

      do i=1,nlevel
         aa(i,i)=1.d0
         enddo

c  initial conditions
c==========================================================
      do 15 io=1,nlevel
c	write(6,*)'io=',io*1.

      do j1=1,nlevel
      psi(j1)=0.d0
      psi0(j1)=0.d0
      psi1(j1)=0.d0
      phi0(j1)=0.d0
      enddo

        if(io.eq.1) then
        do io2=1,nlevel
        ech(io2)=e-v(rmin)-cpot(io2,io2,0)
        ech2(io2)=e-cpot(io2,io2,ngrid)
        enddo
        endif
        
	if(ech(io).gt.0.d0) then
        k=sqrt(2.d0*rmass/hbar**2*ech(io))
        psi0(io)=exp(-ai*k*rmin)
        phi0(io)=-ai*k*psi0(io)
	else
        k=sqrt(2.d0*rmass/hbar**2*abs(ech(io)))
        psi0(io)=exp(k*rmin)
        phi0(io)=k*psi0(io)
        endif

      call rkutta00(psi0,phi0,psi1)

      do i0=1,nlevel
	xi0(i0,io)=(1.d0-fac/12.d0*(v(rmin)-e))*psi0(i0)
	xi1(i0,io)=(1.d0-fac/12.d0*(v(rmin+dr)-e))*psi1(i0)
      do ic=1,nlevel
      xi0(i0,io)=xi0(i0,io)-fac/12.d0*cpot(i0,ic,0)*psi0(ic)
      xi1(i0,io)=xi1(i0,io)-fac/12.d0*cpot(i0,ic,1)*psi1(ic)
      enddo
      enddo

 15   continue

c===============
      do 29 ir=2,ngrid+1
c----------------------------------------------  iterations start
      r=rmin+dr*ir
      r0=rmin+dr*(ir-2)
      r1=rmin+dr*(ir-1)

      do i0=1,nlevel
      do ic=1,nlevel
         dd0(i0,ic)=fac/sqrt(12.d0)*cpot(i0,ic,ir-1)
      if(i0.eq.ic) dd0(i0,ic)=dd0(i0,ic)+fac/sqrt(12.d0)*(v(r1)-e)
     &                        +sqrt(3.d0)
      enddo
      enddo

      do i0=1,nlevel
      do ic=1,nlevel
      dd1(i0,ic)=0.d0
      if(i0.eq.ic) dd1(i0,ic)=dd1(i0,ic)-1.d0
      do ik=1,nlevel
         dd1(i0,ic)=dd1(i0,ic)+dd0(i0,ik)*dd0(ik,ic)
         enddo
         enddo
         enddo
	
      do ich=1,nlevel
      do i0=1,nlevel
         xi(i0,ich)=-xi0(i0,ich)
      do ic=1,nlevel
         xi(i0,ich)=xi(i0,ich)+dd1(i0,ic)*xi1(ic,ich)
         enddo
         enddo
         enddo

      if(ir.eq.ngrid+1) goto 66
      if(ir.eq.ibarrier) call stabilize(xi1,xi,aa,ir)

      do ich=1,nlevel
      do i0=1,nlevel
         xi0(i0,ich)=xi1(i0,ich)
         xi1(i0,ich)=xi(i0,ich)
      enddo
      enddo

 29   continue
c--------------------------------------------------------------
c  matching to the coulomb wave function at rmax

 66   do 17 io=1,nlevel

      do i0=1,nlevel
      do ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*cpot(i0,ic,ngrid-1)
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)-fac/12.d0*(v(rmax-dr)-e)+1.d0
      enddo
      enddo

      call matinv(nlevel,cc,cin)
      do i0=1,nlevel
         psi0(i0)=0.d0
      do ic=1,nlevel
         psi0(i0)=psi0(i0)+cin(i0,ic)*xi0(ic,io)
      enddo
      enddo

      do i0=1,nlevel
      do ic=1,nlevel
      cc(i0,ic)=-fac/12.d0*cpot(i0,ic,ngrid+1)
      if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)-fac/12.d0*(v(rmax+dr)-e)+1.d0
      enddo
      enddo

      call matinv(nlevel,cc,cin)
      do i0=1,nlevel
         psi(i0)=0.d0
      do ic=1,nlevel
         psi(i0)=psi(i0)+cin(i0,ic)*xi(ic,io)
      enddo
      enddo

      do 60 ii=1,nlevel
c   coulomb wave function
      ec=e-cpot(ii,ii,ngrid-1)

      if(ec.lt.0.) then
         r1=rmax-dr
         r2=rmax+dr
         ak=sqrt(2.d0*rmass*abs(ec)/hbar**2)
         
         bb0(ii,io)=(exp(-ak*r2)*psi0(ii)-exp(-ak*r1)*psi(ii))
     &               /(exp(ak*(r1-r2))-exp(-ak*(r1-r2)))
         bb20(ii,io)=-(exp(ak*r2)*psi0(ii)-exp(ak*r1)*psi(ii))
     &               /(exp(ak*(r1-r2))-exp(-ak*(r1-r2)))
         else
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax-dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup0=(gcw(l)+ai*fcw(l))
      cwdown0=gcw(l)-ai*fcw(l)

      ec=e-cpot(ii,ii,ngrid+1)
      rho=sqrt(2.d0*rmass*ec)/hbar*(rmax+dr)
      eta=(zp*zt/137.d0)*sqrt(rmass/2.d0/ec)
      call dfcoul(eta,rho,fcw,fpcw,gcw,gpcw,sigmad,l,iexp)
      cwup1=(gcw(l)+ai*fcw(l))
      cwdown1=gcw(l)-ai*fcw(l)

      bb0(ii,io)=(cwup0*psi(ii)-cwup1*psi0(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
      bb20(ii,io)=(cwdown1*psi0(ii)-cwdown0*psi(ii))
     &                         /(cwup0*cwdown1-cwup1*cwdown0)
      endif

 60   continue
c---------------------------------------------------------------
 17   continue

c===============================================================
c                                        penetration probability
      do i=1,nlevel
      do j=1,nlevel
         bb(i,j)=0.d0
         bb2(i,j)=0.d0
         do j2=1,nlevel
            bb(i,j)=bb(i,j)+bb0(i,j2)*aa(j2,j)
            bb2(i,j)=bb2(i,j)+bb20(i,j2)*aa(j2,j)
            enddo
            enddo
            enddo

c  calculate the inversion of the matrix bb
      call matinv(nlevel,bb,bin)
c
      p=0.d0
      ref=0.d0

      do 85 io=1,nlevel
      if(ech(io).lt.0.d0) goto 85
      k=sqrt((2.d0*rmass/hbar**2*ech(io)))
      kk=sqrt(2.d0*rmass/hbar**2*e)
      p=p+(abs(bin(io,1)))**2*k/kk
 85   continue

      do 86 io=1,nlevel
      if(ech2(io).lt.0.d0) goto 86
      k2=sqrt((2.d0*rmass/hbar**2*ech2(io)))
      kk=sqrt(2.d0*rmass/hbar**2*e)
             dummy=0.d0
             dummy2=0.d0
             do io2=1,nlevel
                dummy2=dummy2+bb2(io,io2)*bin(io2,1)
                enddo
      ref=ref+(abs(dummy2))**2*k2/kk
 86   continue

c      if(p+ref.gt.1.1d0) then
c         write(6,*)'penetrability > 1.1 !!!!!'
c         write(6,*)'something must be strange (accuracy problem) :^('
c         stop
c      endif

c      if(p.gt.0.99d0) p=1.d0-ref

      return
      end

c**************************************************************
      subroutine stabilize(xi1,xi,aa,ir1)
c
c a subroutine to stabilize the solutions by imposing the linear 
c independence. 
c
c refs.
c   z.h. levine, pra30('84) 1120
c   t.n. rescigno and a.e. orel, pra25('82)2402
c
c**************************************************************
      implicit real*8 (a-h,o-z)
      parameter (nlevelmax=30,nrmax=5000)

      dimension psi(nlevelmax,nlevelmax)
      dimension xi(nlevelmax,nlevelmax),xi1(nlevelmax,nlevelmax)
      dimension cc(nlevelmax,nlevelmax),
     &          cin(nlevelmax,nlevelmax)

	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,0:nrmax)

      complex*16 psi,xi,xi1
      complex*16 ai
      complex*16 cc,cin
      complex*16 aa(nlevelmax,nlevelmax)
      complex*16 aa0(nlevelmax,nlevelmax)
      complex*16 xid(nlevelmax,nlevelmax),xid1(nlevelmax,nlevelmax)

      do i=1,nlevel
      do j=1,nlevel
         aa0(i,j)=aa(i,j)
         xid(i,j)=xi(i,j)
         xid1(i,j)=xi1(i,j)
         enddo
         enddo

      ai=(0.d0,1.d0)
      fac=dr**2*(2.d0*rmass/hbar**2) 

c     transforme to the physical wave functions psi from xi 
      r=rmin+ir1*dr
         do i0=1,nlevel
         do ic=1,nlevel
          cc(i0,ic)=-fac/12.d0*cpot(i0,ic,ir1)
          if(i0.eq.ic) cc(i0,ic)=cc(i0,ic)-fac/12.d0*(v(r)-e)+1.d0
          enddo
          enddo
        call matinv(nlevel,cc,cin)
      do ich=1,nlevel
        do i0=1,nlevel
          psi(i0,ich)=0.d0
          do ic=1,nlevel
           psi(i0,ich)=psi(i0,ich)+cin(i0,ic)*xi(ic,ich)
           enddo
           enddo
           enddo

         do i=1,nlevel
         do j=1,nlevel
              aa(i,j)=0.d0
              do k=1,nlevel
                 aa(i,j)=aa(i,j)+psi(i,k)*aa0(k,j)
                 enddo
                 enddo
                 enddo

        call matinv(nlevel,psi,cin)
        
        do i=1,nlevel
        do j=1,nlevel
           xi(i,j)=0.d0
           xi1(i,j)=0.d0
           do k=1,nlevel
              xi(i,j)=xi(i,j)+xid(i,k)*cin(k,j)
              xi1(i,j)=xi1(i,j)+xid1(i,k)*cin(k,j)
              enddo
              enddo
              enddo

      return
      end

c**************************************************************
      subroutine rkutta00(psi0,phi0,psi1)
c
c wave functions at r=rmin+dr
c
c**************************************************************
      implicit real*8 (a-h,o-z)
      external v
      parameter (nlevelmax=30,nrmax=5000)

      dimension psi0(nlevelmax),phi0(nlevelmax),psi1(nlevelmax)
      dimension ak1(nlevelmax),ak2(nlevelmax),
     &          ak3(nlevelmax),ak4(nlevelmax)
      dimension bk1(nlevelmax),bk2(nlevelmax),
     &          bk3(nlevelmax),bk4(nlevelmax)
 
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/dimen/nlevel
	common/dyn/e,rmass
        common/angular/l
	common/grid/rmin,rmax,dr
        common/coup/cpot(nlevelmax,nlevelmax,0:nrmax)
        common/couph/cpoth(nlevelmax,nlevelmax)                  

      complex*16 psi0,phi0,psi1
      complex*16 ai
      complex*16 ak1,ak2,ak3,ak4
      complex*16 bk1,bk2,bk3,bk4

      ai=(0.d0,1.d0)
c define two constants used in this subroutine
      fac=dr*(2.d0*rmass/hbar**2) 

      r=rmin
      rpp=rmin+dr
      rh=rmin+dr/2.d0

      do 19 j1=1,nlevel
      ak1(j1)=0.d0
      ak2(j1)=0.d0
      ak3(j1)=0.d0
      ak4(j1)=0.d0
      bk1(j1)=0.d0
      bk2(j1)=0.d0
      bk3(j1)=0.d0
      bk4(j1)=0.d0
 19   continue

      do 58 i0=1,nlevel
      do 57 ic=1,nlevel
      ak1(i0)=ak1(i0)+fac*cpot(i0,ic,0)*psi0(ic)
 57   continue
      ak1(i0)=ak1(i0)-fac*(e-v(r))*psi0(i0)
      bk1(i0)=dr*phi0(i0)
 58   continue

      do 56 i0=1,nlevel
      do 55 ic=1,nlevel
      ak2(i0)=ak2(i0)+fac*cpoth(i0,ic)*(psi0(ic)+1.d0/2.d0*bk1(ic))
 55   continue
      ak2(i0)=ak2(i0)-fac*(e-v(rh))*(psi0(i0)+1.d0/2.d0*bk1(i0))
      bk2(i0)=dr*(phi0(i0)+1./2.*ak1(i0))
 56   continue

      do 54 i0=1,nlevel
      do 53 ic=1,nlevel
      ak3(i0)=ak3(i0)+fac*cpoth(i0,ic)*(psi0(ic)+1.d0/2.d0*bk2(ic))
 53   continue
      ak3(i0)=ak3(i0)-fac*(e-v(rh))*(psi0(i0)+1.d0/2.d0*bk2(i0))
      bk3(i0)=dr*(phi0(i0)+1.d0/2.d0*ak2(i0))
 54   continue

      do 42 i0=1,nlevel
      do 41 ic=1,nlevel
      ak4(i0)=ak4(i0)+fac*cpot(i0,ic,1)*(psi0(ic)+bk3(ic))
 41   continue
      ak4(i0)=ak4(i0)-fac*(e-v(rpp))*(psi0(i0)+bk3(i0))
      bk4(i0)=dr*(phi0(i0)+ak3(i0))
 42   continue

c wave functions at rmin+dr
      do 59 is=1,nlevel
      psi1(is)=psi0(is)+1.d0/6.d0*(bk1(is)+2.d0*bk2(is)
     &        +2.d0*bk3(is)+bk4(is))  
 59   continue
      return
      end

c******************************************************************
	subroutine cmat0
c
c preparation for the nuclear coupling matrix
c
c*****************************************************************
      implicit real*8 (a-h,o-z)
      parameter (nlevelmax=30)
      dimension a(nlevelmax,nlevelmax)
      common/hion/ap,zp,rp,at,zt,rt
      common/const/hbar,pi
      common/intr/ivibrotp,ivibrott,nproj,ntarg
      common/oscp/betap,betapn,omegap,lambdap,nphononp
      common/osct/betat,betatn,omegat,lambdat,nphonont
      common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
      common/rotp/beta2p,beta4p,e2p,nrotp
      common/rott/beta2t,beta4t,e2t,nrott
      common/dimen/nlevel
      common/mut/imutual(0:nlevelmax,0:nlevelmax)
      common/eigen/ev(nlevelmax),evec(nlevelmax,nlevelmax)
      common/trans/ftr,qtrans,iq,ntrans
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)
        common/grott/bett(0:nlevelmax,0:nlevelmax),erott(0:nlevelmax)
        common/grotp/betp(0:nlevelmax,0:nlevelmax),erotp(0:nlevelmax)

      external cg

        fvibt=rt*betatn/sqrt(4.d0*pi)
	frot2t=rt*beta2t
	frot4t=rt*beta4t
        fvibt2=rt*betat2n/sqrt(4.d0*pi)
        fvibp=rp*betapn/sqrt(4.d0*pi)
	frot2p=rp*beta2p
	frot4p=rp*beta4p

	do 10 i=1,nlevelmax
	do 20 j=1,nlevelmax
	a(i,j)=0.d0
 20	continue
 10	continue

	ndim=nlevel-ntrans
	if(ndim.eq.1) return
	
	   i=0
           do 31 ip=0,nproj
           do 32 it=0,ntarg
           do 33 it2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(it,it2).eq.0) goto 33
			endif
	   i=i+1
           j=0
           do 41 jp=0,nproj
           do 42 jt=0,ntarg
           do 43 jt2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(jt,jt2).eq.0) goto 43
			endif
	   j=j+1

              if(i.gt.j) then
                 a(i,j)=a(j,i)
                 goto 43
                 endif

	   c=0.d0
	   if(ip.eq.jp.and.it2.eq.jt2) then
             if(ivibrott.eq.0) then
c               if(it.eq.jt-1) c=sqrt(jt*1.d0)*fvibt
c               if(jt.eq.it-1) c=sqrt(it*1.d0)*fvibt

               c=rt*betnahv(it,jt)/sqrt(4.d0*pi)

	     else
c		c=frot2t*sqrt((2.d0*2*it+1.d0)*5.d0*(2.d0*2*jt+1.d0)
c     &    /4.d0/pi)*cg(2*it,0,2,0,2*jt,0)**2/(2.d0*2*jt+1.d0)
		c=rt*bett(it,jt)
     &             *sqrt((2.d0*2*it+1.d0)*5.d0*(2.d0*2*jt+1.d0)
     &             /4.d0/pi)*cg(2*it,0,2,0,2*jt,0)**2/(2.d0*2*jt+1.d0)
		c=c+frot4t*sqrt((2.d0*2*it+1.d0)*9.d0*(2.d0*2*jt+1.d0)
     &             /4.d0/pi)*cg(2*it,0,4,0,2*jt,0)**2/(2.d0*2*jt+1.d0)
            endif
	   endif
	   if(ip.eq.jp.and.it.eq.jt) then
c               if(it2.eq.jt2-1) c=c+sqrt(jt2*1.d0)*fvibt2
c               if(jt2.eq.it2-1) c=c+sqrt(it2*1.d0)*fvibt2

               c=c+rt*betnahv2(it2,jt2)/sqrt(4.d0*pi)

	   endif

c-------
	   if(it.eq.jt.and.it2.eq.jt2) then
             if(ivibrotp.eq.0) then
c               if(ip.eq.jp-1) c=c+sqrt(jp*1.d0)*fvibp
c               if(jp.eq.ip-1) c=c+sqrt(ip*1.d0)*fvibp

               c=c+rp*betnahvp(ip,jp)/sqrt(4.d0*pi)

	     else
c		c=c+frot2p*sqrt((2.d0*2*ip+1.d0)*5.d0*(2.d0*2*jp+1.d0)
c     &     /4.d0/pi)*cg(2*ip,0,2,0,2*jp,0)**2/(2.d0*2*jp+1.d0)
		c=c+rp*betp(ip,jp)
     &              *sqrt((2.d0*2*ip+1.d0)*5.d0*(2.d0*2*jp+1.d0)
     &              /4.d0/pi)*cg(2*ip,0,2,0,2*jp,0)**2/(2.d0*2*jp+1.d0)
		c=c+frot4p*sqrt((2.d0*2*ip+1.d0)*9.d0*(2.d0*2*jp+1.d0)
     &              /4.d0/pi)*cg(2*ip,0,4,0,2*jp,0)**2/(2.d0*2*jp+1.d0)
             endif
            endif


	a(i,j)=c

 43              continue
 42              continue
 41              continue
 33              continue
 32              continue
 31              continue

	call mdiag(a,ndim)

      return
      end
c******************************************************************
	subroutine cmat(r,cpot)
c
c coupling matrix
c
c*****************************************************************
      implicit real*8 (a-h,o-z)
      parameter (nlevelmax=30)
      dimension cpot(nlevelmax,nlevelmax)
      common/hion/ap,zp,rp,at,zt,rt
      common/const/hbar,pi
      common/intr/ivibrotp,ivibrott,nproj,ntarg
      common/oscp/betap,betapn,omegap,lambdap,nphononp
      common/osct/betat,betatn,omegat,lambdat,nphonont
      common/osct2/betat2,betat2n,omegat2,lambdat2,nphonont2
      common/rotp/beta2p,beta4p,e2p,nrotp
      common/rott/beta2t,beta4t,e2t,nrott
      common/dimen/nlevel
      common/mut/imutual(0:nlevelmax,0:nlevelmax)
      common/eigen/ev(nlevelmax),evec(nlevelmax,nlevelmax)
      common/trans/ftr,qtrans,iq,ntrans
        common/ahv/betnahv(0:nlevelmax,0:nlevelmax),
     &             betcahv(0:nlevelmax,0:nlevelmax),
     &             omeahv(0:nlevelmax)
        common/ahv2/betnahv2(0:nlevelmax,0:nlevelmax),
     &             betcahv2(0:nlevelmax,0:nlevelmax),
     &             omeahv2(0:nlevelmax)
        common/ahvp/betnahvp(0:nlevelmax,0:nlevelmax),
     &             betcahvp(0:nlevelmax,0:nlevelmax),
     &             omeahvp(0:nlevelmax)
        common/grott/bett(0:nlevelmax,0:nlevelmax),erott(0:nlevelmax)
        common/grotp/betp(0:nlevelmax,0:nlevelmax),erotp(0:nlevelmax)
      external fct,fcp,fact,vncc,cg,fct2,fct4,fcp2,fcp4,fctt
      external fct2v,fcp2v
      external ftrans

	do 10 i=1,nlevelmax
	do 20 j=1,nlevelmax
	cpot(i,j)=0.d0
 20	continue
 10	continue

	ndim=nlevel-ntrans
	if(nlevel.eq.1) return
	if(ndim.eq.1) goto 55
	
	   i=0
           do 31 ip=0,nproj
           do 32 it=0,ntarg
           do 33 it2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(it,it2).eq.0) goto 33
			endif
	   i=i+1
           j=0
           do 41 jp=0,nproj
           do 42 jt=0,ntarg
           do 43 jt2=0,nphonont2
			if(nphonont2.ne.0) then
			if(imutual(jt,jt2).eq.0) goto 43
			endif
	   j=j+1

              if(i.gt.j) then
                 cpot(i,j)=cpot(j,i)
                 goto 43
                 endif
	c=0.d0
c nuclear coupling

	do 50 k=1,ndim
	c=c+vncc(r,ev(k))*evec(i,k)*evec(j,k)
 50	continue

c coulomb coupling 

	   if(ip.eq.jp.and.it2.eq.jt2) then
             if(ivibrott.eq.0) then
c               if(it.eq.jt-1) c=c+sqrt(jt*1.d0)*fct(r)
c               if(jt.eq.it-1) c=c+sqrt(it*1.d0)*fct(r)

               if(it.ne.jt) then
               cc=betcahv(it,jt)*fct(r)
               if(betat.ne.0.d0) cc=cc/betat
               else
               cc=betcahv(it,jt)*fct2v(r)/sqrt(4.d0*pi)
               endif

               c=c+cc

	     else

        	bet=bett(it,jt)+2.d0*sqrt(5.d0/pi)*bett(it,jt)**2/7.d0
c		c=c+sqrt((2.d0*2*it+1.d0)*5.d0*(2.d0*2*jt+1.d0)
c     &     /4.d0/pi)*cg(2*it,0,2,0,2*jt,0)**2/(2.d0*2*jt+1.d0)*fct2(r)
		c=c+sqrt((2.d0*2*it+1.d0)*5.d0*(2.d0*2*jt+1.d0)
     &     /4.d0/pi)*cg(2*it,0,2,0,2*jt,0)**2/(2.d0*2*jt+1.d0)
     &              *fct2(r)*bet
     &		   +sqrt((2.d0*2*it+1.d0)*9.d0*(2.d0*2*jt+1.d0)
     &     /4.d0/pi)*cg(2*it,0,4,0,2*jt,0)**2/(2.d0*2*jt+1.d0)*fct4(r)
            endif
	   endif

	   if(it.eq.jt.and.it2.eq.jt2) then
             if(ivibrotp.eq.0) then
c               if(ip.eq.jp-1) c=c+sqrt(jp*1.d0)*fcp(r)
c               if(jp.eq.ip-1) c=c+sqrt(ip*1.d0)*fcp(r)

               if(ip.ne.jp) then
               cc=betcahvp(ip,jp)*fcp(r)
               if(betap.ne.0.d0) cc=cc/betap
               else
               cc=betcahvp(ip,jp)*fcp2v(r)/sqrt(4.d0*pi)
               endif
               c=c+cc

	     else
        	bet=betp(ip,jp)+2.d0*sqrt(5.d0/pi)*betp(ip,jp)**2/7.d0
c		c=c+sqrt((2.d0*2*ip+1.d0)*5.d0*(2.d0*2*jp+1.d0)
c     &     /4.d0/pi)*cg(2*ip,0,2,0,2*jp,0)**2/(2.d0*2*jp+1.d0)*fcp2(r)
		c=c+sqrt((2.d0*2*ip+1.d0)*5.d0*(2.d0*2*jp+1.d0)
     &     /4.d0/pi)*cg(2*ip,0,2,0,2*jp,0)**2/(2.d0*2*jp+1.d0)
     &          *fcp2(r)*bet
     &            +sqrt((2.d0*2*ip+1.d0)*9.d0*(2.d0*2*jp+1.d0)
     &     /4.d0/pi)*cg(2*ip,0,4,0,2*jp,0)**2/(2.d0*2*jp+1.d0)*fcp4(r)
             endif
            endif

	   if(ip.eq.jp.and.it.eq.jt) then
c               if(it2.eq.jt2-1) c=c+sqrt(jt2*1.d0)*fctt(r)
c               if(jt2.eq.it2-1) c=c+sqrt(it2*1.d0)*fctt(r)

               if(it2.ne.jt2) then
               cc=betcahv2(it2,jt2)*fctt(r)
               if(betat2.ne.0.d0) cc=cc/betat2
               else
               cc=betcahv2(it2,jt2)*fct2v(r)/sqrt(4.d0*pi)
               endif
               c=c+cc

	   endif

c excitation energy
           if(it.eq.jt.and.ip.eq.jp.and.it2.eq.jt2) then
		if(ivibrott.eq.0) then
c		c=c+it*omegat
		c=c+omeahv(it)
		else
c		c=c+2.d0*it*(2.d0*it+1.d0)/6.d0*e2t
		c=c+erott(it)
		endif
		if(ivibrotp.eq.0) then
c		c=c+ip*omegap
		c=c+omeahvp(ip)
		else
c		c=c+2.d0*ip*(2.d0*ip+1.d0)/6.d0*e2p
		c=c+erotp(ip)
		endif
c	        c=c+it2*omegat2
	        c=c+omeahv2(it2)
 	   endif

	cpot(i,j)=c

 43              continue
 42              continue
 41              continue
 33              continue
 32              continue
 31              continue

      c0=cpot(1,1)
      do 51 i=1,ndim
          cpot(i,i)=cpot(i,i)-vn(r)
c         cpot(i,i)=cpot(i,i)-c0
 51      continue

c    transfer channel
 55      if(ntrans.eq.1) then
            cpot(1,nlevel)=ftrans(r)
            cpot(nlevel,1)=ftrans(r)
	    cpot(nlevel,nlevel)=cpot(1,1)-qtrans
     &          +(zp+iq)*(zt-iq)/r*hbar/137.d0-zp*zt/r*hbar/137.d0
            endif

      return
      end

c*****************************************************************
	function v(r)
c
c  potential for the relative motion
c
c*****************************************************************
	implicit real*8(a-h,o-z)
	external vc,vn,vcent
	v=vc(r)+vn(r)+vcent(r)
	return
	end
c*****************************************************************
	function vc(r)
c
c coulomb potential
c
c*****************************************************************
	implicit real*8(a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
c  coulomb radius
c	rc=1.2d0*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
c	if(r.gt.rc) then
	vc=zp*zt/r*hbar/137.d0
c	else
c	vc=zp*zt*hbar/137.d0*(3.d0*rc**2-r**2)/2.d0/rc**3
c	endif
	return
	end
c*****************************************************************
	function dvc(r)
c
c first derivative of the coulomb potential
c
c*****************************************************************
	implicit real*8(a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
c  coulomb radius
c	rc=1.2d0*(ap**(1.d0/3.d0)+at**(1.d0/3.d0))
c	if(r.gt.rc) then
	dvc=-zp*zt/r**2*hbar/137.d0
c	else
c	vc=zp*zt*hbar/137.d0*(3.d0*rc**2-r**2)/2.d0/rc**3
c	endif
	return
	end
c****************************************************************
	function vn(r)
c
c (power of) woods-saxon potential for nuclear potential
c
c****************************************************************
	implicit real*8 (a-h,o-z)
	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/pot/v0,r0,a0,pow

	if(r.lt.50.d0) then
        vn=-v0/(1.d0+exp((r-r0)/a0/pow))**pow
	else
	vn=0.d0
        endif

	return
	end
c****************************************************************
	function vncc(r,xt)
c
c (power of) woods-saxon potential for nuclear potential 
c
c****************************************************************
	implicit real*8 (a-h,o-z)

	common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/pot/v0,r0,a0,pow

	if(r.lt.50.d0) then
           rr=r-r0-xt
        vncc=-v0/(1.d0+exp(rr/a0/pow))**pow
	else
	vncc=0.d0
        endif

	return
	end

c***************************************************************
      function dvn(r)
c     
c the first derivative of vn(r)
c
c**************************************************************
      implicit real*8(a-h,o-z)
      common/pot/v0,r0,a,pow
      external vn

      if(r.lt.50.d0) then
c      dvn=v0/a*exp((r-r0)/a)/(1.d0+exp((r-r0)/a))**2
      dvn=(vn(r+1.d-5)-vn(r-1.d-5))/2.d-5
      else
      dvn=0.d0
      endif

      return
      end
c************************************************************
	function vcent(r)
c
c centrifugal potential
c
c************************************************************
	implicit real*8(a-h,o-z)
	common/const/hbar,pi
	common/angular/l
	common/dyn/e,rmass
	vcent=l*(l+1.d0)*hbar**2/2.d0/rmass/r**2
	return
	end
c************************************************************
	function dvcent(r)
c
c the first derivative of the centrifugal potential
c
c************************************************************
	implicit real*8(a-h,o-z)
	common/const/hbar,pi
	common/angular/l
	common/dyn/e,rmass
	dvcent=-2.d0*l*(l+1.d0)*hbar**2/2.d0/rmass/r**3
	return
	end

c***************************************************************
	function dv(r)
c
c the first derivative of v(r)
c
c***************************************************************
	implicit real*8(a-h,o-z)
	external dvn,dvc,dvcent

	dv=dvn(r)+dvc(r)+dvcent(r)

	return
	end
c*****************************************************************
	function fcp(r)
c
c coulomb coupling form factor for projectile phonon excitation
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/oscp/betap,betapn,omega,lambda,nphonon

	rc=rp
	if(r.gt.rc) then
	fcp=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

        fcp=fcp*betap/sqrt(4.d0*pi)

	return
	end

c*****************************************************************
	function fcp2v(r)
c
c coulomb coupling form factor for projectile excitation
c (for anharmonic vibration)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi

	lambda=2
	rc=rp
	if(r.gt.rc) then
	fcp2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

	return
	end
c*****************************************************************
	function fcp2(r)
c
c coulomb coupling form factor for projectile excitation
c (rotational e2 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rotp/beta2p,beta4p,e2p,nrotp

	lambda=2
	rc=rp
	if(r.gt.rc) then
	fcp2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

c	fcp2=fcp2*(beta2p+2.d0*sqrt(5.d0/pi)*beta2p**2/7.d0)

	return
	end
c*****************************************************************
	function fcp4(r)
c
c coulomb coupling form factor for projectile excitation
c (rotational e4 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rotp/beta2p,beta4p,e2p,nrotp

	lambda=4
	rc=rp
	if(r.gt.rc) then
	fcp4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rp/r)**lambda
	else
	fcp4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rp)**lambda
	endif

	fcp4=fcp4*(beta4p+9.d0*beta2p**2/7.d0/sqrt(pi))

	return
	end
c*****************************************************************
	function fct(r)
c
c coulomb coupling form factor for the target phonon excitation
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/osct/betat,betatn,omega,lambda,nphonon

	rc=rt
	if(r.gt.rc) then
	fct=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

        fct=fct*betat/sqrt(4.d0*pi)

	return
	end

c*****************************************************************
	function fctt(r)
c
c coulomb coupling form factor for the second target phonon excitation
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
        common/osct2/betat2,betat2n,omega2,lambda2,nphonon2

	rc=rt
	if(r.gt.rc) then
	fctt=3.d0/(2.d0*lambda2+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda2
	else
	fctt=3.d0/(2.d0*lambda2+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda2
	endif

        fctt=fctt*betat2/sqrt(4.d0*pi)

	return
	end
c*****************************************************************
	function fct2v(r)
c
c coulomb coupling form factor for target excitation
c (for anharmonic vibration)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi

	lambda=2
	rc=rt
	if(r.gt.rc) then
	fct2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct2v=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

	return
	end


c*****************************************************************
	function fct2(r)
c
c coulomb coupling form factor for target excitation
c (rotational e2 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rott/beta2t,beta4t,e2t,nrott

	lambda=2
	rc=rt
	if(r.gt.rc) then
	fct2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct2=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

c	fct2=fct2*(beta2t+2.d0*sqrt(5.d0/pi)*beta2t**2/7.d0)

	return
	end

c*****************************************************************
	function fct4(r)
c
c coulomb coupling form factor for target excitation
c (rotational e4 coupling)
c
c*****************************************************************
	implicit real*8(a-h,o-z)
        common/hion/ap,zp,rp,at,zt,rt
	common/const/hbar,pi
	common/rott/beta2t,beta4t,e2t,nrott

	lambda=4
	rc=rt
	if(r.gt.rc) then
	fct4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/r
     &                                         *(rt/r)**lambda
	else
	fct4=3.d0/(2.d0*lambda+1.d0)*zp*zt/137.d0*hbar/rc
     &                                         *(r/rt)**lambda
	endif

	fct4=fct4*(beta4t+9.d0*beta2t**2/7.d0/sqrt(pi))

	return
	end

c***************************************************************
      function ddvn(r)
c     
c the second derivative of vn(r)
c
c**************************************************************
      implicit real*8(a-h,o-z)
      common/pot/v0,r0,a,pow
      external dvn

      if(r.lt.50.d0) then
c      ddvn=v0/a**2*exp((r-r0)/a)/(1.d0+exp((r-r0)/a))**2
c     &    -2.d0*v0/a**2*exp(2.d0*(r-r0)/a)
c     &                           /(1.d0+exp((r-r0)/a))**3
      ddvn=(dvn(r+1.d-5)-dvn(r-1.d-5))/2.d-5
      else
      ddvn=0.d0
      endif

      return
      end
c****************************************************************
	function ftrans(r)
c
c coupling form factor for transfer reactions
c
c***************************************************************
	implicit real*8(a-h,o-z)
	common/const/hbar,pi
        common/trans/ftr,qtrans,iq,ntrans
	external dvn

	ftrans=ftr*dvn(r)
	return
	end

c*************************************************************
	subroutine mdiag(a,n)
c
c the jacobi method to compute all eigenvalues and eigenvectors 
c of a real symmetric matrix a, which is of size n by n, stored 
c in a physical np by np array. on output, the elements of a 
c above the diagonal are destroyed. d returns the eigenvalues of 
c a in its first n elements. v is a matrix with the same logical 
c and physical dimensions as a whose columns contain, on output, 
c the normalized eigenvectors of a. nrot returns the number of 
c jacobi rotations which were required. 
c
c n.b.; the i-th component of the eigenvector for the k-th 
c       eigenvalue is given by v(i,k). 
c
c************************************************************
	implicit real*8(a-h,o-z)
	parameter(nmax=100)
        parameter (nlevelmax=30)
        dimension a(nlevelmax,nlevelmax)
	dimension b(nmax),z(nmax)
        common/eigen/d(nlevelmax),v(nlevelmax,nlevelmax)

	do 12 ip=1,n
	do 11 iq=1,n
	v(ip,iq)=0.d0
 11	continue
	v(ip,ip)=1.d0
 12	continue

	do 13 ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.d0
 13	continue

	nrot=0

	do 24 i=1,50
	sm=0.d0

	do 15 ip=1,n-1
	do 14 iq=ip+1,n
	sm=sm+abs(a(ip,iq))
 14	continue
 15     continue

	if(sm.eq.0.d0) return

	if(i.lt.4) then
	tresh=0.2d0*sm/n**2
	else
	tresh=0.d0
	endif

	do 22 ip=1,n-1
	do 21 iq=ip+1,n
	g=100.d0*abs(a(ip,iq))
	
	if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     &     .and.(abs(d(iq))+g.eq.abs(d(iq)))) then
	   a(ip,iq)=0.d0
	elseif(abs(a(ip,iq)).gt.tresh) then
	   h=d(iq)-d(ip)

	   if(abs(h)+g.eq.abs(h)) then
	     t=a(ip,iq)/h
	   else
	     theta=0.5d0*h/a(ip,iq)
	     t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
	     if(theta.lt.0.d0) t=-t
	   endif

	c=1.d0/sqrt(1.d0+t**2)
	s=t*c
	tau=s/(1.d0+c)
	h=t*a(ip,iq)
	z(ip)=z(ip)-h
	z(iq)=z(iq)+h
	d(ip)=d(ip)-h
	d(iq)=d(iq)+h
	a(ip,iq)=0.d0

	do 16 j=1,ip-1
	g=a(j,ip)
	h=a(j,iq)
	a(j,ip)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 16	continue

	do 17 j=ip+1,iq-1
	g=a(ip,j)
	h=a(j,iq)
	a(ip,j)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 17	continue

	do 18 j=iq+1,n
	g=a(ip,j)
	h=a(iq,j)
	a(ip,j)=g-s*(h+g*tau)
	a(iq,j)=h+s*(g-h*tau)
 18	continue

	do 19 j=1,n
	g=v(j,ip)
	h=v(j,iq)
	v(j,ip)=g-s*(h+g*tau)
	v(j,iq)=h+s*(g-h*tau)
 19	continue
	
	nrot=nrot+1
	endif

 21	continue
 22	continue

	do 23 ip=1,n
	b(ip)=b(ip)+z(ip)
	d(ip)=b(ip)
	z(ip)=0.d0

 23	continue
 24	continue

	pause '50 iterations should never happen'

	return
	end

c****************************************************************
      function cg(j1,m1,j2,m2,j3,m3)
c
c  clebsh-gordan coefficient <j1 m1 j2 m2 | j3 m3>
c****************************************************************
      implicit real*8(a-h,o-z)
	external fact

	if(m1+m2.ne.m3) then
	cg=0.d0
	return
	endif

	if(j3.lt.abs(j1-j2)) then
	cg=0.d0
	return
	endif

	if(j3.gt.j1+j2) then
	cg=0.d0
	return
	endif

      ka=j1+j2-j3
      kb=j3+j1-j2
      kc=j2+j3-j1
      kd=j1+j2+j3+1
      del=sqrt(fact(ka)*fact(kb)*fact(kc)/fact(kd))

      cg=0.d0

      do 10 n=0,max(j1+j2-j3,j1-m1,j2+m2)

      ka1=j1+j2-j3-n
      if(ka1.lt.0.d0) go to 10
      ka2=j3-j2+m1+n
      if(ka2.lt.0.d0) go to 10
      ka3=j3-j1-m2+n
      if(ka3.lt.0.d0) go to 10
      ka4=j1-m1-n
      if(ka4.lt.0.d0) go to 10
      ka5=j2+m2-n
      if(ka5.lt.0.d0) go to 10
      an=n*1.d0

      cg=cg+(-1.d0)**n
     &  /(fact(n)*fact(ka1)*fact(ka2)*fact(ka3)*fact(ka4)*fact(ka5))

 10   continue

 20   cg=cg*sqrt(fact(j1+m1)*fact(j1-m1))
      cg=cg*sqrt(fact(j2+m2)*fact(j2-m2))
      cg=cg*sqrt(fact(j3+m3)*fact(j3-m3))

      cg=cg*sqrt(2.d0*j3+1.d0)*del

      return
      end

c***************************************************************
      subroutine matinv(nmax,c,d)
c     
c subroutine for calculating the inversion of a matrix c, 
c whose dimension is nmax
c
c***************************************************************
      implicit complex*16 (a-h,o-z)
      parameter (nlevelmax=30)
      dimension u(nlevelmax),v(nlevelmax),
     &          c(nlevelmax,nlevelmax),d(nlevelmax,nlevelmax)
      deter=1.0
      do 1 m=1,nmax
      do 1 n=1,nmax
      d(n,m)=0.0
      if(n.eq.m) d(n,m)=1.0
 1    continue
      do 10 n=1,nmax
      t=c(n,n)
      if (abs(t).lt.1.e-10) go to3
      go to 6
 3    j=n
 4    j=j+1
      if (j.gt.nmax) go to 11
      t=c(n,j)
      deter=-deter
      if (abs(t).le.1.e-10) go to 4
      do 5 k=1,nmax
      u(k)=c(n,k)
      v(k)=d(n,k)
      c(n,k)=c(j,k)
      d(n,k)=d(j,k)
      c(j,k)=u(k)
 5    d(j,k)=v(k)
 6    do 7 k=1,nmax
      if (k.eq.n) go to 7
      a=c(k,n)/c(n,n)
      do 8 l=1,nmax
      c(k,l)=c(k,l)-a*c(n,l)                              
      d(k,l)=d(k,l)-a*d(n,l)
 8    continue
 7    continue
      b=c(n,n)
      deter=b*deter
      do 10 m=1,nmax
      c(n,m)=c(n,m)/b
      d(n,m)=d(n,m)/b
 10   continue
      return
 11   print 12
 12   format('matrix not invertible')
      return
      end

                         function fact(n)
                         implicit real*8(a-h,o-z)
			 if(n.lt.0) then
			 fact=0.d0
			 return
			 endif
                         if (n.eq.0) then
                            fact=1.d0
                            go to 10
                            endif
                            fact=1.d0
                            do 20 i=1,n
                               fact=fact*i*1.d0
 20                            continue
 10                            return
                               end

C*****************************************************************
C Subroutine for the Coulomb wave function
C*****************************************************************
      SUBROUTINE JFLGAM(XD,YD,PAR,PAI,NBCHIF)
      DOUBLE PRECISION XD,YD,PAR,PAI,TEST,C,HLO2PI,PI,PIS2,PIS4
      DOUBLE PRECISION X,Y,U,V,TRA,TRA1,ALO2PI,RAC2,COSI,SINI
      DOUBLE PRECISION COS2I,SIN2I,ZMOD,DEPI
      DOUBLE PRECISION XX
      DOUBLE PRECISION ALOPI
      DOUBLE PRECISION SUPINT
      DIMENSION TEST(7),C(6)
      DATA TEST/2.9152D+7,2.2958D+3,1.4124D+2,3.9522D+1,19.6611D0,
     &12.791D0,-10.D0/
      DATA C/8.333333333333333D-2,-2.777777777777777D-3,
     &7.936507936507937D-4,-5.952380952380952D-4,8.417508417508418D-4,
     &-1.917526917526918D-3/
      DATA HLO2PI/0.918938533204672/,PI/3.141592653589793/
      DATA PIS2/1.570796326794897/,PIS4/0.785398163397448/
      DATA ALO2PI/1.837877066409345/,RAC2/0.3465735902799726/
      DATA DEPI/6.283185307179586/,ALOPI/1.1447298858494002/
      DATA SUPINT/2147483647.D0/
      NBCHIF=15
      X=DABS(XD)
      XX=X
      IF(YD)1,2,1
    1 Y=DABS(YD)
      KR=1
      I=DMOD(10.99D0-X,SUPINT)
C     TRANSLATION
      IF(I)3,3,4
    4 TRA=I
      X=X+TRA
C     LOGARITHME(X+IY) (X,Y POSITIFS)
    3 IF(X-Y)5,6,7
    5 TRA1=X/Y
      IF(TRA1)8,8,9
    8 U=DLOG(Y)
      V=PIS2
      GO TO 10
    6 U=RAC2+DLOG(X)
      V=PIS4
      GO TO 10
    9 TRA=Y*DSQRT(1.D0+TRA1*TRA1)
      TRA1=Y/X
   11 U=DLOG(TRA)
      V=DATAN(TRA1)
   10 GO TO (12,19,23),KR
    7 TRA1=Y/X
      TRA=X*DSQRT(1.D0+TRA1*TRA1)
      GO TO 11
   12 KR=2
C     DEVELOPPEMENT ASYMPTOTIQUE ( X SUPERIEUR A 10 )
      TRA=X-0.5D0
      PAI=V*TRA+Y*(U-1.D0)
      PAR=-X+HLO2PI+U*TRA-V*Y
      ZMOD=X+Y
      IF(ZMOD-TEST(1))13,13,14
   13 TRA=X*X+Y*Y
      COSI=X/TRA
      SINI=Y/TRA
      SIN2I=(SINI*COSI)+(SINI*COSI)
      COS2I=(COSI+SINI)*(COSI-SINI)
      K=1
      GO TO 15
   16 TRA=COSI*COS2I-SINI*SIN2I
      SINI=SINI*COS2I+COSI*SIN2I
      COSI=TRA
   15 PAR=PAR+C(K)*COSI
      PAI=PAI-C(K)*SINI
      K=K+1
      IF(ZMOD-TEST(K))16,16,14
C     TRANSLATION INVERSE
   17 I=I-1
      X=I
      X=XX+X
      GO TO 3
   19 PAR=PAR-U
      PAI=PAI-V
   14 IF(I-1)18,60,17
   60 IF(XD)17,61,17
C     CONTROLE DU QUADRANT
   18 IF(XD)20,61,21
   61 TRA=PI*Y
      IF(TRA-1.D-2)300,300,301
  300 PAR= TRA*(2.D0+TRA*(-2.D0+TRA*(1.333333333333333D0+TRA*(
     &-0.6666666666666666D0+TRA*(0.2666666666666666D0+TRA*(
     &-0.08888888888888888D0+TRA*0.02539682539682540D0))))))
      TRA1=-DLOG(Y)-DLOG(PAR)
      GO TO 302
  301 PAR=1.D0-DEXP(-TRA-TRA)
      TRA1=-DLOG(Y*PAR)
  302 PAR=0.5D0*(ALO2PI-TRA+TRA1)
      PAI=PAI-PIS2
   21 IF(YD)28,100,100
C     X+IY CHANGE EN -X-IY
   20 TRA=PI*Y
      PAR=ALO2PI-U-PAR-TRA
      PAI=PI-V-PAI
      TRA=DEXP(-TRA-TRA)
      X=PI*DMOD(X,2.D0)
      SINI=(1.D0-TRA)*DCOS(X)
      COSI=(1.D0+TRA)*DSIN(X)
      KR=3
      X=DABS(COSI)
      Y=DABS(SINI)
      GO TO 3
   23 IF(COSI)24,25,25
   24 V=PI-DSIGN(V,SINI)
      GO TO 26
   25 IF(SINI)27,26,26
   27 V=-V
   26 PAR=PAR-U
      PAI=PAI-V
      IF(YD)100,100,28
   28 PAI=-PAI
C     ARGUMENT DANS -PI,PI
  100 TRA=DABS(PAI/DEPI)
      IF(TRA-1.D+15)203,204,204
  204 NBCHIF=0
      PAI=0.D0
      GO TO 29
  203 IF(TRA-1.D0)205,205,206
  206 NBCHIF=DLOG10(TRA)
      NBCHIF=14-NBCHIF
      TRA=DMOD(TRA,SUPINT)
      PAI=DMOD(TRA,1.D0)*DSIGN(DEPI,PAI)
  205 IF(DABS(PAI)-PI)29,29,207
  207 PAI=PAI-DSIGN(DEPI,PAI)
      GO TO 29
C     JFLGAM REEL
    2 PAI=0.D0
      IF(XD)31,32,33
C     CONDITIONS D EXISTENCE
   32 WRITE (6,1000)
 1000 FORMAT (21H JFLGAM(0) EST INFINI)
      GO TO 50
   31 IF(X-4503599627370496.D0)34,35,35
   35 WRITE (6,1001)
 1001 FORMAT (30H ARGUMENT DE JFLGAM TROP GRAND)
      GO TO 50
   34 Y=DMOD(X,SUPINT)
      IF(Y)400,36,400
  400 IF(Y-0.99D0)33,33,405
  405 TRA=IDINT(Y+0.1D0)
      IF(DABS(Y-TRA)-5.D-15)36,36,33
   36 WRITE (6,1002)
 1002 FORMAT (28H JFLGAM (-ENTIER) EST INFINI)
   50 PAR=1.D+74
      NBCHIF=0
      GO TO 29
C     TRANSLATION
   33 I=DMOD(10.99D0-X,SUPINT)
      IF(I)37,37,38
   38 TRA=I
      X=X+TRA
C     DEVELOPPEMENT ASYMPTOTIQUE
   37 Y=DLOG(X)
      PAR=-X+HLO2PI+Y*(X-0.5D0)
      IF(X-TEST(1))39,39,43
   39 COSI=1.D0/X
      COS2I=COSI*COSI
      K=1
      GO TO 41
   42 COSI=COSI*COS2I
   41 PAR=PAR+C(K)*COSI
      K=K+1
      IF(X-TEST(K))42,42,40
C     TRANSLATION INVERSE
   44 X=X-1.D0
   48 Y=DLOG(X)
      PAR=PAR-Y
      I=I-1
   40 IF(I-1)43,49,44
   49 X=DABS(XD)
      GO TO 48
C     X NEGATIF
   43 IF(XD)45,29,29
   45 PAR=ALOPI-PAR-Y
      Y=PI*DMOD(X,2.D0)
      Y=-DSIN(Y)
      IF(Y)46,36,47
   46 Y=-Y
      PAI=PI
   47 PAR=PAR-DLOG(Y)
      ENTRY JFLGV1
   29 RETURN
      END

      SUBROUTINE YFCLEN(ETA,RO,U,UP,V,VP,SIGMA0,IDIV,NN)
      IMPLICIT COMPLEX*16(A-D,F-H),REAL*8(E,P-Z)
C
      IF(NN.EQ.1)GO TO 20
C
      ETA2=ETA*ETA
      FA=DCMPLX(1.D0,ETA)
      M=0.25*ETA+4.D1
C
C          POLYNOMES DE TCHEBICHEV JUSQU'AU RANG M
C
      K=M+2
      X=(ETA+ETA)/RO
      XX=X+X-1.D0
      T0=1.D0
      T1=XX
      XX=XX+XX
      DO 6 J=2,K
      TJ=XX*T1-T0
      T0=T1
      T1=TJ
    6 CONTINUE
      TM=T1
      TL=T0
C
C          INITIALISATION
C
      AM=(0.D0,0.D0)
      AL=(1.D-40,1.D-40)
      BN=(0.D0,0.D0)
      BM=(1.D-40,1.D-40)
      BL=(0.D0,0.D0)
      BK=DCMPLX(4.D0*DFLOAT(M+1),0.D0)*AL+BM
      F=(0.D0,0.D0)
      G =(0.D0,0.D0)
      GP=(0.D0,0.D0)
C
C          RECURRENCE DESCENDANTE
C
      K=M
   10 R=K
      TK=XX*TL-TM
      TM=TL
      TL=TK
      HK=DCMPLX(TK,0.D0)
      C1=DCMPLX(R*(R+1.D0)-ETA2,ETA*(R+R+1.D0))
      C2=(4.D0,0.D0)*DCMPLX(R+1.D0,0.D0)
      C2=C2*DCMPLX(-R-1.D0,ETA*3.D0)
      C3=FA*DCMPLX(-R-R-4.D0,ETA)
      C4=DCMPLX((7.D0*R+5.D0)/4.D0,0.D0)
      C5=DCMPLX(R+R+2.D0,0.D0)
      C6=DCMPLX((R+3.D0)/4.D0,0.D0)
      AK=(C2*AL+C3*AM-C4*BL-C5*BM-C6*BN)/C1
      J=K/2
      J=K-J-J
      IF(J)1,2,1
    1 F=F-AK
      GO TO 3
    2 F=F+AK
    3 Z=ABS(AK)
      G=G+HK*AK
      GP=GP+HK*BK
C
C          F=A0/2-A1+A2-A3+A4-A5+...
C
C          CONGRUENCE MODULO 10**60
C
      IF(Z-1.D60)4,5,5
    5 D=(1.D60,0.D0)
      AK=AK/D
      AL=AL/D
      AM=AM/D
      BK=BK/D
      BL=BL/D
      BM=BM/D
      BN=BN/D
      F=F/D
      G=G/D
      GP=GP/D
    4 IF(K)8,8,9
    9 D=DCMPLX(4.D0*R,0.D0)
      BJ=D*AK+BL
      AM=AL
      AL=AK
      BN=BM
      BM=BL
      BL=BK
      BK=BJ
      K=K-1
      GO TO 10
C
C          NORMALISATION ET CALCUL DE Z(RO)
C
    8 D=(0.5D0,0.D0)*AK
      F=F-D
      G=G-D
      GP=GP-(0.5D0,0.D0)*BK
      D=DCMPLX(0.D0,-ETA*DLOG(2.D0)+SIGMA0)
      AXPO=EXP(D)
      F=F/AXPO
      G=G/F
      GP=GP/F
C
C          CALCUL DE F ET G
C
      D=DCMPLX(0.D0,RO-ETA*DLOG(RO))
      AXPO=EXP(D)
      D=G*AXPO
      GP=AXPO*(DCMPLX(0.D0,1.D0-ETA/RO)*G-DCMPLX(X/RO,0.D0)*GP)
      V=D
      D=(0.D0,-1.D0)*D
      U=D
      VP=GP
      GP=(0.D0,-1.D0)*GP
      UP=GP
      IDIV=0
      RETURN
C
C          SERIE ORIGINE
C
   20 PI=3.141592653589793D0
      XA=0.577215664901533D0
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      Z=138.15510557964276D0
      IDIV=0
      IF(DABS(PIETA)-Z)21,21,22
   22 INDG=IDINT(PIETA/Z)
      IDIV=60*INDG
      IF(ETA.LT.0) IDIV=0
      PIETA=PIETA-Z*DFLOAT(INDG)
   21 PIETA2=0.5D0*PIETA
      P=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      CALL JFDELG(1.D0,ETA,PAR,PAI,NB)
      Z1=ETAP*(XA+XA+DLOG(2.D0)-1.D0+PAR)
      U0=0.D0
      U1=RO
      V0=1.D0
      V1=Z1*RO
      U=U0+U1
      V=V0+V1
      UP=1.D0
      VP=Z1
      XN=2.D0
      DO 104N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1
      V=V+V2
      UP=UP+XN*U2/RO
      VP=VP+XN*V2/RO
      IF(DABS(U2/U).GT.1.D-14)GOTO18
      IF(DABS(V2/V).LE.1.D-14)GOTO19
   18 U0=U1
      U1=U2
      V0=V1
      V1=V2
      XN=XN+1.D0
  104 CONTINUE
   19 PP=V+ETAP*U*DLOG(RO)
      W=U/P
      WP=UP/P
      V=P*PP
      VP=P*(VP+ETAP*(UP*DLOG(RO)+U/RO))
      U=W
      UP=WP
      RETURN
      END

      SUBROUTINE YFASYM(ETA,RAU,FO,FPO,GO,GPO,SIGO,IEXP)
      IMPLICIT REAL*8 (A-H,O-Z)
      IEXP=0
      TRB=0.D0
      RAU2=RAU+RAU
      ETAC=ETA*ETA
      CALL JFLGAM(1.D0,ETA,TRA,SIGO,NTRUC)
   40 N=0
      PS=1.D0
      GS=0.D0
      PT=0.D0
      GT=1.D0-ETA/RAU
      SF=PS
      SG=GS
      SPF=PT
      SPG=GT
   45 DENOM= DFLOAT(N+1)*RAU2
      AN= DFLOAT(N+N+1)*ETA/DENOM
      BN= (ETAC-DFLOAT(N*(N+1)))/DENOM
      PS1=AN*PS-BN*PT
      GS1=AN*GS-BN*GT-PS1/RAU
      PT1=AN*PT+BN*PS
      GT1=AN*GT+BN*GS-PT1/RAU
   42 SF=SF+PS1
      SG=SG+GS1
      SPF=SPF+PT1
      SPG=SPG+GT1
      N=N+1
      IF(N-17)46,48,44
   48 TRA=PS*PS+PT*PT
   44 TRB=PS1*PS1+PT1*PT1
      TEST=TRA-TRB
      IF(TEST)47,47,46
   46 PS=PS1
      GS=GS1
      PT=PT1
      GT=GT1
      TRA=TRB
      GOTO 45
   47 TETAO= RAU-ETA*DLOG (RAU2)+SIGO
      TRA= DSIN(TETAO)
      TRB=DCOS(TETAO)
      GO=SF*TRB-SPF*TRA
      GPO=SG*TRB-SPG*TRA
      FO=SPF*TRB+SF*TRA
      FPO=SPG*TRB+SG*TRA
      RETURN
      END

      SUBROUTINE DFCOUL(ETA,RO,F,FP,G,GP,SIGMA,L,IEXP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(*),FP(*),G(*),GP(*),IEXP(*),SIGMA(*)
      DATA DEPI/6.283185307179586D0/
      ETAC=ETA*ETA
      L1=L+1
      CALL DFCZ0(ETA,RO,F0,FP0,G0,GP0,S,I)
      F(1)=F0
      FP(1)=FP0
      G(1)=G0
      GP(1)=GP0
      IEXP(1)=I
      SIGMA(1)=S
      IF(L)1,1,2
    1 RETURN
    2 LINF=0
      IND=0
      IF((ETA.GT.0).AND.(RO.LT.(ETA+ETA)))GO TO 21
      Z=ETA+DSQRT(ETAC+DFLOAT(L*(L+1)))
      IF(RO.GE.Z)GO TO 20
    7 ROINF=ETA+DSQRT(ETAC+DFLOAT(LINF*(LINF+1)))
      IF(RO-ROINF)3,4,4
    4 IF(LINF-L)5,6,6
    5 LINF=LINF+1
      GO TO 7
    3 IND=1
    6 LIN=LINF+1
   20 XM=1.D0
      IF(IND.EQ.0)LIN=L1
      DO 8 J=2,LIN
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      F(J)=(ZAG*F(J-1)-FP(J-1))/ZIG
      FP(J)=ZIG*F(J-1)-ZAG*F(J)
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG
      GP(J)=ZIG*G(J-1)-ZAG*G(J)
      IEXP(J)=I
      SIG=SIGMA(J-1)+DATAN(ETA/(J-1))
      IPI=SIG/DEPI
      SIG=SIG-IPI*DEPI
      IF(SIG)60,50,70
   60 IF(SIG.LT.(-DEPI/2.D0))SIG=SIG+DEPI
      GO TO 50
   70 IF(SIG.GT.(DEPI/2.D0))SIG=SIG-DEPI
   50 SIGMA(J)=SIG
    8 XM=XM+1.D0
      IF(IND.EQ.0)RETURN
      GO TO 22
   21 LIN=1
   22 FTEST=F(LIN)
      FPTEST=FP(LIN)
      LMAX=LINF+25+IDINT(5.D0*DABS(ETA))
      IF(LMAX-L)9,10,10
    9 LMAX=L
   10 FI=1.D0
      FPI=1.D0
   13 XM=LMAX+1
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      FL=(ZAG*FI+FPI)/ZIG
      FPL=ZAG*FL-ZIG*FI
      IF(DABS(FL)-1.D15)26,27,27
   26 IF(DABS(FPL)-1.D15)28,27,27
   27 FL=FL*1.D-15
      FPL=FPL*1.D-15
   28 FI=FL
      FPI=FPL
      IF(LMAX-L)11,11,12
   12 LMAX=LMAX-1
      GO TO 13
   11 F(LMAX+1)=FL
      FP(LMAX+1)=FPL
      IF(LMAX-LINF)15,15,14
   14 GO TO 12
   15 FACT=FTEST/F(LIN)
      FACTP=FPTEST/FP(LIN)
      INDICE=I/60
      XM=LINF
      DO 16 J=LIN,L1
      F(J)=F(J)*FACT
      FP(J)=FP(J)*FACTP
   25 IF(J.EQ.1)GO TO 16
      ZIG=(DSQRT(ETAC+XM*XM))/XM
      ZAG=ETA/XM+XM/RO
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG
      GP(J)=ZIG*G(J-1)-ZAG*G(J)
      IF(DABS(G(J))-1.D60)17,18,18
   17 IF(DABS(GP(J))-1.D60)19,18,18
   18 G(J)=G(J)/1.D60
      GP(J)=GP(J)/1.D60
      INDICE=INDICE+1
   19 IEXP(J)=INDICE*60
      A=FP(J)*G(J)
      B=-F(J)*GP(J)
      IF(A-1.D0)29,30,30
   29 I1=IDINT(DLOG10(A))
      I2=IDINT(DLOG10(B))
      GO TO 31
   30 I1=IDINT(DLOG10(A))+1
      I2=IDINT(DLOG10(B))+1
   31 F(J)=F(J)*1.D1**(-I2)
      FP(J)=FP(J)*1.D1**(-I1)
      SIG=SIGMA(J-1)+DATAN(ETA/(J-1))
      IPI=SIG/DEPI
      SIG=SIG-IPI*DEPI
      IF(SIG)61,51,71
   61 IF(SIG.LT.(-DEPI/2.D0))SIG=SIG+DEPI
      GO TO 51
   71 IF(SIG.GT.(DEPI/2.D0))SIG=SIG-DEPI
   51 SIGMA(J)=SIG
   16 XM=XM+1.D0
      RETURN
      END

      SUBROUTINE YFIREG(ETA,RO,G0,GP0)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(ETA.LE.0.D0)GOTO250
      IF(ETA.LE.3.D0)GOTO251
      IF(ETA.LE.1.D1)GOTO252
      IF(ETA.LE.18.D0)GOTO253
      IF(ETA.LE.22.D0)GOTO254
      IF(RO.LE.0.3D0+(3.D1-ETA)/8.D1)GOTO200
C   SERIE DE TAYLOR DEPART RAU0
  300 CONTINUE
      RAU0=1.666666666666667D0*DABS(ETA)+7.5D0
      CALL YFASYM(ETA,RAU0,F0,FP0,G0,GP0,SIGMA0,IEXP)
      X=RAU0-RO
      X2=X*X
      X3=X*X2
      UNR=1.D0/RAU0
      ETR0=1.D0-2.D0*ETA*UNR
      U0=G0
      U1=-X*GP0
      U2=-0.5D0*ETR0*X2*U0
      S=U0+U1+U2
      V1=U1/X
      V2=2.D0*U2/X
      T=V1+V2
      XN=3.D0
      DO10N=3,10000
C     N=N
      XN1=XN-1.D0
      XN1=XN*XN1
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1
      S=S+U3
      V3=XN*U3/X
      T=T+V3
   16 IF(DABS(U3/S).GT.1.D-11)GO TO 11
      IF(DABS(V3/T).LE.1.D-11)GO TO 20
   11 U0=U1
      U1=U2
      U2=U3
      XN=XN+1.D0
   10 CONTINUE
   20 G0=S
      GP0=-T
      RETURN
C   SERIE  ORIGINE
  200 CONTINUE
      PI=3.141592653589793D0
      GA=0.577215664901533D0
      ETA2=ETA*ETA
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      PIETA2=0.5D0*PIETA
      B=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      CALL JFDELG(1.D0,ETA,PAR,PAI,NB)
      C1=ETAP*(GA+GA+DLOG(2.D0)-1.D0+PAR)
      U0=0.D0
      U1=RO
      V0=1.D0
      V1=C1*RO
      U=U0+U1
      V=V0+V1
      UP=1.D0
      VP=C1
      XN=2.D0
      DO 104N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1
      V=V+V2
      UP=UP+XN*U2/RO
      VP=VP+XN*V2/RO
   17 IF(DABS(U2/U).GT.1.D-14)GOTO18
      IF(DABS(V2/V).LE.1.D-14)GOTO19
   18 U0=U1
      U1=U2
      V0=V1
      V1=V2
      XN=XN+1.D0
  104 CONTINUE
   19 GP=V+ETAP*U*DLOG(RO)
      G0=B*GP
      GP0=B*(VP+ETAP*(UP*DLOG(RO)+U/RO))
      RETURN
  250 IF(RO.LE.0.5D0*ETA+9.D0)GOTO200
      GOTO 300
  251 IF(RO.LE.2.25D0+7.35D0*(3.D0-ETA))GOTO200
      GOTO 300
  252 IF(RO.LE.1.2D0+1.5D-1*(1.D1-ETA))GOTO200
      GOTO 300
  253 IF(RO.LE.0.6D0+0.75D-1*(18.D0-ETA))GOTO200
      GOTO 300
  254 IF(RO.LE.0.4D0+0.5D-1*(22.D0-ETA))GOTO200
      GOTO 300
      END

      SUBROUTINE YFRICA(ETA,RAU,FO,FPO,GO,GPO,SIGMA0,IDIV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(5),QP(5)
C
C        COEFFICIENTS RICCATI
C
      DATA G61,G62,G63,G64,G65,G66,G67,G68,G69,G610,
     &G611/ 0.1159057617187498D-01,0.3863525390624998D-01,
     &0.4660034179687498D-01,0.4858398437499998D-01,
     &0.1156514485677080D 01,0.5687475585937496D 01,
     &0.1323888288225445D 02,0.1713083224826384D 02,
     &0.1269003295898436D 02,0.5055236816406248D 01,
     &0.8425394694010415D 00/
      DATA G81,G82,G83,G84,G85,G86,G87,G88,G89,G810,G811,G812,G813,G814,
     &G815/ 0.1851092066083633D-01,0.8638429641723630D-01,
     &0.1564757823944092D 00,0.1430139541625977D 00,
     &0.1924622058868408D 00,0.8500803152720129D 01,
     &0.7265429720878595D 02,0.3057942376817972D 03,
     &0.7699689544836672D 03,0.1254157054424285D 04,
     &0.1361719536066055D 04,0.9831831171035763D 03,
     &0.4547869927883148D 03,0.1222640538215636D 03,
     &0.1455524450256709D 02/
      DATA GP61,GP62,GP63,GP64,GP65,GP66/ 0.2897644042968748D-01,
     &0.2318115234375000D 00,0.8056640625000000D 00,
     &0.1601562499999998D 01,0.3046875000000000D 00,
     &0.5624999999999998D 01/
      DATA GP81,GP82,GP83,GP84,GP85,GP86,GP87,
     &GP88/ 0.6478822231292720D-01,0.6910743713378906D 00,
     &0.3322952270507811D 01,0.9483032226562498D 01,
     &0.1769653320312499D 02,0.3478710937499998D 02,
     &0.5020312499999999D 02,0.7874999999999999D 02/
      DATA Q /0.4959570165D-1,0.8888888889D-2,0.2455199181D-2,
     &0.9108958061D-3,0.2534684115D-3/
      DATA QP /0.1728260369D0,0.3174603174D-3,0.3581214850D-2,
     &0.3117824680D-3,0.9073966427D-3/
      CALL JFLGAM(1.D0,ETA,TRA,SIGMA0,IND)
      RAU2=RAU+RAU
      RAUC=RAU*RAU
      ETAC=ETA*ETA
      ETA2=ETA+ETA
      ETARO=ETA*RAU
      ETARO2=ETARO+ETARO
      PIETA=3.141592653589793*ETA
      IND=0
      JND=0
      IG=0
      IF(ETA)20,20,21
   20 IF(-ETARO-14.0625D0)32,22,22
   22 INDICE=1
 
C             RICCATI 3
      IDIV=0
      GO TO 2
   21 IF(DABS(RAU-ETA2).LE.1.D-9)GO TO 18
      IF(RAU-ETA2)30,18,31
   31 IF(RAU-ETA2-2.D1*(ETA**0.25D0))34,33,33
   33 INDICE=0
C             RICCATI  2
      IDIV=0
      GO TO 2
   32 NN=1
      GO TO 35
   34 NN=0
   35 CALL YFCLEN(ETA,RAU,FO,FPO,GO,GPO,SIGMA0,IDIV,NN)
      RETURN
   30 IF(ETARO-12.D0)32,32,23
   23 TRA=ETA2-6.75D0*(ETA**0.4D0)
      IF(RAU-TRA)6,6,24
   24 IND=1
      JND=1
      RO=RAU
      RAU=TRA
      RAU0=TRA
C             RICCATI  1
 
    6 X=RAU/ETA2
      U= (1.D0-X)/X
      X2=X*X
      RU= DSQRT(U)
      RX= DSQRT(X)
      TRE= 1.D0/(U*RU*ETA2)
      TRB=TRE*TRE
      FI= (DSQRT((1.D0-X)*X)+DASIN(RX)-1.570796326794897D0)*ETA2
      TR1= -0.25D0*DLOG(U)
  602 TR2= -((9.D0*U+6.D0)*U+5.D0)/48.D0
  603 TR3= ((((-3.D0*U-4.D0)*U+6.D0)*U+12.D0)*U+5.D0)/64.D0
  604 TR4=- ((((((U+2.D0)*945.D0*U+1395.D0)*U+12300.D0)*U+25191.D0)
     &*U+19890.D0)*U+5525.D0)/46080.D0
  605 TR5= ((((((((-27.D0*U-72.D0)*U-68.D0)*U+360.D0)*U+2190.D0)
     &*U+4808.D0)*U+5148.D0)*U+2712.D0)*U+565.D0)/2048.D0
  606 TR6=- (((((((((G61*U+G62)*U+G63)*U+G64)*U+G65)*U+G66)*U+G67)
     &*U+G68)*U+G69)*U+G610)*U+G611
  607 TR7= ((((((((((((-81.*U-324.)*U-486.)*U-404.)*U+4509.)*U+52344.)
     &*U+233436.)*U+567864.)*U+838521.)*U+775884.)*U+441450.)
     &*U+141660.)*U+19675.) /6144.
  608 TR8= (((((((((((((G81*U+G82)*U+G83)*U+G84)*U+G85)*U+G86)*U+G87)
     &*U+G88)*U+G89)*U+G810)*U+G811)*U+G812)*U+G813)*U+G814)*U+G815
      PSIP=PSIP+TRA
      XXX=138.1551055796428D0
      FI= FI+TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))
      PSI=-FI
      INDG=IDINT(PSI/XXX)
      IDIV=60*INDG
      TRA= TR1+TRB*(TR3+TRB*(TR5+TRB*TR7))
      FI=FI+TRA
      PSI=PSI+TRA
 
      FIP=RU*ETA2
      TRA=1.D0/(X2*U)
      TR1=0.25D0
      TRE= TRE/(X2*X2*U)
      TRB=TRB/(X2*X2)
  702 TR2= -(8.D0*X-3.D0)/32.D0
  703 TR3= ((24.D0*X-12.D0)*X+3.D0)/64.D0
  704 TR4= (((-1536.D0*X+704.D0)*X-336.D0)*X+63.D0)/2048.D0
  705 TR5= ((((1920.D0*X-576.D0)*X+504.D0)*X-180.D0)*X+27.D0)/1024.D0
  706 TR6 = ((((-GP66*X+GP65)*X-GP64)*X+GP63)*X-GP62)*X+GP61
  707 TR7= - ((((((-40320.D0*X-10560.D0)*X-13248.D0)*X+7560.D0)
     &*X-3132.D0)*X+756.D0)*X-81.D0) / 2048.D0
  708 TR8 =- ((((((GP88*X+GP87)*X+GP86)*X-GP85)*X+GP84)*X-GP83)*X+GP82)
     &*X-GP81
      FIP=FIP +TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))
      TRA= TRA*(TR1+TRB*(TR3+TRB*(TR5+TRB*TR7)))
      FIP=FIP+TRA
      PSIP=-FIP
      IF(INDG.EQ.0)GO TO 8
      PSI=PSI-XXX*DFLOAT(INDG)
      FI =FI +XXX*DFLOAT(INDG)
    8 FO=0.5D0*DEXP(FI)
      GO= DEXP(PSI)
      FPO= FO*FIP/ETA2
      GPO=GO*PSIP/ETA2
      IF(JND.EQ.0)RETURN
      RAU=RO
      GO=FO
      GPO=FPO
   27 X=RAU0-RO
      X2=X*X
      X3=X*X2
      UNR=1.D0/RAU0
      ETR0=1.D0-2.D0*ETA*UNR
      U0=GO
      U1=-X*GPO
      U2=-0.5D0*ETR0*X2*U0
      S=U0+U1+U2
      V1=U1/X
      V2=2.D0*U2/X
      T=V1+V2
      XN=3.D0
 
      DO10N=3,10000
C     N=N
      XN1=XN-1.D0
      XN1=XN*XN1
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1
      S=S+U3
      V3=XN*U3/X
      T=T+V3
   16 IF(DABS(U3/S).GT.1.D-10)GO TO 11
      IF(DABS(V3/T).LE.1.D-10)GO TO 17
   11 U0=U1
      U1=U2
      U2=U3
      XN=XN+1.D0
   10 CONTINUE
   17 IF(IG)25,26,25
   25 GO=S
      GPO=-T
      FO=HO
      FPO=HPO
      RETURN
   26 HO=S
      HPO=-T
   18 ET0=ETA**(0.166666666666667)
      ETAD=ETAC*ETAC
      ET=ETA**(0.6666666666666667)
      ET1=ET*ET
      ET2=ET1*ET1
      ET3=ET2*ET
      ET4=ETAD*ET
      ET5=ET4*ET
      FO=1.D0-Q(1)/ET1-Q(2)/ETAC-Q(3)/ET3-Q(4)/ETAD-Q(5)/ET5
      GO=1.D0+Q(1)/ET1-Q(2)/ETAC+Q(3)/ET3-Q(4)/ETAD+Q(5)/ET5
      FPO=1.D0+QP(1)/ET+QP(2)/ETAC+QP(3)/ET2+QP(4)/ETAD+QP(5)/ET4
      GPO=1.D0-QP(1)/ET+QP(2)/ETAC-QP(3)/ET2+QP(4)/ETAD-QP(5)/ET4
      FO=0.7063326373D0*ET0*FO
      GO=1.223404016D0*ET0*GO
      FPO=0.4086957323D0*FPO/ET0
      GPO=-0.7078817734D0*GPO/ET0
      IDIV=0
      IF(IND.EQ.0)RETURN
      IG=1
      RAU0=ETA2
      GO TO 27
    2 X=ETA2/RAU
      X2=X*X
      U=1.D0-X
      RU= DSQRT(U)
      U3=U*U*U
      TRD= 1.D0/(U3*ETA2*ETA2)
      TRC=X2*TRD
      TRE=1.D0/(U*RU*ETA2)
      FI= -0.25D0*DLOG(U)
      TRB=TRD/64.D0
      TR3= (((3.D0*U-4.D0)*U-6.D0)*U+12.D0)*U-5.D0
  501 TR5= ((((((((-27.D0*U+72.D0)*U-68.D0)*U-360.D0)*U+2190.D0)
     &*U-4808.D0)*U+5148.D0)*U-2712.D0)*U+565.D0)/32.D0
  502 TR7= ((((((((((((81.D0*U-324.D0)*U+486.D0)*U-404.D0)*U-4509.D0)
     &*U+52344.D0)*U-233436.D0)*U+567864.D0)*U-838521.D0)*U+775884.D0)
     &*U-441450.D0)*U+141660.D0)*U-19675.D0)/96.D0
      FI= FI+TRB*(TR3+TRD*(TR5+TRD*TR7))
 
      FIP=0.25D0/U
      TRB=3.D0*TRC/(64.D0*U)
      TR3= (X-4.D0)*X+8.D0
  511 TR5= ((((9.D0*X-60.D0)*X+168.D0)*X-192.D0)*X+640.D0)/16.D0
  512 TR7= ((((((-27.D0*X+252.D0)*X-1044.D0)*X+2520.D0)*X-4416.D0)
     &*X-3520.D0)*X-13440.D0)/32.D0
      FIP =FIP+TRB*(TR3+TRC*(TR5+TRC*TR7))
      TRA= DABS((RU-1.D0)/(RU+1.D0))
      PSI= (0.5D0*DLOG(TRA)+RU/X)*ETA2+0.785398163397448D0
      TR2= -((9.D0*U-6.D0)*U+5.D0)/48.D0
  521 TR4= ((((((U-2.D0)*945.D0*U+1395.D0)*U-12300.D0)*U+25191.D0)
     &*U-19890.D0)*U+5525.D0)/46080.D0
  522 TR6 = (((((((((-G61*U+G62)*U-G63)*U+G64)*U-G65)*U+G66)*U-G67)
     &*U+G68)*U-G69)*U+G610)*U-G611
  523 TR8= (((((((((((((G81*U-G82)*U+G83)*U-G84)*U+G85)*U-G86)*U+G87)
     &*U-G88)*U+G89)*U-G810)*U+G811)*U-G812)*U+G813)*U-G814)*U+G815
      PSI= PSI+TRE*(TR2+TRD*(TR4+TRD*(TR6+TRD*TR8)))
      PSIP = -RU*ETA2/X2
      TRB=TRE*X/U
      TR2= (3.D0*X-8.D0)/32.D0
  531 TR4= - (((63.D0*X-336.D0)*X+704.D0)*X-1536.D0)/2048.D0
  532 TR6 = ((((GP61*X-GP62)*X+GP63)*X-GP64)*X+GP65)*X-GP66
  533 TR8 = ((((((-GP81*X+GP82)*X-GP83)*X+GP84)*X-GP85)*X+GP86)*X+GP87)
     &*X+GP88
      PSIP =PSIP+ TRB*(TR2+TRC*(TR4+TRC*(TR6+TRC*TR8)))
      TRA= DEXP(FI)
      FO= TRA*DSIN(PSI)
      GO= TRA*DCOS(PSI)
      IF(INDICE)535,536,535
  535 TRA=FO
      FO=GO
      GO=-TRA
  536 TRA=-ETA2/RAUC
      FPO=(FIP*FO+PSIP*GO)*TRA
      GPO=(FIP*GO-PSIP*FO)*TRA
      RETURN
      END

      SUBROUTINE DFCZ0(ETA,RO,F0,FP0,G0,GP0,SIGMA0,IEXP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A1(110),A2(110),B1(110),B2(110)
      IF(RO)2,2,1
    2 WRITE (6,1000)
 1000 FORMAT (21H RO NEGATIF OU NUL **)
      RETURN
    1 IF(ETA-30.D0)3,3,4
    4 IF(DABS(ETA)-5.D2)28,28,29
   28 CALL YFRICA(ETA,RO,F0,FP0,G0,GP0,SIGMA0,IEXP)
      RETURN
   29 WRITE (6,1001)
 1001 FORMAT (42H VALEUR ABSOLUE DE ETA SUPE-&EU-E A 500 **)
      RETURN
    3 IF(ETA+8.D0)4,5,5
    5 IF(ETA)6,7,6
    7 F0=DSIN(RO)
      G0=DCOS(RO)
      FP0=G0
      GP0=-F0
      IEXP=0
      SIGMA0=0.D0
      RETURN
    6 BORNE=1.666666666666667D0*DABS(ETA)+7.5D0
      IF(RO-BORNE)8,9,9
    9 CALL YFASYM(ETA,RO,F0,FP0,G0,GP0,SIGMA0,IEXP)
      RETURN
    8 IF(ETA-10.D0)10,11,11
   10 IF(ETA)12,12,13
   13 IF(RO-2.D0)14,12,12
   11 IF(ETA-(5.D0*RO+6.D1)/7.D0)12,12,14
   12 CALL YFASYM(ETA,BORNE,F0,FP0,G0,GP0,SIGMA0,IEXP)
      H=BORNE
      DH=F0/H
      IF(ETA)20,21,21
   20 N=-0.5D0*ETA+5.D0
      GO TO 22
   21 N=ETA/5.D0+5.D0
   22 N=5*(N+1)
      Z=4.D0/H
      Y=1.D0-(ETA+ETA)*Z
      A1(N+2)=1.D-55
      A1(N+3)=0.D0
      A1(N+4)=1.D-64
      B1(N+3)=1.D-50
      B1(N+4)=1.D-68
      A2(N+2)=0.D0
      A2(N+3)=1.D-74
      A2(N+4)=1.D-53
      B2(N+3)=0.D0
      B2(N+4)=1.D-66
      M=N+2
      DI=N
      DO 23 II=2,M
      I=M-II+2
      B1(I)=B1(I+2)+Z*(DI+1.D0)*A1(I+1)
      S=A1(I+2)+Y*(A1(I+1)-A1(I))
      Q=(DI+2.D0)*B1(I)+(DI-1.D0)*B1(I+1)
      A1(I-1)=S-Z*Q
      B2(I)=B2(I+2)+Z*(DI+1.D0)*A2(I+1)
      S=A2(I+2)+Y*(A2(I+1)-A2(I))
      Q=(DI+2.D0)*B2(I)+(DI-1.D0)*B2(I+1)
      A2(I-1)=S-Z*Q
      IF(I.GE.N)GO TO 23
      D=-(B2(I+2)+B2(I))/(B1(I+2)+B1(I))
      DO 24 J=I,M
      A2(J)=A2(J)+D*A1(J)
      B2(J)=B2(J)+D*B1(J)
   24 CONTINUE
      A2(I-1)=A2(I-1)+D*A1(I-1)
   23 DI=DI-1.D0
      Q=A1(3)-A1(1)
      C=A2(3)-A2(1)
      C=Q/C
      X1=0.5D0*(A1(2)-C*A2(2))
      DO 25 I=3,M
      X1=X1+A1(I)-C*A2(I)
   25 CONTINUE
      X1=DH/X1
      X2=-C*X1
      DO 26 I=2,M
      B1(I)=X1*B1(I)+X2*B2(I)
      A1(I)=X1*A1(I)+X2*A2(I)
   26 CONTINUE
      A1(1)=X1*A1(1)+X2*A2(1)
      B1(1)=0.D0
      X=RO/H
      Y=2.D0*X-1.D0
      T1=1.D0
      T2=Y
      RESULT=0.5D0*A1(2)+Y*A1(3)
      DERIVE=0.5D0*B1(2)+Y*B1(3)
      DO 27 I=2,N
      TI=2.D0*Y*T2-T1
      T1=T2
      T2=TI
      RESULT=RESULT+TI*A1(I+2)
      DERIVE=DERIVE+TI*B1(I+2)
   27 CONTINUE
      F0=RESULT*RO
      FP0=DERIVE*RO+RESULT
      GO TO 30
C   SERIE ORIGINE REGULIERE
   14 PI=3.141592653589793D0
      CALL JFLGAM(1.D0,ETA,TRA,SIGMA0,NTRUC)
      IEXP=0
      RO2=RO*RO
      ETAP=ETA+ETA
      PIETA=PI*ETA
      PIETA2=0.5D0*PIETA
      B=DEXP(PIETA2)*DSQRT(DSINH(PIETA)/PIETA)
      U0=0.D0
      U1=RO
      U=U0+U1
      UP=1.D0
      XN=2.D0
      DO 15 N=2,10000
      XN1=XN*(XN-1.D0)
      U2=(ETAP*RO*U1-RO2*U0)/XN1
      U=U+U2
      UP=UP+XN*U2/RO
   17 IF(DABS(U2/U).LT.1.D-10)GO TO 19
   18 U0=U1
      U1=U2
      XN=XN+1.D0
   15 CONTINUE
   19 F0=U/B
      FP0=UP/B
   30 CALL YFIREG(ETA,RO,G0,GP0)
      RETURN
      END

      SUBROUTINE JFDELG (XD,YD,PAR,PAI,NBCHIF)
      DOUBLE PRECISION XD,YD,PAR,PAI,TEST,C,PI
      DOUBLE PRECISION X,Y,U,V,TRA,TRA1,COSI,SINI
      DOUBLE PRECISION COS2I,SIN2I,ZMOD,DEPI
      DOUBLE PRECISION TRB,XX
      DOUBLE PRECISION RAC2,PIS4
      DOUBLE PRECISION SUPINT
      DIMENSION TEST(7),C(6)
      DATA TEST/2.9152D+7,2.2958D+3,1.4124D+2,3.9522D+1,19.6611D0,
     &12.791D0,-10.D0/
      DATA RAC2/0.3465735902799726/,PIS4/0.785398163397448/
      DATA C/8.333333333333333D-2,-8.33333333333333D-3,
     &3.968253968253968D-3,-4.166666666666667D-3,7.575757575757576D-3,
     &-2.109279609279609D-2/
      DATA PI/3.141592653589793/
      DATA DEPI/6.283185307179586/
      DATA SUPINT/2147483647.D0/
      X=DABS(XD)
      XX=X
      NBCHIF=15
      IF(YD)1,2,1
    1 Y=DABS(YD)
      KR=1
      I=DMOD(10.99D0-X,SUPINT)
C     TRANSLATION
      IF(I)3,3,4
    4 TRA=I
      X=X+TRA
C     LOGARITHME(X+IY) (X,Y POSITIFS)
    3 IF(X-Y)5,6,7
    5 TRA1=X/Y
      TRB=1.D0+TRA1*TRA1
      TRA=Y*DSQRT(TRB)
      SINI=1./(TRB*Y)
      COSI=SINI*TRA1
      TRA1=Y/X
      GO TO 11
    6 U=RAC2+DLOG(X)
      V=PIS4
      SINI=0.5D0/X
      COSI=SINI
      GO TO 10
    7 TRA1=Y/X
      TRB=1.D0+TRA1*TRA1
      TRA=X*DSQRT(TRB)
      COSI=1./(TRB*X)
      SINI=COSI*TRA1
   11 U=DLOG(TRA)
      V=DATAN(TRA1)
C     DEVELOPPEMENT ASYMPTOTIQUE ( X SUPERIEUR A 10 )
   10 PAR=U-0.5*COSI
      PAI=V+0.5*SINI
      ZMOD=X+Y
      IF(ZMOD-TEST(1))13,13,14
   13 SIN2I=(SINI*COSI)+(SINI*COSI)
      COS2I=(COSI+SINI)*(COSI-SINI)
      SINI=SIN2I
      COSI=COS2I
      K=1
      GO TO 15
   16 TRA=COSI*COS2I-SINI*SIN2I
      SINI=SINI*COS2I+COSI*SIN2I
      COSI=TRA
   15 PAR=PAR-C(K)*COSI
      PAI=PAI+C(K)*SINI
      K=K+1
      IF(ZMOD-TEST(K))16,16,14
C     TRANSLATION INVERSE
   17 I=I-1
      X=I
      X=XX+X
   56 IF(X-Y)55,55,57
   55 TRA1=X/Y
      TRB=X*TRA1+Y
      SINI=1.D0/TRB
      COSI=TRA1/TRB
      GO TO 19
   57 TRA1=Y/X
      TRB=X+Y*TRA1
      COSI=1.D0/TRB
      SINI=TRA1/TRB
   19 PAR=PAR-COSI
      PAI=PAI+SINI
   14 IF(I)18,18,17
 
C     CONTROLE DU QUADRANT
   18 IF(XD)20,61,21
   61 TRA=PI*Y
      IF(TRA-1.D-2)300,300,301
  300 TRB= TRA*(2.D0+TRA*(-2.D0+TRA*(1.333333333333333D0+TRA*(
     &-0.6666666666666666D0+TRA*(0.2666666666666666D0+TRA*(
     &-0.08888888888888888D0+TRA*0.02539682539682540D0))))))
      TRB=(2.D0-TRB)/TRB
      GO TO 302
  301 TRB= DEXP(-TRA-TRA)
      TRB=(1.D0+TRB)/(1.D0-TRB)
  302 PAI=0.5D0*(1.D0/Y+PI*TRB)
   21 IF(YD)28,100,100
C     X+IY CHANGE EN -X-IY
   20 TRA=DEXP(-DEPI*Y)
      TRB=TRA*TRA
      COS2I=DEPI*DMOD(X,1.D0)
      SIN2I=-2.D0*TRA*DCOS(COS2I)+1.D0+TRB
      PAR=PAR+COSI+DEPI*TRA*DSIN(COS2I)/SIN2I
      PAI=PAI-SINI+PI*(TRB-1.D0)/SIN2I
      IF(YD)100,100,28
   28 PAI=-PAI
   35 WRITE (6,1001)
 1001 FORMAT (30H ARGUMENT DE JFDELG TROP GRAND)
      GO TO 50
   34 Y=DMOD(X,SUPINT)
      IF(Y) 400,36,400
  400 IF(Y-0.99D0) 33,33,405
  405 TRA= IDINT(Y+0.1D0)
      IF(DABS(Y-TRA)-5.D-15)36,36,33
   31 IF(X-4503599627370496.D0)34,35,35
 
C     ARGUMENT DANS -PI,PI
  100 TRA=DABS(PAI/DEPI)
      IF(TRA-1.D+15)203,204,204
  204 NBCHIF=0
      PAI=0.D0
      GO TO 29
  203 IF(TRA-1.D0)205,205,206
  206 NBCHIF=DLOG10(TRA)
      NBCHIF=14-NBCHIF
      TRA=DMOD(TRA,SUPINT)
      PAI=DMOD(TRA,1.D0)*DSIGN(DEPI,PAI)
  205 IF(DABS(PAI)-PI)29,29,207
  207 PAI=PAI-DSIGN(DEPI,PAI)
      GO TO 29
C        DELGAMMA REEL
    2 PAI=0.D0
      IF(XD)31,32,33
C     CONDITIONS D EXISTENCE
   32 WRITE (6,1000)
 1000 FORMAT (21H JFDELG(0) EST INFINI)
      GO TO 50
   36 WRITE (6,1002)
 1002 FORMAT (28H JFDELG (-ENTIER) EST INFINI)
   50 PAR=1.D+74
      NBCHIF=0
      GO TO 29
C     TRANSLATION
   33 I=DMOD(10.99D0-X,SUPINT)
      IF(I)37,37,38
   38 TRA=I
      X=X+TRA
C     DEVELOPPEMENT ASYMPTOTIQUE
   37 Y=DLOG(X)
      PAR=Y-0.5D0/X
      IF(X-TEST(1))39,39,43
   39 COS2I=1.D0/(X*X)
      COSI=COS2I
      K=1
      GO TO 41
   42 COSI=COSI*COS2I
   41 PAR=PAR-C(K)*COSI
      K=K+1
      IF(X-TEST(K))42,42,40
C     TRANSLATION INVERSE
   44 I=I-1
      X=I
      X=XX+X
      PAR=PAR-1.D0/X
   40 IF(I)43,43,44
C     X NEGATIF
   43 IF(XD)45,29,29
   45 PAR=PAR+1.D0/X
      Y=PI*DMOD(X,2.D0)
      PAR=PAR+PI*DCOS(Y)/DSIN(Y)
      ENTRY JFDEV1
   29 RETURN
      END




