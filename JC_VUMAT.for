! ########################################################################
! User subroutine to model the Jonhson-Cook plasticity, damage and the 
! Taylor-Quinney conversion of mechanical work into heat during plastic
! deformation.
!
! Abaqus version: Abaqus 2022
! Intel oneAPI Compiler 2022.0.2 intel64
! Visual Studio 2019
!
! Author: Mauro Francisco Arcidiacono
! ########################################################################
! ########################################################################
! 
! State Variable (SV) Definitions
! SV1: initiation flag. If 0, first step, otherwise, 1.
! SV2: equivalent plastic strain.
! SV3: temperature.
! SV4: Von Mises yield stress.
! SV5: plastic strain increment.
! SV6: total number of iterations.
!
! ########################################################################


      subroutine vumat(
! Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
! Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
!
      include 'vaba_param.inc'
      parameter (i_info_AnnealFlag = 1, 
     *     i_info_Intpt    = 2, ! Integration station number
     *     i_info_layer  = 3, ! Layer number
     *     i_info_kspt   = 4, ! Section point number in current layer
     *     i_info_effModDefn = 5, ! =1 if Bulk/ShearMod need to be defined
     *     i_info_ElemNumStartLoc   = 6) ! Start loc of user element number
!
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), dtArray(2*(nblock)+1), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock), 
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock), jInfoArray(*)
!
      character*80 cmname
      integer num_iter, i, j
      real*8 E, nu, A, B, n, m, Tm, Tr, C, epsilon_dot_zero, D1, D2,
     1 D3, D4, D5, beta, Cp, unit_conversion_factor, 
     2 convergence_tolerance_factor, tolerance, mu, lambda, 
     3 eps, Temp, sigmaY, hyd_stress_trial 
     4 C_mat(ndir+nshr, ndir+nshr), stress_old(ndir+nshr),
     5 stress_trial(ndir+nshr), dev_stress(ndir+nshr), 
     6 pl_strain_dir(ndir+nshr), C_PSD(ndir+nshr)
!
      pointer (ptrjElemNum, jElemNum)
      dimension jElemNum(nblock)
!
      lAnneal = jInfoArray(i_info_AnnealFlag) 
      iLayer = jInfoArray(i_info_layer)
      kspt   = jInfoArray(i_info_kspt)
      intPt  = jInfoArray(i_info_Intpt)
      iUpdateEffMod = jInfoArray(i_info_effModDefn)
      iElemNumStartLoc = jInfoArray(i_info_ElemNumStartLoc)
      ptrjElemNum = loc(jInfoArray(iElemNumStartLoc))
!

! ############################ Input parameters ############################

      ! Elastic properties
      E = props(1)                          ! Young's Modulus
      nu = props(2)                         ! Poisson's Ratio

      ! Johnson-Cook (JC) plasticity parameters
      A = props(3)                          ! Parameter A
      B = props(4)                          ! Parameter B
      n = props(5)                          ! Parameter n 
      m = props(6)                          ! Parameter m
      Tm = props(7)                         ! Melting temperature
      Tr = props(8)                         ! Reference/Room temperature
      C = props(9)                          ! Parameter C
      epsilon_dot_zero = props(10)          ! Reference strain rate

      ! Johnson-Cook damage parameters
      D1 = props(11)
      D2 = props(12)
      D3 = props(13)
      D4 = props(14)
      D5 = props(15)

      ! Taylor-Quinney (TQ) parameters
      beta = props(16)                      ! Taylor-Quinney coefficient
      Cp = props(17)                        ! Heat capacity
      ! Factor for unit conversion to compute the temperature using the 
      ! TQ equation. Use only when needed according to your units,
      ! otherwise, assign 1. 
      unit_conversion_factor = props(18)    

      ! Convergence tolerance for the increment of effective strain calculation
      convergence_tolerance_factor = props(19)
      tolerance = E*convergence_tolerance_factor

      ! Maximum number of iterations
      num_iter = props(20)



! ############################ Elasticity Matrix ############################

      ! Lame parameters
      ! Lame first parameter
      lambda = (E*nu)/((1 + nu)*(1 - 2.0d0*nu))  
      ! Lame second parameter  
      mu = E/(2*(1 + nu))

      ! Elasticity matrix (C_mat)
      ! Linear elastic, homogeneous and isotropic material 
      ! Initialize the elasticity matrix to zero
      C_mat = 0.0d0

	do i = 1, ndir
	    do j = 1, ndir
              if (i == j) then
                  C_mat(i, j) = lambda + 2.d0*mu
              else
	            C_mat(i, j) = lambda
              end if
	    end do 
	end do 
	do i = ndir + 1, ndir + nshr
	    C_mat(i, i) = 2.d0*mu
	end do 



! ###########################################################################################################      
! ############################ Loop through each element to perform computations ############################
! ###########################################################################################################  

      ! Global loop starting point -> loops through each element of the model
      ! The index of the elements is i
      do i = 1, nblock
          
          ! ############################## Initial state ###############################   
          if (stateOld(i, 1) == 0) then
              
              ! Compute Hooke's Law as a function of the strain increment
              call elastic_stress(strainInc, stressNew, stressOld, 
     1                            lambda, mu, ndir, nshr)

              stateNew(i, 1) = 1.d0           ! Initiation flag
              stateNew(i, 2) = 0.d0           ! Equivalent plastic strain
              stateNew(i, 3) = Tr             ! Initial temperature
              stateNew(i, 4) = 0.d0           ! Yield stress
              stateNew(i, 5) = 0.d0		    ! Plastic strain increment in the last step
              stateNew(i, 6) = 1	          ! Number of iterations

          ! ############################ From step 2 onward ############################    
          else

              eps = stateOld(i, 2)                                      ! Equivalent plastic strain
              Temp = stateOld(i, 3)                                     ! Temperature
              sigmaY = stateOld(i, 4)                                   ! Yield stress
              stress_old(1:ndir+nshr) = stressOld(k, 1:ndir+nshr)       ! Stress old
              
              stress_trial = 0.0d0          

              ! Compute Hooke's Law as a function of the strain increment
              ! Calculates the stress trial tensor assuming a pure elastic response
              call elastic_stress(strainInc, stress_trial, stressOld, 
     1                            lambda, mu, ndir, nshr) 

              
              ! Start the calculations to obtain the direction and magnitude of the
              ! plastic strain increment
              ! Direction = df/dSigma
              ! Magnitude = dLambda (plastic multiplier)
              ! dEpsilon^p = dLambda*df/dSigma = dp*3/2*DevStress/EquivalentStress
              
              ! Calculate the hydrostatic component of the stress trial tensor
              hyd_stress_trial = sum(stress_trial(1:ndir))/3.d0
              ! Calculate the deviatoric tensor
              dev_stress(1:ndir) = stress_trial(1:ndir) - hyd_stress_trial
              dev_stress(ndir+1:ndir+nshr) = stress_trial(ndir+1:ndir+nshr)

              ! Compute the equivalent stress
              equiv_stress = sqrt(3.d0/2.d0 *(dev_stress(1)**2.d0 + dev_stress(2)**2.d0 +
     1        dev_stress(3)**2.d0 + 2.d0*dev_stress(4)**2.d0 + 2.d0*dev_stress(5)**2.d0 +
     2        2.d0*dev_stress(6)**2.d0 ))

              ! Plastic strain direction
              if (equiv_stress /= 0) then
                  pl_strain_dir = (3.d0*dev_stress)/(2.d0*equiv_stress)
              else
                  pl_strain_dir = 0

              ! Elasticity matrix C times the plastic strain direction to calculate
              ! later the plastic corrector
              C_PSD = matmul(C_mat, pl_strain_dir)

              
              ! Initial value of the strain increment
              pl_strain_inc = 0.d0						      ! Initial plastic strain increment
              pl_strain_inc_min = 0.d0                                  ! Minimum plastic strain increment
              pl_strain_inc_max = sigeqv/(2.d0*mu)                      ! Maximum plastic strain increment
          
          
            endif

      enddo


! Subroutines
      include 'utils.for'


      return
      end