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
! SV2: effective plastic strain.
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
      integer num_iter, max_num_iter, i, j
      real*8 E, nu, A, B, n, m, Tm, Tr, C, epsilon_dot_zero, D1, D2,
     1 D3, D4, D5, beta, Cp, unit_conversion_factor, 
     2 convergence_tolerance_factor, tolerance, mu, lambda, 
     3 eps, Temp, sigmaY, equiv_stress, pl_strain_inc, 
     4 eps_iter,
     5 C_mat(ndir+nshr, ndir+nshr), stress_old(ndir+nshr),
     6 trial_stress_inc(ndir+nshr), dev_stress(ndir+nshr), 
     7 pl_strain_dir(ndir+nshr), C_PSD(ndir+nshr), 
     8 stress_inc_iter(ndir+nshr)
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
      max_num_iter = props(20)



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

              eps = stateOld(i, 2)                                      ! Previous effective plastic strain
              Temp = stateOld(i, 3)                                     ! Temperature
              sigmaY = stateOld(i, 4)                                   ! Yield stress
              stress_old(1:ndir+nshr) = stressOld(k, 1:ndir+nshr)       ! Stress tensor in previous step
              
              trial_stress_inc = 0.0d0         

              ! Compute Hooke's Law as a function of the strain increment
              ! Calculates the stress trial increment tensor assuming a pure elastic response
              call elastic_stress(strainInc, trial_stress_inc, stressOld, 
     1                            lambda, mu, ndir, nshr) 

              ! Start the calculations to obtain the direction and magnitude of the
              ! plastic strain increment
              ! Direction = df/dSigma
              ! Magnitude = dLambda (plastic multiplier)
              ! dEpsilon^p = dLambda*df/dSigma = dp*3/2*DevStress/EquivalentStress
              
              call equivalent_stress(equiv_stress, trial_stress_inc, 
     1                               dev_stress, ndir, nshr)

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
              pl_strain_inc_max = equiv_stress/(2.d0*mu)                ! Maximum plastic strain increment


              !##################################################################################################
              ! ######## Iteration to obtain the initial plastic strain increment (pl_strain_inc) starts ########
              ! Return mapping algorithm
              
              num_iter = 0    ! number of iterations
              
              do while (num_iter < max_num_iter) 
                  ! Iteration effective plastic strain = effective plastic strain + increment
                  eps_iter = eps + pl_strain_inc  
                  ! Effective plastic strain increment rate  
                  eps_rate = pl_strain_inc/dtArray(1)
                  
                  ! Compute the stress increment corrected with the plastic corrector for this
                  ! iteration (stress_inc_iter)
                  stress_inc_iter = trial_stress_inc - pl_strain_inc*C_PSD

                  ! Calculate the equivalent stress with the corrected stress increment in this
                  ! iteration
                  call equivalent_stress(equiv_stress, stress_inc_iter, 
     1                               dev_stress, ndir, nshr)


        
        c evaluating new yield stress for the increment 'ep'
                                 call func_syield(A, B, epbarN, n, Tr, Tm, tempN, m, C, eprateN,
             1			rate0, sigyN)
                if (sigyN .lt. sigyT) sigyN = sigyT
        
                f = sigeqv - sigyN
                if (abs(f) .lt. tol) exit 			! out of the iteration
        
        c elastic criteria, ep = 0 & f < 0
                if ((ep .eq. 0.d0) .and. (f .lt. 0.d0)) go to 12
        
        c rearranging the margins and new plastic strain increm.
                if ((f .ge. 0.d0) .and. (ep .ge. epmin)) epmin = ep
                if ((f .lt. 0.d0) .and. (ep .lt. epmax)) epmax = ep
                ep = 0.5d0 * (epmax + epmin)
        
        c restoring the plastic strain increment from the last step
                if (iter .eq. 1) ep = stateOld(k, 5)
              end do 
        
         12   stressNew(k, 1:6) = sNew(1:6)

         num_iter = num_iter + 1

         if (iter .gt. Niter-1.) then
            print*, 'too many iterations, iter = ', iter
            call XPLB_EXIT
          end if
        









        c work, plastic work and temp
              dWork = dot_product( 0.5d0*(sOld(1:6) + sNew(1:6)),
             1  strainInc(k, 1:6))
              dPwork = 0.5d0 * ep * sigeqv
              enerInternNew(k) = enerInternOld(k) + dWork/rho
              enerInelasNew(k) = enerInelasOld(k) + dPwork/rho
        
        c updating state variables
                           stateNew(k, 1) = 1.d0							! to flag the initial step
                          stateNew(k, 2) = epbarN						! plastic strain
                          stateNew(k, 3) = tempT + WH*TQ*dPwork/rho/Cp	! temperature
                          stateNew(k, 4) = sigyN						! yield stress
                          stateNew(k, 5) = ep 							! plastic strain incre.
                          stateNew(k, 6) = stateOld(k, 6) + iter ! number of iterations
        

          
          
            endif

      enddo


! Subroutines
      include 'utils.for'


      return
      end