! ########################################################################
! User subroutine to model the Johnson-Cook plasticity, damage and the 
! Taylor-Quinney conversion of mechanical work into heat during plastic
! deformation.
!
! Abaqus version: Abaqus 2022
! Intel Fortran Compiler 2021.11
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
! SV6: parameter D of the Johnson-Cook damage model.
! SV7: damage evolution parameter.
! SV8: total number of iterations.
! SV9: status of the element. If 1, the element is active, otherwise, the
! element was deleted.
! SV10 to SV16: stress tensor in the previous step.
!
! Note: this code was made to run the job using double precision due to
! variable declaration inside the subroutines (real*8).
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
     1 D3, D4, D5, beta, Cp, D, pl_disp_failure, dmg_evol,
     2 convergence_tolerance_factor, tolerance, mu, lambda, 
     3 eps, Temp, sigmaY, equiv_stress, pl_strain_inc, trace_strainInc, 
     4 eps_iter, equiv_stress_jc, f, pl_strain_inc_min, 
     5 pl_strain_inc_max, dWork, dPwork, rho, equiv_strain_fracture,
     6 C_mat(ndir+nshr, ndir+nshr), stress_old(ndir+nshr),
     7 trial_stress(ndir+nshr), dev_stress(ndir+nshr), 
     8 pl_strain_dir(ndir+nshr), C_PSD(ndir+nshr), 
     9 corrected_stress_iter(ndir+nshr)
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

      ! Plastic displacement at failure
      pl_disp_failure = props(16)

      ! Taylor-Quinney (TQ) parameters
      beta = props(17)                      ! Taylor-Quinney coefficient
      Cp = props(18)                        ! Heat capacity
      rho = props(19)                       ! Density

      ! Convergence tolerance for the increment of effective strain calculation
      convergence_tolerance_factor = props(20)
      tolerance = E*convergence_tolerance_factor

      ! Maximum number of iterations
      max_num_iter = props(21)


! ############################ Elasticity Matrix ############################

      ! Lame parameters
      ! Lame first parameter
      lambda = (E*nu)/((1.d0 + nu)*(1.d0 - 2.0d0*nu))  
      ! Lame second parameter  
      mu = E/(2.0d0*(1.d0 + nu))

      ! Elasticity matrix (C_mat)
      ! Linear elastic, homogeneous and isotropic material 
      ! Initialize the elasticity matrix to zero
      C_mat = 0.d0

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
! ####################### Loop through each integration point to perform computations #######################
! ###########################################################################################################  

      ! Global loop starting point -> loops through each integration point of the model
      ! The index of the integration points is i
      do i = 1, nblock
          
          ! ############################## Initial state ###############################   
          if (stateOld(i, 1) == 0.d0) then
              
              ! Compute Hooke's Law as a function of the strain increment
              ! Trace calculation (EPSILON_kk)
              trace_strainInc = sum(strainInc(i, 1:ndir))
            
              ! Calculate the direct components
              stressNew(i, 1:ndir) = stressOld(i, 1:ndir) +  
     1        lambda*trace_strainInc + 2.d0*mu*strainInc(i, 1:ndir)

              ! Calculate the shear components
              stressNew(i, ndir+1:ndir+nshr) = stressOld(i, ndir+1:ndir+nshr) +
     1        2.d0*mu*strainInc(i, ndir+1:ndir+nshr)  
              
              stateNew(i, 1) = 1.d0           ! Initiation flag
              stateNew(i, 2) = 0.d0           ! Equivalent plastic strain
              stateNew(i, 3) = Tr             ! Initial temperature
              stateNew(i, 4) = 0.d0           ! Yield stress
              stateNew(i, 5) = 0.d0           ! Plastic strain increment
              stateNew(i, 6) = 0.d0           ! Parameter D
              stateNew(i, 7) = 0.d0           ! Damage evolution parameter
              stateNew(i, 8) = 1              ! Total number of iterations
              stateNew(i, 9) = 1              ! Element status
              
              ! Store the stress tensor in SDVs for damage evolution computation
              do j = 1, ndir+nshr
                  stateNew(i, j + 9) = stressNew(i, j)
              end do

          ! ############################ From step 2 onward ############################    
          else

              eps = stateOld(i, 2)                                      ! Previous effective plastic strain
              Temp = stateOld(i, 3)                                     ! Temperature
              sigmaY = stateOld(i, 4)                                   ! Yield stress
              D = stateOld(i, 6)                                        ! Parameter D
              dmg_evol = stateOld(i, 7)                                 ! Damage evolution parameter
              trial_stress = 0.d0                                       ! Trial total stress with increment

              ! stress_old stores the stress tensor in the previous step
              if (D < 1.d0) then
                  stress_old(1:ndir+nshr) = stressOld(i, 1:ndir+nshr)       
              else
                  stress_old(1:ndir+nshr) = stateOld(i, 10:ndir+nshr+10)    
              end if
              
              ! Compute Hooke's Law as a function of the strain increment
              ! Calculates the stress trial increment tensor assuming a pure elastic response
              call elastic_stress(strainInc(i, 1:ndir+nshr), trial_stress,  
     1                            stress_old(1:ndir+nshr), lambda, mu, ndir, nshr) 

              ! Start the calculations to obtain the direction and magnitude of the
              ! plastic strain increment
              ! Direction = df/dSigma
              ! Magnitude = dLambda (plastic multiplier)
              ! dEpsilon^p = dLambda*df/dSigma = dp*3/2*DevStress/EquivalentStress
              call equivalent_stress(equiv_stress, trial_stress, 
     1                               dev_stress, ndir, nshr)

              ! Plastic strain direction
              if (equiv_stress /= 0) then
                  pl_strain_dir = (3.d0*dev_stress)/(2.d0*equiv_stress)
              else
                  pl_strain_dir = 0
              end if

              ! Elasticity matrix C times the plastic strain direction to calculate
              ! later the plastic corrector
              C_PSD = matmul(C_mat, pl_strain_dir)

              ! Initial value of the strain increment
              pl_strain_inc = 0.d0                                      ! Initial plastic strain increment
              pl_strain_inc_min = 0.d0                                  ! Minimum plastic strain increment
              pl_strain_inc_max = equiv_stress/(2.d0*mu)                ! Maximum plastic strain increment


              !##################################################################################################
              ! ######## Iteration to obtain the initial plastic strain increment (pl_strain_inc) starts ########
              ! Return mapping algorithm
              
              num_iter = 0    ! number of iterations
              
              do while (num_iter <= max_num_iter) 

                  ! Iteration control
                  num_iter = num_iter + 1
                  if (num_iter == max_num_iter) then
                      print*, 'ERROR - too many iterations | iter = ', num_iter
                      call XPLB_EXIT 
                  end if

                  ! Iteration effective plastic strain = effective plastic strain + increment
                  eps_iter = eps + pl_strain_inc  
                  ! Effective plastic strain increment rate  
                  eps_rate = pl_strain_inc/dtArray(1)
                  
                  ! Compute the stress using the plastic corrector for this iteration 
                  ! (corrected_stress_iter)
                  corrected_stress_iter = trial_stress - pl_strain_inc*C_PSD

                  ! Calculate the equivalent stress with the corrected stress in this
                  ! iteration
                  call equivalent_stress(equiv_stress, corrected_stress_iter, 
     1                                   dev_stress, ndir, nshr)

                  ! Calculate the Johnson-Cook equivalent stress
                  call johnson_cook_plasticity(eps_iter, A, B, n, m, 
     1                                         Tm, Tr, Temp, C,
     2                                         epsilon_dot_zero, 
     3                                         eps_rate, equiv_stress_jc)

                  ! If the calculated JC is less than the previous one, take the previous one
                  ! as equivalent yield stress. 
                  if (equiv_stress_jc < sigmaY) then
                      equiv_stress_jc = sigmaY
                  end if
                  
                  ! Compare the calculated JC stress with the equivalent stress resulting from
                  ! the corrected stress state. These two values should be as close as possible
                  ! within the specified tolerance to obtain the plastic strain increment for a
                  ! JC material
                  f = equiv_stress - equiv_stress_jc

                  if (abs(f) < tolerance) then
                      exit
                  end if

                  ! If there is no plastic strain increment and the JC yield stress is bigger
                  ! than the equivalent stress, the stress state is elastic
                  if ((pl_strain_inc == 0.d0) .and. (f < 0.d0)) then
                      exit
                  end if

                  ! Update the min and max plastic strain increment limits
                  if ((f >= 0.d0) .and. (pl_strain_inc >= pl_strain_inc_min)) then
                      pl_strain_inc_min = pl_strain_inc  
                  end if    

                  if ((f < 0.d0) .and. (pl_strain_inc < pl_strain_inc_max)) then
                      pl_strain_inc_max = pl_strain_inc  
                  end if

                  ! Update the plastic strain increment
                  pl_strain_inc = 0.5d0 * (pl_strain_inc_max + pl_strain_inc_min)

                  ! Update the pl_strain_inc to ensure continuity across increments and to 
                  ! improve the initial guess
                  if (num_iter == 1) then
                      pl_strain_inc = stateOld(i, 5)
                  end if

              end do

              ! Save the newly calculated stress
              stressNew(i, 1:ndir+nshr) = corrected_stress_iter(1:ndir+nshr)

              ! Store the stress tensor in SDVs for damage evolution computation
              do j = 1, ndir+nshr
                  stateNew(i, j + 9) = stressNew(i, j)
              end do

              ! Calculate the work increment
              dWork = dot_product(0.5d0 * (corrected_stress_iter(1:ndir+nshr) + stress_old(1:ndir+nshr)), 
     1                            strainInc(i, 1:ndir+nshr))

              ! Calculate the plastic work increment
              dPwork = 0.5d0 * pl_strain_inc * equiv_stress

              ! Calculate the internal energy per unit mass
              enerInternNew(i) = enerInternOld(i) + dWork/rho

              ! Calculate the dissipated inelastic energy per unit mass
              enerInelasNew(i) = enerInelasOld(i) + dPwork/rho

              ! Evaluate the equivalent strain at fracture
              call johnson_cook_damage(equiv_strain_fracture, equiv_stress, Tm, Tr, D1, D2,
     1                                 D3, D4, D5, Temp, epsilon_dot_zero, eps_rate, 
     2                                 corrected_stress_iter, ndir, nshr)
              
              ! Update the parameter D of the Johnson-Cook damage model
              if (equiv_strain_fracture /= 0.d0) then
                  D = D + abs(pl_strain_inc / equiv_strain_fracture)
              end if

              ! Fracture is allowed to occur when D = 1.0
              if (D >= 1.d0) then
                  D = 1.d0
                  dmg_evol = dmg_evol + 
     1            abs(pl_strain_inc * charLength(i) / pl_disp_failure)
                  if (dmg_evol >= 1.d0) then
                      dmg_evol = 1.d0
                  end if
                  stressNew(i, 1:ndir+nshr) = (1.d0 - dmg_evol)*stressNew(i, 1:ndir+nshr)
              end if

              ! Update the state variables
              stateNew(i, 1) = 1.d0                                                           ! Initiation flag	
              stateNew(i, 2) = eps_iter                                                       ! Equivalent plastic strain
              stateNew(i, 3) = Temp + beta*dPwork/rho/Cp                                      ! Temperature
              stateNew(i, 4) = equiv_stress_jc                                                ! Yield stress
              stateNew(i, 5) = pl_strain_inc                                                  ! Plastic strain increment in the step
              stateNew(i, 6) = D                                                              ! Parameter D
              stateNew(i, 7) = dmg_evol                                                       ! Damage evolution parameter
              stateNew(i, 8) = stateOld(i, 8) + num_iter                                      ! Total number of iterations

              if (dmg_evol >= 1.d0) then
                  stateNew(i, 9) = 0    ! Delete element
              else
                  stateNew(i, 9) = 1    ! Active element
              end if

          end if

      end do

      return
      end

! Additional subroutines
      include 'utils.for'
