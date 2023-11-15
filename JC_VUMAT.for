! #######################################################################
! User subroutine to model the Jonhson-Cook plasticity, damage and the 
! Taylor-Quinney conversion of mechanical work into heat during plastic
! deformation.
!
! Abaqus version: Abaqus 2022
! Intel oneAPI Compiler 2022.0.2 intel64
! Visual Studio 2019
!
! Author: Mauro Francisco Arcidiacono
! #######################################################################

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
      integer num_iter
      real*8 E, nu, A, B, n, m, Tm, Tr, C, epsilon_dot_zero, D1, D2,
    1 D3, D4, D5, beta, Cp, unit_conversion_factor, convergence_tolerance_factor,
    2 tolerance, mu, lambda 
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
      lambda = (E*nu)/((1 + nu)*(1 - 2*nu))  
      ! Lame second parameter  
      mu = E/(2*(1 + nu))

      
      
      


      return
      end