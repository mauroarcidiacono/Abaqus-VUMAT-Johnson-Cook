      subroutine elastic_stress (strainInc, stressNew, stressOld,
     1 lambda, mu, ndir, nshr)
      ! This subroutine calculates the elastic stress tensor in Voigt
      ! notation using a strain tensor. The calculation is done according 
      ! to Hooke's law.
      ! SIGMA_ij = lambda*EPSILON_kk*KRONECKERDELTA_ij + 2*mu*EPSILON_ij

          integer ndir, nshr
          real*8 strainInc(ndir+nshr), stressNew(ndir+nshr),  
     1    stressOld(ndir+nshr), lambda, mu, trace_strainInc 
          
          ! Trace calculation (EPSILON_kk)
          trace_strainInc = sum(strainInc(1:ndir))
            
          ! Calculate the direct components
          stressNew(1:ndir) = stressOld(1:ndir) +  
     1    lambda*trace_strainInc + 2*mu*strainInc(1:ndir)

          ! Calculate the shear components
          stressNew(ndir+1:ndir+nshr) = stressOld(ndir+1:ndir+nshr) +
     1    2*mu*strainInc(ndir+1:ndir+nshr)
     
      return
      end


      subroutine equivalent_stress (equiv_stress, stress_tensor, 
     1 dev_stress, ndir, nshr)
      ! This subroutine calculates the equivalent Von Mises stress       
      ! from a Voigt notation stress tensor.

          integer ndir, nshr
          real*8 stress_tensor(ndir+nshr), hyd, 
     1    dev_stress(ndir+nshr), equiv_stress
          
          ! Calculate the hydrostatic component of the stress tensor
          hyd = sum(stress_tensor(1:ndir))/3.d0
          ! Calculate the deviatoric tensor
          dev_stress(1:ndir) = stress_tensor(1:ndir) - hyd
          dev_stress(ndir+1:ndir+nshr) = stress_tensor(ndir+1:ndir+nshr)

          ! Compute the equivalent Von Mises stress
          ! stressVM = sqrt(3/2 * dev_stress : dev_stress)
          equiv_stress = sqrt(3.d0/2.d0 *(dev_stress(1)**2.d0 + 
     1    dev_stress(2)**2.d0 + dev_stress(3)**2.d0 + 
     2    2.d0*dev_stress(4)**2.d0 + 2.d0*dev_stress(5)**2.d0 +     
     3    2.d0*dev_stress(6)**2.d0))

      return
      end


      subroutine johnson_cook_plasticity(eps_iter, A, B, n, m, Tm, Tr, 
     1 T, C, epsilon_dot_zero, eps_rate, equiv_stress_jc)
      ! This subroutine calculates the Von Mises stress using the 
      ! Johnson Cook equation. 

          real*8 eps_iter, A, B, n, m, Tm, Tr, T, C,  
     1    epsilon_dot_zero, eps_rate, homologous_Temp, 
     2    equiv_stress_jc  

          if (T < Tr) then
	       homologous_Temp = 0.d0
	   else if (T > Tm) then
	       homologous_Temp = 1.d0
	   else
		homologous_Temp = (T - Tr)/(Tm - Tr)
	   end if
	  
	   if (eps_rate == 0.d0) then
	       eps_rate = epsilon_dot_zero
	   end if
		
          equiv_stress_jc = (A + B*eps_iter**n)*
     1    (1 + C*log(eps_rate/epsilon_dot_zero))*
     2    (1 - homologous_Temp**m)

      return
      end


      subroutine johnson_cook_damage(equiv_strain_frac, equiv_stress,
     1 Tm, Tr, D1, D2, D3, D4, D5, T, epsilon_dot_zero, eps_rate,
     2 stress_tensor, ndir, nshr)
      ! This subroutine calculates the equivalent strain for fracture
      ! for the Johnson-Cook damage model. 

          integer ndir, nshr
          real*8 D, Tm, Tr, D1, D2, D3, D4, D5, T,  
     1    epsilon_dot_zero, eps_rate, homologous_Temp, 
     2    pressure_stress_ratio, equiv_strain_frac,
     3    stress_tensor(ndir+nshr), hyd, equiv_stress

          if (T < Tr) then
	       homologous_Temp = 0.d0
	   else if (T > Tm) then
	       homologous_Temp = 1.d0
	   else
	       homologous_Temp = (T - Tr)/(Tm - Tr)
	   end if
	  
	   if (eps_rate == 0.d0) then
	       eps_rate = epsilon_dot_zero
	   end if

          ! Calculate the hydrostatic component of the stress tensor
          hyd = sum(stress_tensor(1:ndir))/3.d0

          ! Calculate the pressure stress ratio
          if (equiv_stress == 0) then
              pressure_stress_ratio = 0
          else
              pressure_stress_ratio = hyd/equiv_stress
          end if
		
          equiv_strain_frac = (D1 + 
     1    D2*exp(-D3*pressure_stress_ratio))*
     2    (1 + D4*log(eps_rate/epsilon_dot_zero))*
     3    (1 + D5*homologous_Temp)

      return
      end


