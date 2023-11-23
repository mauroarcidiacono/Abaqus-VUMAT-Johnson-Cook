      subroutine elastic_stress (strainInc, stressNew, stressOld,
     1 lambda, mu, ndir, nshr)
       ! This subroutine calculates the elastic stress tensor in Voigt
       ! notation using a strain tensor. The calculation is done according 
       ! to Hooke's law.
       ! SIGMA_ij = lambda*EPSILON_kk*KRONECKERDELTA_ij + 2*mu*EPSILON_ij

            integer ndir, nshr
            real*8 strainInc(ndir+nshr), stressNew(ndir+nshr),  
     1      stressOld(ndir+nshr), lambda, mu 
          
            ! Trace calculation (EPSILON_kk)
            trace_strainInc = sum(strainInc(1:ndir))
            
            ! Calculate the direct components
            stressNew(1:ndir) = stressOld(1:ndir) +  
     1      lambda*trace_strainInc + 2*mu*strainInc(1:ndir)

            ! Calculate the shear components
            stressNew(ndir+1:ndir+nshr) = stressOld(ndir+1:ndir+nshr) +
     1      2*mu*strainInc(ndir+1:ndir+nshr)
     
      return
      end


      subroutine equivalent_stress (equiv_stress, stress_tensor, 
     1 dev_stress, ndir, nshr)
       ! This subroutine calculates the equivalent stress from a Voigt     
       ! notation stress tensor.

            integer ndir, nshr
            real*8 stress_tensor(ndir+nshr), hyd, dev_stress(ndir+nshr),
     1      equiv_stress
          
            ! Calculate the hydrostatic component of the stress tensor
            hyd = sum(stress_tensor(1:ndir))/3.d0
            ! Calculate the deviatoric tensor
            dev_stress(1:ndir) = stress_tensor(1:ndir) - hyd
            dev_stress(ndir+1:ndir+nshr) = stress_tensor(ndir+1:ndir+nshr)

            ! Compute the equivalent stress
            equiv_stress = sqrt(3.d0/2.d0 *(dev_stress(1)**2.d0 + dev_stress(2)**2.d0 +
     1      dev_stress(3)**2.d0 + 2.d0*dev_stress(4)**2.d0 + 2.d0*dev_stress(5)**2.d0 +
     2      2.d0*dev_stress(6)**2.d0 ))
     
      return
      end


      subroutine johnson_cook(eps_iter, A, B, n, m, Tm, Tr, C, 
     1 epsilon_dot_zero, eps_rate, equiv_stress_jc)
      ! This subroutine calculates the Von Mises stress using the 
      ! Johnson Cook equation. 

      real*8 A, B, n, m, Tm, Tr, C, epsilon_dot_zero, eps_rate, 
     1 equiv_stress_jc, T_  

      equiv_stress = (A + B*eps_iter**n)*
     1 (1 + C*log(eps_rate/epsilon_dot_zero))*
     2 (1 - )

      return
      end

            