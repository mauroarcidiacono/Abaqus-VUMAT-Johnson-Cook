      subroutine elastic_stress (strainInc, stressNew, stressOld,
     1 lambda, mu, ndir, nshr)
       ! This subroutine calculates the elastic stress using a strain
       ! tensor. The calculation is done according to Hooke's law.
       ! SIGMA_ij = lambda*EPSILON_kk*KRONECKERDELTA_ij + 2*mu*EPSILON_ij

            include 'vaba_param.inc'
            integer ndir, nshr
            real*8 strainInc(ndir+nshr), stressNew(ndir+nshr),  
     1      stressOld(ndir+nshr), lambda, mu 
          
            ! Trace calculation (EPSILON_kk)
            trace_strainInc = sum(strainInc(1:ndir))
            
            ! Calculate the direct components
            stressNew(1:ndir) = stressOld(1:ndir) + lambda*trace_strainInc + 
     1      2*mu*strainInc(1:ndir)

            ! Calculate the shear components
            stressNew(ndir+1:ndir+nshr) = stressOld(ndir+1:ndir+nshr) +
     1      2*mu*strainInc(ndir+1:ndir+nshr)
     
            return
            end