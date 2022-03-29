	program main
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 24
	integer, parameter :: nl = 3, nb = 2, nm = 3 
	integer, parameter :: nt = 474
      integer, parameter :: iiniter = 1000, initer = 20000
      integer, parameter :: niter = 100000, niter1 = 20000
	
	real*8 yold, yoldstd, time
	
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, zeta(ni), sigmaz
      real*8 sigmai(ni), sigmay
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl), omega(nb,nb,nl)
      real*8 gi(ni), theta(nl)
      real*8 sigmaa, sigmab
      
      real*8 sigma_mu, sigma_theta, d0, s0
      real*8 a0, b0, a1, b1, a2, b2, a3, b3, a4, b4
            
      integer ngi(nl)
      real*8 tijA(nt), tijB(nt)
      real*8 weight(nl), sumw, ystar, api, sumh
      real*8 temp1, temp2, temp3, temp4, temp5
      real*8 harmon(nt,nl), Eyij(ni,nt)
      
      real*8 cmaxir(ni,niter1), cpoir(ni,niter1)
      real*8 probi(ni,nl), cmaxi(ni), cpoi(ni)
      real*8 prob(ni,nl), alpml
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vzeta/zeta
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      common /vsigmay/sigmay
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      common /vomega/omega
      
      common /vgi/gi
      common /vtheta/theta

      common /vsigmaa/sigmaa
      common /vsigmab/sigmab
      
      common /vsigma_mu/sigma_mu
      common /vsigma_theta/sigma_theta
      common /vd0/d0
      common /vs0/s0
      common /va0/a0
      common /vb0/b0
      common /va1/a1
      common /vb1/b1
      common /va2/a2
      common /vb2/b2
      common /va3/a3
      common /vb3/b3
      common /va4/a4
      common /vb4/b4
                 
      common /vprobi/probi
      common /vcmaxi/cmaxi
      common /vcpoi/cpoi
      
    	open(unit =  5, file = 'CircadianDataP1.txt') 
    	open(unit =  6, file = 'CircadianDataP1_Time.txt') 
    	open(unit =  7, file = 'CircadianP1_Initial.txt')

    	open(unit = 11, file = 'CircadianP1_Output1.txt') 
    	open(unit = 12, file = 'CircadianP1_Output2.txt') 
    	open(unit = 13, file = 'CircadianP1_Output3.txt') 
    	open(unit = 14, file = 'CircadianP1_Output4.txt') 
    	open(unit = 15, file = 'CircadianP1_Output5.txt') 
    	open(unit = 16, file = 'CircadianP1_Output6.txt') 
    	open(unit = 17, file = 'CircadianP1_Output7.txt') 
    	open(unit = 18, file = 'CircadianP1_Output8.txt') 
    	open(unit = 19, file = 'CircadianP1_Output9.txt') 
      
    	open(unit = 21, file = 'CircadianP1_Output10.txt') 
    	open(unit = 22, file = 'CircadianP1_Output11.txt') 
    	open(unit = 23, file = 'CircadianP1_Output12.txt') 
    	open(unit = 24, file = 'CircadianP1_LPML.txt') 

      iseed = 9999999

      do i = 1, ni
          do j = 1, nj
              yij(i,j) = -999.d0
          enddo
      enddo

	do ii = 1, 15912
		
	    read(5,*) i, j, yold, time
                		  
	    yij(i,j) = dlog(yold)
        
          if ( (time .gt. 24.d0) 
     +         .and. (time .le. 48.d0) ) then
              
              time = time - 24.d0
              
          else if ( (time .gt. 48.d0) 
     +         .and. (time .le. 72.d0) ) then
     
              time = time - 48.d0
              
          else if ( (time .gt. 72.d0) 
     +         .and. (time .le. 96.d0) ) then
     
              time = time - 72.d0
              
          endif    
                
	    tij(i,j) = time/24.d0
          
      enddo  
      
      do i = 1, ni
          ij(i) = 0
          do j = 1, nj
              if (yij(i,j) .ne. -999.d0) then
                  ij(i) = ij(i) + 1
              endif
          enddo
      enddo   
            
	do jj = 1, nt
		
	    read(6,*) j, time
                		  
          tijA(j) = time
        
          if ( (time .gt. 24.d0) 
     +         .and. (time .le. 48.d0) ) then
              
              time = time - 24.d0
              
          else if ( (time .gt. 48.d0) 
     +         .and. (time .le. 72.d0) ) then
     
              time = time - 48.d0
              
          else if ( (time .gt. 72.d0) 
     +         .and. (time .le. 96.d0) ) then
     
              time = time - 72.d0
              
          endif    
                
	    tijB(j) = time/24.d0
                          
      enddo  
          
c     Hyperparameter for pior distribution  
                  
      sigma_mu = 100.d0; sigma_theta = 100.d0
      d0 = dfloat(nb) + 0.1d0; s0 = 0.1d0      
      a0 = 2.0d0; b0 = 1.0d0  
      a1 = 1.0d0; b1 = 0.1d0        
      a2 = 2.0d0; b2 = 1.0d0        
      a3 = 2.0d0; b3 = 1.0d0        
      a4 = 2.0d0; b4 = 1.0d0        
                 
      api = dconst('PI')
            
c     Initial values    
             
      read(7,*) mu
      read(7,*) sigmaz 
      read(7,*) sigmay
      read(7,*) sigmai
      read(7,*) sigmaa
      read(7,*) sigmab
      read(7,*) betam
      read(7,*) alpham
      read(7,*) omega
      read(7,*) theta
                          
      do i = 1, ni
                              
          zeta(i) = 0.d0
          
          do k = 1, nb
              do l = 1, nl
                  bi(i,k,l) = 0.d0
              enddo
          enddo
          
          gi(i) = 1.d0
                 
      enddo
                      
c     Run MCMC    
      
      call rnset(iseed)   
                 
      do ir = 1, iiniter
          
          call gibbs_gi(iseed)          
          call gibbs_bi2(iseed)          
          call gibbs_bi1(iseed)          
          call gibbs_zeta(iseed)          
     
          do l = 1, nl
              ngi(l) = 0
              do i = 1, ni
                  if (gi(i) .eq. dfloat(l)) then
                      ngi(l) = ngi(l) + 1
                  endif
              enddo
          enddo
          
          sumw = 0.d0
          do l = 1, nl              
              sumw = sumw + dexp(theta(l))
          enddo
          do l = 1, nl              
              weight(l) = dexp(theta(l))/sumw
          enddo
                    
      enddo
      
      icount = 0
      do ir = 1, initer
                        
          call gibbs(iseed)
          
          do l = 1, nl
              ngi(l) = 0
              do i = 1, ni
                  if (gi(i) .eq. dfloat(l)) then
                      ngi(l) = ngi(l) + 1
                  endif
              enddo
          enddo
          
          sumw = 0.d0
          do l = 1, nl              
              sumw = sumw + dexp(theta(l))
          enddo
          do l = 1, nl              
              weight(l) = dexp(theta(l))/sumw
          enddo
                                                
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1

              write(11,1) icount, ngi
              write(12,1) icount, nint(gi)
              write(13,3) icount, mu, sigmaz, sigmay 
              write(14,3) icount, betam
              write(15,3) icount, alpham
              write(16,3) icount, omega
              write(17,3) icount, theta
              write(18,3) icount, sigmaa, sigmab 
              write(19,3) icount, sigmai
              
              call Gen_Prob_CPO(iseed)  
              
              do i = 1, ni
                  cmaxir(i,icount) = -cmaxi(i)
                  cpoir(i,icount) = cpoi(i)
              enddo                                   
                            
          endif
                                   
      enddo  
               
      do i = 1, ni
          cmaxi(i) = cmaxir(i,1)    
          do ir = 2, icount
              if (cmaxi(i) .lt. cmaxir(i,ir)) then
                  cmaxi(i) = cmaxir(i,ir)
              endif
          enddo  
      enddo
      
      do i = 1, ni
          cpoi(i) = 0.d0
      enddo    
      do ir = 1, icount      
          do i = 1, ni
              cpoi(i) = cpoi(i) 
     +                + dexp(cmaxir(i,ir) - cmaxi(i))/cpoir(i,ir)
          enddo                 
      enddo
      
      alpml = 0.0d0
      do i = 1, ni    
          alpml = alpml 
     +          + dlog(dfloat(icount)) - cmaxi(i) - dlog(cpoi(i))
      enddo
            
      do j = 1, nt
          do l = 1, nl
              harmon(j,l) = 0.d0
          enddo    
      enddo
      do i = 1, ni
          do j = 1, nt
              Eyij(i,j) = 0.d0
          enddo
          do l = 1, nl
              prob(i,l) = 0.d0
          enddo
      enddo
          
      icount = 0
      do ir = 1, niter
              
          call gibbs(iseed)

          do l = 1, nl
              ngi(l) = 0
              do i = 1, ni
                  if (gi(i) .eq. dfloat(l)) then
                      ngi(l) = ngi(l) + 1
                  endif
              enddo
          enddo
          
          sumw = 0.d0
          do l = 1, nl              
              sumw = sumw + dexp(theta(l))
          enddo
          do l = 1, nl              
              weight(l) = dexp(theta(l))/sumw
          enddo
                    
          if (ir - ir/5*5 .eq. 1) then
        
              icount = icount + 1

              write(11,1) icount, ngi
              write(12,1) icount, nint(gi)
              write(13,3) icount, mu, sigmaz, sigmay 
              write(14,3) icount, betam
              write(15,3) icount, alpham
              write(16,3) icount, omega
              write(17,3) icount, theta
              write(18,3) icount, sigmaa, sigmab 
              write(19,3) icount, sigmai
              
              do l = 2, nl
                  do j = 1, nt
                      
                      sumh = 0.d0
                      do m = 1, nm
                          
                          temp1 = alpham(m,l)
                  
                          temp2 = dexp(temp1)
     +                            /(1.d0 + dexp(temp1))
                  
                          temp3 = 2.d0*api*dfloat(m)
     +                            *(tijB(j) + temp2)
                  
                          temp4 = betam(m,l) 
                      
                          temp5 = dexp(temp4)
                  
                          sumh = sumh + temp5*dcos(temp3)
                      
                      enddo
                                                                      
                      harmon(j,l) = harmon(j,l) 
     +                            + dexp(sumh)/dfloat(niter1)
                      
                  enddo
              enddo
                            
              do i = 1, ni
                              
                  l = nint(gi(i))
                  
                  if (l .eq. 1) then

                      do j = 1, nt
                      
                          ystar = mu + zeta(i)
                                            
                          Eyij(i,j) = Eyij(i,j) 
     +                              + dexp(ystar)/dfloat(niter1)
                          
                      enddo
                              
                  else 
                      
                      do j = 1, nt
                      
                          sumh = 0.d0
                          do m = 1, nm
                          
                              temp1 = alpham(m,l) + bi(i,2,l)
                  
                              temp2 = dexp(temp1)
     +                                /(1.d0 + dexp(temp1))
                  
                              temp3 = 2.d0*api*dfloat(m)
     +                                *(tijB(j) + temp2)
                  
                              temp4 = betam(m,l) + bi(i,1,l)
                      
                              temp5 = dexp(temp4)
                  
                              sumh = sumh + temp5*dcos(temp3)
                      
                          enddo
                          
                          ystar = mu + zeta(i) + sumh
                                            
                          Eyij(i,j) = Eyij(i,j) 
     +                              + dexp(ystar)/dfloat(niter1)
                          
                      enddo
                                            
                  endif
                  
              enddo
                       
              call Gen_Prob_CPO(iseed)  
              
              do i = 1, ni
                                    
                  do l = 1, nl
                      prob(i,l) = prob(i,l) 
     +                          + probi(i,l)/dfloat(niter1)
                  enddo

                  cmaxir(i,icount) = -cmaxi(i)
                  cpoir(i,icount) = cpoi(i)
                  
              enddo
                                   
          endif
          
      enddo    	  
                
      do j = 1, nt
          write(21,3) j, tijA(j), (harmon(j,l),l=1,nl)
      enddo
      
      do i = 1, ni
          do j = 1, nt
              write(22,5) i, j, tijA(j), Eyij(i,j)
          enddo
      enddo
      
      do i = 1, ni
          write(23,3) i, (prob(i,l),l=1,nl)
      enddo
      
      do i = 1, ni
          cmaxi(i) = cmaxir(i,1)    
          do ir = 2, icount
              if (cmaxi(i) .lt. cmaxir(i,ir)) then
                  cmaxi(i) = cmaxir(i,ir)
              endif
          enddo  
      enddo
      
      do i = 1, ni
          cpoi(i) = 0.d0
      enddo    
      do ir = 1, icount      
          do i = 1, ni
              cpoi(i) = cpoi(i) 
     +                + dexp(cmaxir(i,ir) - cmaxi(i))/cpoir(i,ir)
          enddo                 
      enddo
      
      alpml = 0.0d0
      do i = 1, ni    
          alpml = alpml 
     +          + dlog(dfloat(icount)) - cmaxi(i) - dlog(cpoi(i))
      enddo
          
      write(24,4) alpml
      
    1 format(1000i5)
    2 format(2i5, 1000f10.5)
    3 format(i5, 1000f20.10)
    4 format(1000f20.10)
    5 format(2i5, 1000f20.10)
            
      end program
           	            
      
	subroutine Gen_Prob_CPO(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 24
	integer, parameter :: nl = 3, nb = 2, nm = 3 
	integer, parameter :: nr = 1000 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 omega(nb,nb,nl)
      real*8 gi(ni), theta(nl)
      real*8 probi(ni,nl), cmaxi(ni), cpoi(ni)

      real*8 omegal(nb,nb)
      real*8 ystar, api, sumh
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp 
      real*8 weight(nl), sumw     
      real*8 bstar(nb), tol, rsig(nb,nb)
      real*8 suma, sumb, pdfIR(nr)
      real*8 prob(nl), cpo(nl)
      real*8 amax, summ, pdf
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vomega/omega
      
      common /vgi/gi
      common /vtheta/theta
      
      common /vprobi/probi
      common /vcmaxi/cmaxi
      common /vcpoi/cpoi
      
      external dconst, dlinrg, drnunf, dchfac, drnmvn

      api = dconst('pi')
                           
      sumw = 0.d0
      do l = 1, nl
          sumw = sumw + dexp(theta(l))
      enddo
      do l = 1, nl              
          weight(l) = dexp(theta(l))/sumw
      enddo
            
      do i = 1, ni
                          
          do l = 1, nl
                                                                
              if (l .eq. 1) then
                        
                  suma = 0.d0; sumb = 0.d0
                  do j = 1, ij(i)
                                      
                      ystar = yij(i,j) - mu 

                      suma = suma - ystar**2/(2.d0*sigmai(i))

                      sumb = sumb + ystar/sigmai(i)
                                                            
                  enddo
      
                  temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
                  
                  pdf = dlog(weight(l))
     +                - dfloat(ij(i))
     +                  *dlog(2.d0*api*sigmai(i))/2.d0
     +                + suma 
     +                - dlog(sigmaz)/2.d0
     +                - dlog(temp)/2.d0
     +                + sumb**2/(2.d0*temp)
                                    
              else
                  
                  do j1 = 1, nb
                      do j2 = 1, nb
                          omegal(j1,j2) = omega(j1,j2,l)
                      enddo
                  enddo
                  
	            tol = 100.d0*dmach(4)
	            call dchfac(nb,omegal,nb,tol,irank,rsig,nb)

                  do ir = 1, nr
                  
	                call rnset(iseed)
	                call drnmvn(1,nb,rsig,nb,bstar,1)
	                call rnget(iseed)      
                                    
                      suma = 0.d0; sumb = 0.d0
                      do j = 1, ij(i)
                            
                          sumh = 0.d0                        
                          do m = 1, nm
                          
                              temp1 = alpham(m,l) + bstar(2)
                  
                              temp2 = dexp(temp1)
     +                                /(1.d0 + dexp(temp1))
                  
                              temp3 = 2.d0*api*dfloat(m)
     +                                *(tij(i,j) + temp2)
                  
                              temp4 = betam(m,l) + bstar(1)
                      
                              temp5 = dexp(temp4)
                  
                              sumh = sumh + temp5*dcos(temp3)
                      
                          enddo
                      
                          ystar = yij(i,j) - mu - sumh

                          suma = suma - ystar**2/(2.d0*sigmai(i))

                          sumb = sumb + ystar/sigmai(i)
                                                            
                      enddo
                         
                      temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
                                        
                      pdfIR(ir) = dlog(weight(l))
     +                          - dfloat(ij(i))
     +                            *dlog(2.d0*api*sigmai(i))/2.d0
     +                          + suma 
     +                          - dlog(sigmaz)/2.d0
     +                          - dlog(temp)/2.d0
     +                          + sumb**2/(2.d0*temp)
                  
                  enddo
                                    
                  amax = pdfIR(1)
                  do ir = 2, nr
                      if (amax .lt. pdfIR(ir)) then
                          amax = pdfIR(ir)
                      endif
                  enddo
                  
                  summ = 0.d0
                  do ir = 1, nr                  
                      summ = summ 
     +                     + dexp(pdfIR(ir) - amax)/dfloat(nr)
                  enddo
                  
                  pdf = amax + dlog(summ)                 
                                        
              endif
                                      
              prob(l) = pdf
              cpo(l) = pdf
                            
          enddo
          
          amax = prob(1)
          do l = 2, nl
              if (amax .lt. prob(l)) then
                  amax = prob(l)
              endif
          enddo       
          
          summ = 0.d0
          do l = 1, nl
              prob(l) = dexp(prob(l)-amax)
              summ = summ + prob(l)
          enddo
          
          do l = 1, nl
              probi(i,l) = prob(l)/summ
          enddo
                                        
          cmaxi(i) = cpo(1)
          do l = 2, nl
              if (cmaxi(i) .lt. cpo(l)) then
                  cmaxi(i) = cpo(l)
              endif
          enddo       
          
          cpoi(i) = 0.d0
          do l = 1, nl
              cpoi(i) = cpoi(i) + dexp(cpo(l) - cmaxi(i))
          enddo
                              
      enddo     
                       
      end subroutine	   
      
      include 'CircadianP1_Gibbs.f'
	include 'optim1.f'
            