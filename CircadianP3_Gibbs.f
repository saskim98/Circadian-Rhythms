      subroutine gibbs(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
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
                       
      call gibbs_gi(iseed)          

      
      call gibbs_mu(iseed)          
      
      
      call gibbs_bi2(iseed)          
      
      call gibbs_alpham(iseed) 
      
	call gibbs_sigmaa(iseed)  

      
      call gibbs_bi1(iseed)          

      call gibbs_betam(iseed)         

      call gibbs_theta(iseed)          
      
	call gibbs_sigmab(iseed)    

            
      call gibbs_sigmaz(iseed)          
            
      call gibbs_zeta(iseed)          
      
            
      call gibbs_sigmay(iseed)          
      
      call gibbs_sigmai(iseed)          
            
      call gibbs_omega(iseed)          
      
      end subroutine
          
      
	subroutine gibbs_gi(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl), omega(nb,nb,nl)
      real*8 gi(ni), theta(nl)

      real*8 omegal(nb,nb), omegai(nb,nb)
      real*8 ystar, bstar(nb), api, sumh
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp 
      real*8 weight(nl), sumw     
      real*8 suma, sumb, bsb
      
      integer ipvt(nb)
      real*8 fac(nb,nb), det1, det2, det      
      
      real*8 prob(nl), cprob(nl)
      real*8 amax, u, summ, pdf
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      common /vomega/omega
      
      common /vgi/gi
      common /vtheta/theta
            
      external dconst, dlinrg, dlftrg, dlfdrg, dqdagi, drnunf

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
                  
                  do jj = 1, nb
                      bstar(jj) = bi(i,jj,l)
                  enddo
                      
                  do j1 = 1, nb
                      do j2 = 1, nb
                          omegal(j1,j2) = omega(j1,j2,l)
                      enddo
                  enddo

                  call dlinrg(nb, omegal, nb, omegai, nb)
      
                  call dlftrg(nb, omegal, nb, fac, nb, ipvt)
                  call dlfdrg(nb, fac, nb, ipvt, det1, det2)
                  det = det1*10.d0**det2                  
                                 
                  bsb = dblinf(nb,nb,omegai,nb,bstar,bstar)
                                                      
                  suma = 0.d0; sumb = 0.d0
                  do j = 1, ij(i)
                            
                      sumh = 0.d0                        
                      do m = 1, nm
                          
                          temp1 = alpham(m,l) + bstar(2)
                  
                          temp2 = dexp(temp1)
     +                            /(1.d0 + dexp(temp1))
                  
                          temp3 = 2.d0*api*dfloat(m)
     +                            *(tij(i,j) + temp2)
                  
                          temp4 = betam(m,l) + bstar(1)
                      
                          temp5 = dexp(temp4)
                  
                          sumh = sumh + temp5*dcos(temp3)
                      
                      enddo
                      
                      ystar = yij(i,j) - mu - sumh

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
     +                - dfloat(nb)*dlog(2.d0*api)/2.d0
     +                - dlog(det)/2.d0  
     +                - bsb/2.d0
                  
              endif
                                      
              prob(l) = pdf
                            
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
              prob(l) = prob(l)/summ
          enddo
        
          cprob(1) = prob(1)
          do l = 2, nl
              cprob(l) = cprob(l-1) + prob(l)
          enddo

		call rnset(iseed)
		u = drnunf()
		call rnget(iseed)
                		                        
          if (u .le. cprob(1)) then
              gi(i) = 1.d0
          endif
          do l = 2, nl
              if ((u .gt. cprob(l-1)) 
     +             .and. (u .le. cprob(l))) then
                  gi(i) = dfloat(l)
              endif
          enddo
                                              
      enddo     
                       
      end subroutine	   
      
      
	subroutine gibbs_theta(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
	integer, parameter :: npoint = 5, cpoint = 3, ndim = 3
		      
      real*8 betam(nm,nl)
      real*8 gi(ni), theta(nl)
      real*8 sigma_theta, a4, b4
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 epsilon, astar, error
      real*8 yystar(npoint), xxstar(npoint,ndim)
      real*8 xx(ndim,ndim), xy(ndim)
      real*8 xxi(ndim,ndim), xhat(ndim)
      
      integer ldum
      
      common /vbetam/betam

      common /vgi/gi
      common /vtheta/theta
            
      common /vsigma_theta/sigma_theta
      common /va4/a4
      common /vb4/b4
      
      common /vldum/ldum
      
      external ftheta, drnnof, drnunf
        
      epsilon = 0.1d0
      
      do l = 2, nl

          ldum = l
        
          if (l .eq. 2) then
              
              temp = theta(2) - theta(3)
              
              bold = dlog(temp)
              
          else if (l .eq. nl) then
          
              temp = theta(nl-1) - theta(nl)

              bold = dlog(temp)
              
          else
              
              temp = (theta(l-1) - theta(l+1))
     +               /(theta(l) - theta(l+1))

              bold = dlog(dlog(temp))
              
          endif
                    
          nopt = 1
          reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
          step(1) = 0.2d0 ; start(1) = bold
          call nelmin(ftheta, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
          amean = xmin(1)
          
  11      do kk = 1, npoint
                  
              error = dfloat(kk - cpoint)*0.2d0
              astar = amean + error*epsilon
                  
              yystar(kk) = ftheta(astar)
              xxstar(kk,1) = astar**2
              xxstar(kk,2) = astar
              xxstar(kk,3) = 1.d0
                  
          enddo
              
          do kk = 1, npoint
              if (kk .ne. cpoint) then                  
                  if (yystar(kk) .le. ynewlo) then
                      epsilon = epsilon*1.2d0
                      goto 11
                  endif
              endif
          enddo
        
          do kk1 = 1, ndim
              do kk2 = 1, ndim
                  xx(kk1,kk2) = 0.d0
                  do kk = 1, npoint
                      xx(kk1,kk2) = xx(kk1,kk2) 
     +                            + xxstar(kk,kk1)
     +                              *xxstar(kk,kk2)
                  enddo
              enddo
              xy(kk1) = 0.d0
              do kk2 = 1, npoint
                  xy(kk1) = xy(kk1) + xxstar(kk2,kk1)
     +                                *yystar(kk2)
              enddo
          enddo

          call dlinrg(ndim, xx, ndim, xxi, ndim)	
                        
          do kk1 = 1, ndim
              xhat(kk1) = 0.d0
              do kk2 = 1, ndim
                  xhat(kk1) = xhat(kk1)
     +                      + xxi(kk1,kk2)*xy(kk2)
              enddo
          enddo
                                                       
          asigma = 1.0d0/(xhat(1)*2.d0)*1.5d0            
          if (asigma .lt. 0.0d0) asigma = -asigma

          bpdf = -ftheta(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 20
	  
              call rnset(iseed)
	        rv = drnnof()
	        call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -ftheta(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
		    ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
		        bpdf = apdf
		        bold = anew
 	        else
		        call rnset(iseed)
		        u = drnunf()
		        call rnget(iseed)
		        if (dlog(u) .le. ratio) then
		            bpdf = apdf
		            bold = anew
 		        endif
		    endif
		
          enddo     

          if (l .eq. 2) then
              
              theta(2) = theta(3) + dexp(bold)
              
          else if (l .eq. nl) then
          
              theta(nl) = theta(nl-1) - dexp(bold)
              
          else
              
              theta(l) = theta(l+1) 
     +                 + (theta(l-1) - theta(l+1))
     +                   *dexp(-dexp(bold))
              
          endif
          
      enddo
      
      end subroutine
            
      
      real*8 function ftheta(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		      
      real*8 betam(nm,nl)
      real*8 gi(ni), theta(nl)
      real*8 sigma_theta, a4, b4
            
      integer ldum      
      real*8 star, weight(nl)
      real*8 sump, sumw, temp
      real*8 suma, sumb, sumc
      real*8 summ, pdf
      
      common /vbetam/betam

      common /vgi/gi
      common /vtheta/theta
            
      common /vsigma_theta/sigma_theta
      common /va4/a4
      common /vb4/b4
      
      common /vldum/ldum
      
      l = ldum
      
      if (l .eq. 2) then
              
          theta(2) = theta(3) + dexp(star)
              
      else if (l .eq. nl) then
          
          theta(nl) = theta(nl-1) - dexp(star)
              
      else
              
          theta(l) = theta(l+1) 
     +             + (theta(l-1) - theta(l+1))
     +                *dexp(-dexp(star))
              
      endif
            
      sumw = 0.d0
      do ll = 1, nl
          sumw = sumw + dexp(theta(ll))
      enddo
      do ll = 1, nl              
          weight(ll) = dexp(theta(ll))/sumw
      enddo
      
      sump = 0.d0
      do i = 1, ni
          
          ll = nint(gi(i)) 
          
          if (ll .eq. l) then
              sump = sump + theta(ll)
          endif
          
          sump = sump - dlog(sumw)
          
      enddo
                  
      sumw = 0.d0; suma = a4; sumb = b4
      do ll = 2, nl
          do mm = 1, nm

              sumw = sumw - dlog(weight(ll))/2.d0
              
              suma = suma + 1.d0/2.d0
              
              sumb = sumb + betam(mm,ll)**2
     +                      /(2.d0*weight(ll))
              
          enddo
      enddo
      
      sumc = -theta(l)**2/(2.d0*sigma_theta)
      
      summ = sump + sumw - suma*dlog(sumb) + sumc
      
      if (l .eq. 2) then
              
          pdf = summ + star
                        
      else if (l .eq. nl) then
          
          pdf = summ + star 
              
      else
          
          pdf = summ + star - dexp(star)
              
      endif
      
      ftheta = -pdf
     
      end function
      
      
	subroutine gibbs_betam(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni), theta(nl)
      
      real*8 a3, b3, a4, b4
            
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      integer ldum, mdum
      real*8 weight(nl), sumw      
      real*8 ystar, api, sumh, sumhm
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp     
      real*8 suma, sumb, sumc, sumd     
      real*8 sum1, sum2, sum3, sum4, summ
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
      common /vtheta/theta
                  
      common /va3/a3
      common /vb3/b3
      common /va4/a4
      common /vb4/b4
      
      common /vldum/ldum
      common /vmdum/mdum
      common /vweight/weight
      
      external fbetam, drnnof, drnunf
           
      api = dconst('PI')
      
      sumw = 0.d0
      do l = 1, nl
          sumw = sumw + dexp(theta(l))
      enddo
      do l = 1, nl              
          weight(l) = dexp(theta(l))/sumw
      enddo
      
      do m = 1, nm
          betam(m,1) = 0.d0
      enddo                      
      do l = 2, nl
          
          ldum = l
                
          do m = 1, nm

              mdum = m
          
              bold = betam(m,l)

              nopt = 1
              reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
              step(1) = 0.2d0 ; start(1) = bold
              call nelmin(fbetam, nopt, start, xmin, ynewlo,  
     +                    reqmin, step, konvge, kcount, 
     +                    icount, numres, ifault)
              amean = xmin(1)
              betam(m,l) = amean
                          
              summ = 0.d0
              do i = 1, ni   
      
                  ll = nint(gi(i))

                  if (ll .eq. l) then
          
                  suma = 0.d0; sumb = 0.d0
                  sumc = 0.d0; sumd = 0.d0
                  do j = 1, ij(i)
                    
                      sumh = 0.d0; sumhm = 0.d0
                      do mm = 1, nm
                          
                          temp1 = alpham(mm,ll) + bi(i,2,ll)
                  
                          temp2 = dexp(temp1)
     +                            /(1.d0 + dexp(temp1))
                  
                          temp3 = 2.d0*api*dfloat(mm)
     +                            *(tij(i,j) + temp2)
                  
                          temp4 = betam(mm,ll) + bi(i,1,ll)
                      
                          temp5 = dexp(temp4)
                  
                          sumh = sumh + temp5*dcos(temp3)
                      
                          if (mm .eq. m) then
                              sumhm = temp5*dcos(temp3)
                          endif
                          
                      enddo
                                                        
                      ystar = yij(i,j) - mu - sumh
              
                      suma = suma + sumhm**2/sigmai(i)

                      sumb = sumb + ystar*sumhm/sigmai(i)

                      sumc = sumc + sumhm/sigmai(i)

                      sumd = sumd + ystar/sigmai(i)
          
                  enddo
                                
                  temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
                
                  summ = summ 
     +                 - suma + sumb 
     +                 - (-sumc**2 + sumd*sumc)/temp
          
              endif
              enddo
              
              sum1 = a3; sum2 = b3
              sum3 = a4; sum4 = b4
              do ll = 2, nl
                  do mm = 1, nm

                      sum1 = sum1 + 1.d0/2.d0
                      sum2 = sum2 + alpham(mm,ll)**2
     +                              *dexp(-betam(mm,ll))/2.d0
                      
                      sum3 = sum3 + 1.d0/2.d0
                      sum4 = sum4 + betam(mm,ll)**2
     +                              /(2.d0*weight(ll))
                      
                  enddo
              enddo

              temp = alpham(m,l)**2*dexp(-betam(m,l))/2.d0
              suma = -sum1*temp*(sum2 - temp)/sum2**2
          
              temp = betam(m,l)**2/weight(l)
              sumb = -sum3/weight(l)*(sum4 - temp)/sum4**2
                            
              asigma = summ + suma + sumb
                            
              asigma = -1.d0/asigma*1.5d0
              if (asigma .lt. 0.0d0) asigma = -asigma      
                 
              bpdf = -fbetam(bold) 
     +             + (bold - amean)**2/(2.d0*asigma)
        
              do ii = 1, 20
	  
                  call rnset(iseed)
                  rv = drnnof()
                  call rnget(iseed)
                  anew = amean + rv*dsqrt(asigma)

                  apdf = -fbetam(anew) 
     +                 + (anew - amean)**2/(2.d0*asigma)
                  ratio = apdf - bpdf 
     
                  if (ratio .ge. 0.0d0) then
                      bpdf = apdf
                      bold = anew
                  else
                      call rnset(iseed)
                      u = drnunf()
                      call rnget(iseed)
                      if (dlog(u) .le. ratio) then
                          bpdf = apdf
                          bold = anew
                      endif
                  endif
		
              enddo     

              betam(m,l) = bold
              
          enddo
                
      enddo
            
      end subroutine
      
       
      real*8 function fbetam(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
      
      real*8 a3, b3, a4, b4
      
      integer ldum, mdum
      real*8 weight(nl)    

      real*8 star, ystar, api, sumh
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp
      real*8 suma, sumb, summ
      real*8 sum1, sum2, sum3, sum4, pdf      
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
                  
      common /va3/a3
      common /vb3/b3
      common /va4/a4
      common /vb4/b4
      
      common /vldum/ldum
      common /vmdum/mdum
      common /vweight/weight

      external dconst, dblinf

      api = dconst('pi')
      
      l = ldum
      m = mdum
      
      betam(m,l) = star
      
      summ = 0.d0
      do i = 1, ni   
      
          ll = nint(gi(i))

          if (ll .eq. l) then
          
              suma = 0.d0; sumb = 0.d0
              do j = 1, ij(i)
                    
                  sumh = 0.d0              
                  do mm = 1, nm
                          
                      temp1 = alpham(mm,ll) + bi(i,2,ll)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(mm)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(mm,ll) + bi(i,1,ll)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
                                                        
                  ystar = yij(i,j) - mu - sumh
              
                  suma = suma + ystar**2/sigmai(i)

                  sumb = sumb + ystar/sigmai(i)
          
              enddo
                                
              temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
      
              summ = summ - suma/2.d0 + sumb**2/(2.d0*temp)
          
          endif
              
      enddo
              
      sum1 = a3; sum2 = b3
      sum3 = a4; sum4 = b4
      do ll = 2, nl
          do mm = 1, nm
              
              sum1 = sum1 + 1.d0/2.d0      
              sum2 = sum2 + alpham(mm,ll)**2
     +                      *dexp(-betam(mm,ll))/2.d0
              
              sum3 = sum3 + 1.d0/2.d0
              sum4 = sum4 + betam(mm,ll)**2
     +                      /(2.d0*weight(ll))
              
          enddo
      enddo

      pdf = summ 
     +    - sum1*dlog(sum2)
     +    - sum3*dlog(sum4)
     +    - betam(m,l)/2.d0 
            
      fbetam = -pdf
    
      end function

            
	subroutine gibbs_alpham(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
      
      real*8 a3, b3
            
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      integer ldum, mdum
      real*8 ystar, api, sumh
      real*8 temp1, temp2, temp3, temp
      real*8 temp4, temp5, temp6, temp7
      real*8 dera, derb, derc
      real*8 suma, sumb, sumc
      real*8 sumd, sume, summ
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
                  
      common /va3/a3
      common /vb3/b3
      
      common /vldum/ldum
      common /vmdum/mdum
      
      external falpham, drnnof, drnunf, dconst
        
      api = dconst('PI')
            
      do m = 1, nm
          alpham(m,1) = 0.d0
      enddo                      
      do l = 2, nl
          
          ldum = l
                
          do m = 1, nm

              mdum = m
          
              bold = alpham(m,l)

              nopt = 1
              reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
              step(1) = 0.2d0 ; start(1) = bold
              call nelmin(falpham, nopt, start, xmin, ynewlo,  
     +                    reqmin, step, konvge, kcount, 
     +                    icount, numres, ifault)
              amean = xmin(1)          
              alpham(m,l) = amean
                    
              summ = 0.d0
              do i = 1, ni   
      
                  ll = nint(gi(i))
                  if (ll .eq. l) then
          
                  suma = 0.d0; sumb = 0.d0
                  sumc = 0.d0; sumd = 0.d0
                  sume = 0.d0
                  do j = 1, ij(i)
                    
                      sumh = 0.d0              
                      do mm = 1, nm
                          
                          temp1 = alpham(mm,ll) + bi(i,2,ll)
                  
                          temp2 = dexp(temp1)
     +                            /(1.d0 + dexp(temp1))
                  
                          temp3 = 2.d0*api*dfloat(mm)
     +                            *(tij(i,j) + temp2)
                  
                          temp4 = betam(mm,ll) + bi(i,1,ll)
                      
                          temp5 = dexp(temp4)
                  
                          sumh = sumh + temp5*dcos(temp3)
                      
                      enddo
                                                        
                      ystar = yij(i,j) - mu - sumh
              
                      temp1 = alpham(m,l) + bi(i,2,l) 
                      temp2 = dexp(temp1)/(1.d0 + dexp(temp1))  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                      temp4 = betam(m,l) + bi(i,1,l)  
                      temp5 = dexp(temp4)
                      temp6 = dexp(temp1)/(1.d0 + dexp(temp1))**2
                      temp7 = dexp(temp1)*(1.d0 - dexp(temp1))
     +                        /(1.d0 + dexp(temp1))**3
                                     
                      dera = -temp5*dsin(temp3)
     +                        *2.d0*api*dfloat(m)*temp6
                      
                      derb = -temp5*dcos(temp3)
     +                        *(2.d0*api*dfloat(m)*temp6)**2
                                            
                      derc = -temp5*dsin(temp3)
     +                        *(2.d0*api*dfloat(m))*temp7
                      
                      suma = suma + dera**2/sigmai(i)

                      sumb = sumb + ystar*(derb + derc)/sigmai(i)

                      sumc = sumc + dera/sigmai(i)

                      sumd = sumd + ystar/sigmai(i)

                      sume = sume + (derb + derc)/sigmai(i)
                      
                  enddo
                                
                  temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
      
                  summ = summ 
     +                 - suma + sumb 
     +                 - (-sumc**2 + sumd*sume)/temp
          
              endif
              enddo
              
              suma = a3; sumb = b3
              do ll = 2, nl
                  do mm = 1, nm

                      suma = suma + 1.d0/2.d0
                      sumb = sumb + alpham(mm,ll)**2
     +                              *dexp(-betam(mm,ll))/2.d0
                                     
                  enddo
              enddo
                        
              temp = alpham(m,l)**2*dexp(-betam(m,l))
          
              sumc = -suma*dexp(-betam(m,l))
     +                *(sumb - temp)/sumb**2
                        
              asigma = summ + sumc

              asigma = -1.d0/asigma*1.5d0
              if (asigma .lt. 0.0d0) asigma = -asigma      
                 
              bpdf = -falpham(bold) 
     +             + (bold - amean)**2/(2.d0*asigma)
        
              do ii = 1, 20
	  
                  call rnset(iseed)
                  rv = drnnof()
                  call rnget(iseed)
                  anew = amean + rv*dsqrt(asigma)

                  apdf = -falpham(anew) 
     +                 + (anew - amean)**2/(2.d0*asigma)
                  ratio = apdf - bpdf 
     
                  if (ratio .ge. 0.0d0) then
                      bpdf = apdf
                      bold = anew
                  else
                      call rnset(iseed)
                      u = drnunf()
                      call rnget(iseed)
                      if (dlog(u) .le. ratio) then
                          bpdf = apdf
                          bold = anew
                      endif
                  endif
		
              enddo     

              alpham(m,l) = bold
              
          enddo
                
      enddo
            
      end subroutine
      
       
      real*8 function falpham(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
      
      real*8 a3, b3
      
      integer ldum, mdum

      real*8 star, ystar, api, sumh
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp
      real*8 suma, sumb, summ
      real*8 sum1, sum2, pdf      
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
                  
      common /va3/a3
      common /vb3/b3
      
      common /vldum/ldum
      common /vmdum/mdum

      external dconst, dblinf

      api = dconst('pi')
      
      l = ldum
      m = mdum
      
      alpham(m,l) = star
      
      summ = 0.d0
      do i = 1, ni   
      
          ll = nint(gi(i))

          if (ll .eq. l) then
          
              suma = 0.d0; sumb = 0.d0
              do j = 1, ij(i)
                    
                  sumh = 0.d0              
                  do mm = 1, nm
                          
                      temp1 = alpham(mm,ll) + bi(i,2,ll)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(mm)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(mm,ll) + bi(i,1,ll)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
                                                        
                  ystar = yij(i,j) - mu - sumh
              
                  suma = suma + ystar**2/sigmai(i)

                  sumb = sumb + ystar/sigmai(i)
          
              enddo
                                
              temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
      
              summ = summ - suma/2.d0 + sumb**2/(2.d0*temp)
          
          endif
              
      enddo
              
      sum1 = a3; sum2 = b3
      do ll = 2, nl
          do mm = 1, nm
              
              sum1 = sum1 + 1.d0/2.d0
              sum2 = sum2 + alpham(mm,ll)**2
     +                      *dexp(-betam(mm,ll))/2.d0

          enddo
      enddo

      pdf = summ - sum1*dlog(sum2)
            
      falpham = -pdf
    
      end function

      
	subroutine gibbs_sigmab(iseed)    
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
      real*8 betam(nm,nl), theta(nl), sigmab
      real*8 a4, b4
      
      real*8 weight(nl), sumw, shape, scale, rv
      
      common /vbetam/betam
      
      common /vtheta/theta

      common /vsigmab/sigmab
      
      common /va4/a4
      common /vb4/b4
      
      external drngam
      
      sumw = 0.d0
      do l = 1, nl
          sumw = sumw + dexp(theta(l))
      enddo
      do l = 1, nl              
          weight(l) = dexp(theta(l))/sumw
      enddo
      
      shape = a4; scale = b4
      do l = 2, nl
          do m = 1, nm

              shape = shape + 1.d0/2.d0
              scale = scale + betam(m,l)**2
     +                        /(2.d0*weight(l))
              
          enddo
      enddo
      
      call rnset(iseed)
      call drngam(1,shape,rv)
      call rnget(iseed)

      sigmab = scale/rv
      
      end subroutine  

      
	subroutine gibbs_sigmaa(iseed)    
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
      real*8 betam(nm,nl), alpham(nm,nl), sigmaa
      real*8 a3, b3
      
      real*8 shape, scale, rv
      
      common /vbetam/betam
      common /valpham/alpham

      common /vsigmaa/sigmaa
      
      common /va3/a3
      common /vb3/b3
      
      external drngam
      
      shape = a3; scale = b3
      do l = 2, nl
          do m = 1, nm

              shape = shape + 1.d0/2.d0
              scale = scale + alpham(m,l)**2
     +                        *dexp(-betam(m,l))/2.d0
                       
          enddo
      enddo
      
      call rnset(iseed)
      call drngam(1,shape,rv)
      call rnget(iseed)

      sigmaa = scale/rv
      
      end subroutine  
      

      subroutine gibbs_bi2(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl), omega(nb,nb,nl)
      real*8 gi(ni)
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      integer ldum, idum
      real*8 omegal(nb,nb), omegai(nb,nb)
      real*8 ystar, api, sumh
      real*8 sumh1, sumh2, sumh3
      real*8 temp1, temp2, temp3, temp4
      real*8 temp5, temp6, temp7, temp
      real*8 suma, sumb, sumc, sumd, sume
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      common /vomega/omega
      
      common /vgi/gi
            
      common /vldum/ldum
      common /vidum/idum
      common /vomegai/omegai
      
      external fbi2, drnnof, drnunf, dconst, dlinrg
             
      api = dconst('PI')
      
      do i = 1, ni
          bi(i,2,1) = 0.d0     
      enddo            
      do l = 2, nl
          
          ldum = l

          do j1 = 1, nb
              do j2 = 1, nb
                  omegal(j1,j2) = omega(j1,j2,l)
              enddo
          enddo

          call dlinrg(nb, omegal, nb, omegai, nb)
                              
          do i = 1, ni

              idum = i
                  
              bold = bi(i,2,l)

              nopt = 1
              reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
              step(1) = 0.2d0 ; start(1) = bold
              call nelmin(fbi2, nopt, start, xmin, ynewlo,  
     +                    reqmin, step, konvge, kcount, 
     +                    icount, numres, ifault)
              amean = xmin(1)                          
              bi(i,2,l) = amean
              
              suma = 0.d0; sumb = 0.d0
              sumc = 0.d0; sumd = 0.d0
              sume = 0.d0
              do j = 1, ij(i)
                    
                  sumh = 0.d0 ; sumh1 = 0.d0
                  sumh2 = 0.d0; sumh3 = 0.d0
                  do m = 1, nm
                          
                      temp1 = alpham(m,l) + bi(i,2,l)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(m,l) + bi(i,1,l)
                      
                      temp5 = dexp(temp4)
                  
                      temp6 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))**2
                      
                      temp7 = dexp(temp1)
     +                        *(1.d0 - dexp(temp1))
     +                        /(1.d0 + dexp(temp1))**3
                      
                      sumh = sumh + temp5*dcos(temp3)
                      
                      sumh1 = sumh1 
     +                      - temp5*dsin(temp3)
     +                        *(2.d0*api*dfloat(m))*temp6

                      sumh2 = sumh2 
     +                      - temp5*dcos(temp3)
     +                        *(2.d0*api*dfloat(m)*temp6)**2

                      sumh3 = sumh3 
     +                      - temp5*dsin(temp3)
     +                        *(2.d0*api*dfloat(m))*temp7
                      
                  enddo
                                                        
                  ystar = yij(i,j) - mu - sumh
              
                  suma = suma + sumh1**2/sigmai(i)

                  sumb = sumb + ystar*(sumh2 + sumh3)/sigmai(i)

                  sumc = sumc + sumh1/sigmai(i)
                  
                  sumd = sumd + ystar/sigmai(i)

                  sume = sume + (sumh2 + sumh3)/sigmai(i)
                  
              enddo

              temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
              
              asigma = -suma + sumb 
     +               - (-sumc**2 + sumd*sume)/temp
     +               - omegai(2,2)
                     
              asigma = -1.d0/asigma*1.5d0
              if (asigma .lt. 0.0d0) asigma = -asigma      
                 
              bpdf = -fbi2(bold) 
     +             + (bold - amean)**2/(2.d0*asigma)
        
              do ii = 1, 20
	  
                  call rnset(iseed)
                  rv = drnnof()
                  call rnget(iseed)
                  anew = amean + rv*dsqrt(asigma)

                  apdf = -fbi2(anew) 
     +                 + (anew - amean)**2/(2.d0*asigma)
                  ratio = apdf - bpdf 
     
                  if (ratio .ge. 0.0d0) then
                      bpdf = apdf
                      bold = anew
                  else
                      call rnset(iseed)
                      u = drnunf()
                      call rnget(iseed)
                      if (dlog(u) .le. ratio) then
                          bpdf = apdf
                          bold = anew
                      endif
                  endif
		
              enddo     

              bi(i,2,l) = bold
                                       
          enddo
          
      enddo
            
      end subroutine
      
       
      real*8 function fbi2(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
            
      integer ldum, idum
      real*8 omegai(nb,nb)

      real*8 star, ystar, api 
      real*8 bstar(nb), sumh
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp
      real*8 suma, sumb, sumc, pdf      
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
                  
      common /vldum/ldum
      common /vidum/idum
      common /vomegai/omegai

      external dconst, dblinf

      api = dconst('PI')
      
      l = ldum
      i = idum
      
      bi(i,2,l) = star
      
      suma = 0.d0; sumb = 0.d0
      do j = 1, ij(i)
                    
          sumh = 0.d0              
          do m = 1, nm
                          
              temp1 = alpham(m,l) + bi(i,2,l)
                  
              temp2 = dexp(temp1)
     +                /(1.d0 + dexp(temp1))
                  
              temp3 = 2.d0*api*dfloat(m)
     +                *(tij(i,j) + temp2)
                  
              temp4 = betam(m,l) + bi(i,1,l)
                      
              temp5 = dexp(temp4)
                  
              sumh = sumh + temp5*dcos(temp3)
                      
          enddo
                                                        
          ystar = yij(i,j) - mu - sumh
              
          suma = suma + ystar**2/sigmai(i)

          sumb = sumb + ystar/sigmai(i)
          
      enddo
                        
      do jj = 1, nb          
          bstar(jj) = bi(i,jj,l)          
      enddo
      sumc = dblinf(nb,nb,omegai,nb,bstar,bstar)
        
      temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
      
      pdf = -suma/2.d0 + sumb**2/(2.d0*temp) - sumc/2.d0
            
      fbi2 = -pdf
    
      end function
      
      
      subroutine gibbs_bi1(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl), omega(nb,nb,nl)
      real*8 gi(ni)
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      integer ldum, idum
      real*8 omegal(nb,nb), omegai(nb,nb)
      real*8 ystar, api, sumh
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp
      real*8 suma, sumb, sumc, sumd
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      common /vomega/omega
      
      common /vgi/gi
            
      common /vldum/ldum
      common /vidum/idum
      common /vomegai/omegai
      
      external fbi1, drnnof, drnunf, dconst, dlinrg
             
      api = dconst('PI')
      
      do i = 1, ni
          bi(i,1,1) = 0.d0     
      enddo      
      do l = 2, nl
          
          ldum = l

          do j1 = 1, nb
              do j2 = 1, nb
                  omegal(j1,j2) = omega(j1,j2,l)
              enddo
          enddo

          call dlinrg(nb, omegal, nb, omegai, nb)
                              
          do i = 1, ni

              idum = i
                  
              bold = bi(i,1,l)

              nopt = 1
              reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
              step(1) = 0.2d0 ; start(1) = bold
              call nelmin(fbi1, nopt, start, xmin, ynewlo,  
     +                    reqmin, step, konvge, kcount, 
     +                    icount, numres, ifault)
              amean = xmin(1)                          
              bi(i,1,l) = amean
                                               
              suma = 0.d0; sumb = 0.d0
              sumc = 0.d0; sumd = 0.d0
              do j = 1, ij(i)
                    
                  sumh = 0.d0              
                  do m = 1, nm
                          
                      temp1 = alpham(m,l) + bi(i,2,l)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(m,l) + bi(i,1,l)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
                                                        
                  ystar = yij(i,j) - mu - sumh
              
                  suma = suma + sumh**2/sigmai(i)

                  sumb = sumb + ystar*sumh/sigmai(i)
                  
                  sumc = sumc + sumh/sigmai(i)

                  sumd = sumd + ystar/sigmai(i)
                  
              enddo

              temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
              
              asigma = -suma + sumb 
     +               - (-sumc**2 + sumd*sumc)/temp
     +               - omegai(1,1) 
              
              asigma = -1.d0/asigma*1.5d0
              if (asigma .lt. 0.0d0) asigma = -asigma      
                 
              bpdf = -fbi1(bold) 
     +             + (bold - amean)**2/(2.d0*asigma)
        
              do ii = 1, 20
	  
                  call rnset(iseed)
                  rv = drnnof()
                  call rnget(iseed)
                  anew = amean + rv*dsqrt(asigma)

                  apdf = -fbi1(anew) 
     +                 + (anew - amean)**2/(2.d0*asigma)
                  ratio = apdf - bpdf 
     
                  if (ratio .ge. 0.0d0) then
                      bpdf = apdf
                      bold = anew
                  else
                      call rnset(iseed)
                      u = drnunf()
                      call rnget(iseed)
                      if (dlog(u) .le. ratio) then
                          bpdf = apdf
                          bold = anew
                      endif
                  endif
		
              enddo     

              bi(i,1,l) = bold
                                       
          enddo
          
      enddo
                  
      end subroutine
      
       
      real*8 function fbi1(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz, sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
            
      integer ldum, idum
      real*8 omegai(nb,nb)

      real*8 star, ystar, api
      real*8 bstar(nb), sumh
      real*8 temp1, temp2, temp3
      real*8 temp4, temp5, temp
      real*8 suma, sumb, sumc, pdf      
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
                  
      common /vldum/ldum
      common /vidum/idum
      common /vomegai/omegai

      external dconst, dblinf

      api = dconst('PI')
      
      l = ldum
      i = idum
      
      bi(i,1,l) = star
      
      suma = 0.d0; sumb = 0.d0
      do j = 1, ij(i)
                    
          sumh = 0.d0              
          do m = 1, nm
                          
              temp1 = alpham(m,l) + bi(i,2,l)
                  
              temp2 = dexp(temp1)
     +                /(1.d0 + dexp(temp1))
                  
              temp3 = 2.d0*api*dfloat(m)
     +                *(tij(i,j) + temp2)
                  
              temp4 = betam(m,l) + bi(i,1,l)
                      
              temp5 = dexp(temp4)
                  
              sumh = sumh + temp5*dcos(temp3)
                      
          enddo
                                                        
          ystar = yij(i,j) - mu - sumh
              
          suma = suma + ystar**2/sigmai(i)

          sumb = sumb + ystar/sigmai(i)
          
      enddo
                        
      do jj = 1, nb          
          bstar(jj) = bi(i,jj,l)          
      enddo
      sumc = dblinf(nb,nb,omegai,nb,bstar,bstar)
        
      temp = 1.d0/sigmaz + dfloat(ij(i))/sigmai(i)
      
      pdf = -suma/2.d0 + sumb**2/(2.d0*temp) - sumc/2.d0
            
      fbi1 = -pdf
    
      end function
      
      
	subroutine gibbs_mu(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz
      real*8 sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
                              
      real*8 sigma_mu
      
      real*8 ystar, api, sumh, sumj
      real*8 temp1, temp2, temp3, temp4, temp5
            
      real*8 hmean, hsigma, amean, asigma, rv
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
            
      common /vsigma_mu/sigma_mu
      
      external dconst, drnnof
          
      api = dconst('pi')

      hmean = 0.d0      
      hsigma = 1.d0/sigma_mu
      
      do i = 1, ni
          
          l = nint(gi(i))
          
          sumj = 0.d0
          do j = 1, ij(i)
          
              sumh = 0.d0              
              if (l .ne. 1) then
                  do m = 1, nm
                          
                      temp1 = alpham(m,l) + bi(i,2,l)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(m,l) + bi(i,1,l)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
              endif
              
              ystar = yij(i,j) - sumh
              
              sumj = sumj + ystar/sigmai(i)
                                          
          enddo
          
          temp1 = dfloat(ij(i))/sigmai(i)          
          temp2 = 1.d0/sigmaz + temp1
          
          hmean = hmean + sumj*(1.d0 - temp1/temp2)
          hsigma = hsigma + temp1*(1.d0 - temp1/temp2)
          
      enddo     
          
      amean = hmean/hsigma
      asigma = 1.d0/hsigma
                            
      call rnset(iseed)
      rv = drnnof()
      call rnget(iseed)
		    
      mu = amean + rv*dsqrt(asigma)   
          
      end subroutine	

      
	subroutine gibbs_sigmaz(iseed)   
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, sigmaz
      real*8 sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
                     
      real*8 a2, b2
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      real*8 ystar, api, sumh, sumi(ni)
      real*8 temp1, temp2, temp3, temp4, temp5
      real*8 star, temp, summ
            
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
            
      common /va2/a2
      common /vb2/b2

      common /vsumi/sumi
      
      external fsigmaz, drnnof, drnunf
                
      api = dconst('pi')
      
      do i = 1, ni
          
          l = nint(gi(i))
          
          sumi(i) = 0.d0          
          do j = 1, ij(i)
          
              sumh = 0.d0              
              if (l .ne. 1) then
                  do m = 1, nm
                          
                      temp1 = alpham(m,l) + bi(i,2,l)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(m,l) + bi(i,1,l)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
              endif
              
              ystar = yij(i,j) - mu - sumh
              
              sumi(i) = sumi(i) + ystar/sigmai(i)
                                          
          enddo
                              
      enddo     
                    
      bold = dlog(sigmaz)

      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fsigmaz, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
	amean = xmin(1)
      star = dexp(-amean)
      
      summ = 0.d0
      do i = 1, ni
          
          temp = star + dfloat(ij(i))/sigmai(i) 
          
          summ = summ 
     +         - 1.d0/2.d0*star*(dfloat(ij(i))/sigmai(i))
     +           /temp**2 
     +         - sumi(i)**2/2.d0
     +           *star*(dfloat(ij(i))/sigmai(i) - star)
     +           /temp**3
                                        
      enddo     
                    
      asigma = summ - b2*star
      
      asigma = -1.d0/asigma*1.5d0
      if (asigma .lt. 0.0d0) asigma = -asigma      
                 
      bpdf = -fsigmaz(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
        
      do ii = 1, 20
	  
          call rnset(iseed)
          rv = drnnof()
          call rnget(iseed)
          anew = amean + rv*dsqrt(asigma)

          apdf = -fsigmaz(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
              bpdf = apdf
              bold = anew
          else
              call rnset(iseed)
              u = drnunf()
              call rnget(iseed)
              if (dlog(u) .le. ratio) then
                  bpdf = apdf
                  bold = anew
              endif
          endif
		
      enddo     

      sigmaz = dexp(bold)      
      
      end subroutine	
      
      
      real*8 function fsigmaz(star)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
      real*8 sigmai(ni)                     
      real*8 a2, b2
            
      real*8 sumi(ni)
      
      real*8 star, temp, summ, pdf
            
      common /vij/ij
            
      common /vsigmai/sigmai
                  
      common /va2/a2
      common /vb2/b2

      common /vsumi/sumi
            
      summ = 0.d0
      do i = 1, ni
          
          temp = dexp(-star) + dfloat(ij(i))/sigmai(i) 
          
          summ = summ 
     +         - star/2.d0 
     +         - dlog(temp)/2.d0 
     +         + sumi(i)**2/(2.d0*temp)
                                        
      enddo     
                    
      pdf = summ - a2*star - b2*dexp(-star)

      fsigmaz = -pdf
    
	end function
      
      
	subroutine gibbs_zeta(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, zeta(ni), sigmaz
      real*8 sigmai(ni)
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)
                              
      real*8 ystar, api, sumh
      real*8 temp1, temp2, temp3, temp4, temp5
      
      real*8 hmean, hsigma, amean, asigma, rv
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vzeta/zeta
      common /vsigmaz/sigmaz
      
      common /vsigmai/sigmai
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
            
      external dconst, drnnof
          
      api = dconst('pi')
      
      do i = 1, ni

          l = nint(gi(i))
          
          hmean = 0.d0                
          hsigma = 1.d0/sigmaz
          
          do j = 1, ij(i)
          
              sumh = 0.d0              
              if (l .ne. 1) then
                  do m = 1, nm
                          
                      temp1 = alpham(m,l) + bi(i,2,l)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(m,l) + bi(i,1,l)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
              endif
              
              ystar = yij(i,j) - mu - sumh
              
              hmean = hmean + ystar/sigmai(i)
              hsigma = hsigma + 1.d0/sigmai(i)
                            
          enddo
          
          amean = hmean/hsigma
          asigma = 1.d0/hsigma
                            
          call rnset(iseed)
          rv = drnnof()
          call rnget(iseed)
		    
          zeta(i) = amean + rv*dsqrt(asigma)   
                    
      enddo     
          
      end subroutine	
            
      
	subroutine gibbs_sigmay(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, zeta(ni), sigmay
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)

      real*8 a0, b0, a1, b1

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      real*8 ystar, api, sumh
      real*8 temp1, temp2, temp3, temp4, temp5
      real*8 sumi(ni), summ
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vzeta/zeta
      
      common /vsigmay/sigmay
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
      
      common /va0/a0
      common /vb0/b0
      common /va1/a1
      common /vb1/b1
      
      common /vsumi/sumi
      
      external fsigmay, drnnof, drnunf, dconst

      api = dconst('pi')
          
      do i = 1, ni
                               
          l = nint(gi(i))

          sumi(i) = 0.d0
          do j = 1, ij(i)
                               
              sumh = 0.d0              
              if (l .ne. 1) then
                  do m = 1, nm
                          
                      temp1 = alpham(m,l) + bi(i,2,l)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(m,l) + bi(i,1,l)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
              endif
              
              ystar = yij(i,j) - mu - zeta(i) - sumh

              sumi(i) = sumi(i) + ystar**2/2.d0
                                
          enddo
          
      enddo
      
      bold = dlog(sigmay)

      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcouns = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fsigmay, nopt, start, xmin, ynewlo,  
     +            reqmin, step, konvge, kcount, 
     +            icount, numres, ifault)
	amean = xmin(1)
      sigmay = dexp(amean)
      
      summ = 0.d0
      do i = 1, ni

          summ = summ 
     +         - (a0 + dfloat(ij(i))/2.d0)
     +           *b0*sigmay*sumi(i)/(b0*sigmay + sumi(i))**2
                              
      enddo         
      asigma = summ - b1*dexp(amean)
      
      asigma = -1.d0/asigma*1.5d0
      if (asigma .lt. 0.0d0) asigma = -asigma      
                 
      bpdf = -fsigmay(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)
        
      do ii = 1, 20
	  
          call rnset(iseed)
          rv = drnnof()
          call rnget(iseed)
          anew = amean + rv*dsqrt(asigma)

          apdf = -fsigmay(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
              bpdf = apdf
              bold = anew
          else
              call rnset(iseed)
              u = drnunf()
              call rnget(iseed)
              if (dlog(u) .le. ratio) then
                  bpdf = apdf
                  bold = anew
              endif
          endif
		
      enddo     

      sigmay = dexp(bold)
            
      end subroutine	
  
      
      real*8 function fsigmay(star)      
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
      real*8 a0, b0, a1, b1
      real*8 sumi(ni)
      
      real*8 star, summ, pdf

      common /vij/ij
      
      common /va0/a0
      common /vb0/b0
      common /va1/a1
      common /vb1/b1

      common /vsumi/sumi
                
      sigmay = dexp(star)
      
      summ = 0.d0
      do i = 1, ni

          summ = summ + a0*star 
     +                - (a0 + dfloat(ij(i))/2.d0)
     +                   *dlog(b0*sigmay + sumi(i))
                              
      enddo   
      
      pdf = summ + a1*star - b1*dexp(star)
                    
      fsigmay = -pdf
    
      end function    

      
	subroutine gibbs_sigmai(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 
		
	integer ij(ni)
	real*8 yij(ni,nj), tij(ni,nj)
      real*8 mu, zeta(ni)
      real*8 sigmai(ni), sigmay
      real*8 betam(nm,nl), alpham(nm,nl) 
      real*8 bi(ni,nb,nl)
      real*8 gi(ni)

      real*8 a0, b0

      real*8 shape, scale, rv
      real*8 ystar, api, sumh
      real*8 temp1, temp2, temp3, temp4, temp5
      
      common /vij/ij
      
      common /vyij/yij
      common /vtij/tij
            
      common /vmu/mu
      common /vzeta/zeta
      
      common /vsigmai/sigmai
      common /vsigmay/sigmay
      
      common /vbetam/betam
      common /valpham/alpham

      common /vbi/bi
      
      common /vgi/gi
      
      common /va0/a0
      common /vb0/b0
      
      external dconst, drngam

      api = dconst('pi')
      
      do i = 1, ni
                
          l = nint(gi(i))

          shape = a0; scale = b0*sigmay          
          do j = 1, ij(i)
                               
              sumh = 0.d0              
              if (l .ne. 1) then
                  do m = 1, nm
                          
                      temp1 = alpham(m,l) + bi(i,2,l)
                  
                      temp2 = dexp(temp1)
     +                        /(1.d0 + dexp(temp1))
                  
                      temp3 = 2.d0*api*dfloat(m)
     +                        *(tij(i,j) + temp2)
                  
                      temp4 = betam(m,l) + bi(i,1,l)
                      
                      temp5 = dexp(temp4)
                  
                      sumh = sumh + temp5*dcos(temp3)
                      
                  enddo
              endif
              
              ystar = yij(i,j) - mu - zeta(i) - sumh

              shape = shape + 1.d0/2.d0
              scale = scale + ystar**2/2.d0
                                
          enddo
        
          call rnset(iseed)
          call drngam(1,shape,rv)
          call rnget(iseed)

          sigmai(i) = scale/rv
                    
      enddo     
            
      end subroutine	
            
      
	subroutine gibbs_omega(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 

      real*8 bi(ni,nb,nl), omega(nb,nb,nl), gi(ni)
      real*8 d0, s0
      
	real*8 shape, qq(nb,nb), omegal(nb,nb)
      
      common /vbi/bi
      common /vomega/omega

      common /vgi/gi
      
      common /vd0/d0
      common /vs0/s0
          
      do l = 2, nl
          
	    shape = d0	      
          do j1 = 1, nb
              do j2 = 1, nb
                  qq(j1,j2) = 0.d0
              enddo
              qq(j1,j1) = s0
          enddo

          do i = 1, ni            

              ll = nint(gi(i))
              
              if (ll .eq. l) then
          
	            shape = shape + 1.d0
                               
                  do j1 = 1, nb
                      do j2 = 1, nb
                          qq(j1,j2) = qq(j1,j2) 
     +                              + bi(i,j1,ll)
     +                                *bi(i,j2,ll)
                      enddo
                  enddo
                  
              endif
          
          enddo
      
	    call gen_iwish(omegal,shape,qq,nb,iseed)
          
          do j1 = 1, nb
              do j2 = 1, nb
                  omega(j1,j2,l) = omegal(j1,j2)
              enddo
          enddo
          
      enddo
      
	end subroutine
      
      
	subroutine gen_iwish(sigma,shape,qq,nn,iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 

      real*8 sigma(nn,nn), qq(nn,nn), shape
      
	real*8 qqinv(nn,nn)
	real*8 df, chi(nn), anv(nn,nn), bb(nn,nn)
	real*8 tol, rsig(nn,nn), suma
	real*8 btemp(nn,nn), aainv(nn,nn)
            
	external drnchi, drnnof, dlinrg, 
     +         dchfac, dmach, dmxytf, dmrrrr 
                  
      do j = 1, nn
      
		df = shape - dfloat(j) + 1.d0
		
		call rnset(iseed)
		call drnchi(1, df, chi(j))
		call rnget(iseed)
		
		do i = 2, nn
			call rnset(iseed)		
			anv(j,i) = drnnof()
			call rnget(iseed)
		enddo
		
	enddo

      bb(1,1) = chi(1)
      do j = 2, nn

		suma = 0.0d0
		do i = 1, j-1         
			suma = suma + anv(i,j)**2
		enddo
		bb(j,j) = chi(j) + suma

		bb(1,j) = anv(1,j)*dsqrt(chi(1))
          bb(j,1) = bb(1,j)
        
          if (j .gt. 2) then
		    do i = 2, j-1         
			    suma = 0.0d0
			    do k = 1, i-1
				    suma = suma + anv(k,j)*anv(k,i)
			    enddo
			    bb(i,j) = anv(i,j)*dsqrt(chi(i)) + suma
		        bb(j,i) = bb(i,j)
		    enddo
	    endif
	  
	enddo
      
	call dlinrg(nn, qq, nn, qqinv, nn)	

	tol = 100.0d0*dmach(4)	
	call dchfac(nn, qqinv, nn, tol, irank, rsig, nn)
	
      call dmxtyf(nn,nn,rsig,nn,nn,nn,bb,
     +            nn,nn,nn,btemp,nn)
      call dmrrrr(nn,nn,btemp,nn,nn,nn,rsig,
     +            nn,nn,nn,aainv,nn)
		
	call dlinrg(nn, aainv, nn, sigma, nn)	      
                                                        
      end subroutine      


      subroutine gen_wish(sigma,shape,qq,nn,iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 663, nj = 13
	integer, parameter :: nl = 3, nb = 2, nm = 3 

      real*8 sigma(nn,nn), qq(nn,nn), shape
      
	real*8 qqinv(nn,nn)
	real*8 df, chi(nn), anv(nn,nn), bb(nn,nn)
	real*8 tol, rsig(nn,nn), suma
	real*8 btemp(nn,nn), aainv(nn,nn)
            
	external drnchi, drnnof, dlinrg, 
     +         dchfac, dmach, dmxytf, dmrrrr 
                  
      do j = 1, nn
      
          df = shape - dfloat(j) + 1.d0
		
		call rnset(iseed)
		call drnchi(1, df, chi(j))
		call rnget(iseed)
		
		do i = 2, nn
			call rnset(iseed)		
			anv(j,i) = drnnof()
			call rnget(iseed)
		enddo
		
	enddo

      bb(1,1) = chi(1)
      do j = 2, nn

		suma = 0.0d0
		do i = 1, j-1         
			suma = suma + anv(i,j)**2
		enddo
		bb(j,j) = chi(j) + suma

		bb(1,j) = anv(1,j)*dsqrt(chi(1))
          bb(j,1) = bb(1,j)
        
          if (j .gt. 2) then
              do i = 2, j-1         
			    suma = 0.0d0
			    do k = 1, i-1
				    suma = suma + anv(k,j)*anv(k,i)
			    enddo
			    bb(i,j) = anv(i,j)*dsqrt(chi(i)) + suma
		        bb(j,i) = bb(i,j)
		    enddo
	    endif
	  
	enddo
      
	call dlinrg(nn, qq, nn, qqinv, nn)	

	tol = 100.0d0*dmach(4)	
	call dchfac(nn, qqinv, nn, tol, irank, rsig, nn)
	
      call dmxtyf(nn,nn,rsig,nn,nn,nn,bb,
     +            nn,nn,nn,btemp,nn)
      call dmrrrr(nn,nn,btemp,nn,nn,nn,rsig,
     +            nn,nn,nn,sigma,nn)
		                                                        
      end subroutine