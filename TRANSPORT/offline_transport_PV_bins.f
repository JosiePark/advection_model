c PROGRAM THAT PERFORMS OFFLINE ADVECTION
c RELEASES PARTICLES UNIFORMLY AT DIFFERENT TIMES

      program offline_transport
      
      use mod_advection_input
      use mod_PVbin_advection_constants
      use mod_variables
      use mod_qg2_netcdf
      use mod_traj_netcdf
      use mod_2dcubic
      use mod_bicubic
      use mod_random
      use mod_time_interp
      use mod_1dinterp
      use mod_laplace
      use mod_rk4
      
      implicit none 
      
c ---- DEFINE DIMENSIONAL VARIABLES
      
      uscale=1.
      scale=basinscale/dfloat(ii) 
      tscale=scale/uscale
      SS=(scale/Rd)**2

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      beta_nondim=beta*scale*scale/uscale
      
      BETA_NONDIM_U1=BETA_NONDIM+U_0*S1
      BETA_NONDIM_U2=BETA_NONDIM-U_0*S2
        
        do j = 1,jj1
            y_c(j) = dfloat(j-1)
        enddo

        
 
c --- READ TIME AND TIME AVERAGED STREAM FUNCTION ----
      
      call read_time(file_name,time,t_len)
      call read_ave_file(ave_name,ii,jj,psi1_av,psi2_av,time_av)


c ----- CALCULATE THE TRUE TIME ARRAY (WHICH TAKES INTO ACCOUNT THE CHOSEN OFFLINE TIME STEP)

        release_length_secs = release_length*86400 ! time recording in seconds for a release
        max_time_lag = time(t_len)*86400
        
        t_tot = int(release_length_secs/dt) ! total number of time steps per release
        l_tot_day = int(86400/dt) ! total number of time steps per day

        dt_nondim = dt/tscale
            
        allocate(time_o(t_tot))
        allocate(time_dim(t_tot))
        
c ------ CREATE TRAJECTORY DATA FILE ----------------

      if(i_pvbin.eq.1) then

      if (i_eddy .eq. 1) then
        call create_binned_file(eddy_name,nbins,npoints,release_no)
      endif
      if(i_pseudo .eq. 1) then
         call create_binned_file(pseudo_name,nbins,npoints,release_no)
      endif
      if (i_full .eq. 1) then
          call create_binned_file(full_name,nbins,npoints,release_no)
      endif
      
      else
    
            if (i_eddy .eq. 1) then
        call create_binned_file(eddy_name,nbins,npoints,release_no)
      endif
      if(i_pseudo .eq. 1) then
         call create_binned_file(pseudo_name,nbins,npoints,release_no)
      endif
      if (i_full .eq. 1) then
          call create_binned_file(full_name,nbins,npoints,release_no)
      endif
      
      endif

        
c ----- Initialise trajectory arrays
        
        if((i_full.eq.1).or.(i_pseudo.eq.1)) then
        allocate(x1r(nbins,npoints), 
     &   x2r(nbins,npoints)
     & ,y1r(nbins,npoints),y2r(nbins,npoints)
     & ,x1r_coord(nbins,npoints)
     & ,x2r_coord(nbins,npoints)
     & ,y1r_coord(nbins,npoints)
     & ,y2r_coord(nbins,npoints))
        endif
        
        if (i_eddy.eq.1) then
        allocate(x1r_eddy(nbins,npoints)
     & ,x2r_eddy(nbins,npoints)
     & ,y1r_eddy(nbins,npoints)
     & ,y2r_eddy(nbins,npoints)
     & ,x1r_eddy_coord(nbins,npoints)
     & ,x2r_eddy_coord(nbins,npoints)
     & ,y1r_eddy_coord(nbins,npoints)
     & ,y2r_eddy_coord(nbins,npoints))
        endif
        
        if (i_pseudo.eq.1) then        
        allocate(x1r_pseudo(nbins,npoints)
     & ,x2r_pseudo(nbins,npoints)
     & ,y1r_pseudo(nbins,npoints)
     & ,y2r_pseudo(nbins,npoints)
     & ,x1r_pseudo_coord(nbins,npoints)
     & ,x2r_pseudo_coord(nbins,npoints)
     & ,y1r_pseudo_coord(nbins,npoints)
     & ,y2r_pseudo_coord(nbins,npoints)
     & ,x1r_mean(nbins,npoints)
     & ,y1r_mean(nbins,npoints)
     & ,x2r_mean(nbins,npoints)
     & ,y2r_mean(nbins,npoints))
        endif    
        
        allocate(release_time(release_no) ,nrel(release_no))
        allocate(k_s(release_no))
        
        
        do k = 1,release_no    
            release_time(k) = (k-1)*release_interval
            nrel(k) = 1
            k_s(k) = 0
        enddo
        
        
      if ((i_pseudo.eq.1).or.(i_eddy.eq.1)) then
      

      
c ---- Calculate coefficients for time-averaged stream function
c ---- Needed for spatial interpolation
c ------- ONLY FOR PSEUDO TRAJECTORIES -----------


      if(i_pseudo.eq.1) then


      if (isolve.eq.0) then
      
        
      call A_matrix(ii,jj,psi1_av,M1_av)
      call A_matrix(ii,jj,psi2_av,M2_av)
    

      
      else if (isolve.eq.1) then
    
      
      call cubic_coeff_x(ii,jj,psi1_av
     & ,a1_av,b1_av,c1_av,d1_av)
      call cubic_coeff_x(ii,jj,psi2_av
     & ,a2_av,b2_av,c2_av,d2_av)
     
     
      endif
      endif
      endif
      

       
        iseed = 102
        
        nrec = 0
        do k = 1,release_no
        
        print*, 'Starting particle advection for release number',k
        print*, 't_tot =',t_tot
        
c ------------ MAIN CYCLE -----------------------------

        do t = 1,t_tot
            time_o(t) = (t-1)*dt_nondim + release_time(k)*86400/tscale
            time_dim(t) = time_o(t)*tscale/86400
            
        enddo
 
        do t = 1,t_tot
        
        if (t.eq.1) then
        
      if (i_pvbin.eq.1) then
      
      
c ------- CONSTRUCT BINS UNIFORM IN PV --------------------
      
c Calculate zonally time-averaged full PV across whole domain
c Add extra domain above and blow as the PV snapshot may fall outside of range of PV_bar
      
c      call PV_bar_from_psi(ii,jj,jj1,basinscale,beta
c     & ,Rd,H1,H2,U_0,psi1_av,psi2_av,PV_bar)

c CALCULATE FULL ZONALLY AVERAGED PV

      call rel_from_psi(ii,jj,psi1_av,psi2_av,rel1_av,rel2_av) ! calculate time-averaged relative vorticity
      print*,'relative vorticity done'
      call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1_av,psi2_av,rel1_av
     & ,rel2_av,zeta1_av,zeta2_av) ! time averaged PV anomaly
       print *,'PV anomaly done'
c CALCULATE ZONAL FULL TIME-AVERAGED PV

      do j = 1,jj1
            beta1_y(j) = beta_nondim_u1*(dfloat(j-1)) ! only do for the top layer
      enddo
      print*,'beta done'
      do j=1,jj
      do kk = 1,coord_range
            idx = j+(kk-1)*jj
            PV_bar(idx) = 0.
            do i = 1,ii
                PV_av(i,idx) = zeta1_av(i,j) + beta1_y(idx)
                PV_bar(idx) = PV_bar(idx) + PV_av(i,idx)
            enddo
            PV_bar(idx) = PV_bar(idx)/dfloat(ii)
      enddo
      enddo
     
      print*,'PV_bar calculated'
     
c READ STREAM FUNCTION

c ----- find entry m in time such that time(m) < time_o(k) <= time(m+1)  
        do m = 1,t_len
            if (time(m) >= time_dim(t)) then
            k_new = m
                if(k_new .eq. 2)then
                k_new = 3
                else if (k_new.eq.t_len)then
                k_new = k_new-1
                else if(k_new.eq.1) then
                k_new = 3
                endif
            exit
            endif
        enddo
        
        print*,'found time index'

c ----- cubic interpolate psi in time to find psi
c ----- at time_o(k), time_o(k-1)

c ---- time_o(k)
        
        call read_psi(file_name,ii,jj,k_new-2,psi1,psi2)
        
        print*,'read psi'

        
        time_cubic = time(k_new-2:k_new+1)
     
        call interp_time(ii,jj,time_cubic,time_dim(t),psi1,psi1_new)
        call interp_time(ii,jj,time_cubic,time_dim(t),psi2,psi2_new) ! psi1_new is the stream function
        
        print*,'time interpolated'
c------------------------------
C CALCULATE PV
c ----------------------------------

      call rel_from_psi(ii,jj,psi1_new,psi2_new,rel1,rel2) ! calculate time-averaged relative vorticity
      call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1_new,psi2_new,rel1
     & ,rel2,zeta1,zeta2)
     
c ---------- make sure PV matches with PV_bar in regards to beta*y
      do j = jj+1,jj*2
            beta1_y(j-jj) = beta_nondim_u1*(dfloat(j-1)) ! only do for the top layer
      enddo
      
      do j=1,jj
            do i = 1,ii
                PV(i,j) = zeta1(i,j) + beta1_y(j)
            enddo
      enddo
 

        call cubic_coeff_x(ii,jj
     &    ,PV
     &    ,a,b,c,d)      ! calculate coefficients required for interpolation   
     
     
c Bin the PV

     
      d_bin = (maxval(PV)-minval(PV))/dfloat(nbins)
      
      do j = 1,jj1
        j_bin(j) = dfloat(j-1)
      enddo
      
     
      do p = 1,nbins
    

        PV_bin(p) = minval(PV) + dfloat(p)*d_bin     

      enddo
      
      print*,'pv bins constructed:',PV_bin
      
      
c --------------------------------------------------------
C GENERATE PARTICLES
c -----------------------------
        do p = 1,nbins
            bin_count(p) = 0
            bin_stop(p) = 0
        enddo
        points_count = 0
        
10      continue
        x0 = ran1(iseed)*dfloat(ii)
        y0 = ran1(iseed)*dfloat(jj)
        
c -------------------------------------
C CALCULATE PV AT PARTICLE LOCATION
c -------------------------------------

        call cubic_poly_x(ii,jj,x0 ! construct polynomial PV_x
     &    ,y0,a,b,c,d,PV_x)
     
        print*,'cubic polynomial constructed'
     
        call cubic_interp(ii,jj,PV_x,y0,PV_0) 
        
        print*,'interpolated'
C --------------------------------------  
c SORT PV_0 INTO APPROPRIATE BIN
C ----------------------------------

        do p = 1,nbins
            if (PV_0.lt.PV_bin(p)) then
                bin_index = p
                go to 11
            endif
        enddo
        
        print*,'sorted into bin'
        
11     continue
                bin_count(bin_index) = bin_count(bin_index) + 1
                print*,'bin_count,p=',bin_count(bin_index),bin_index
                if (bin_count(bin_index).le.npoints) then
                    if ((i_full.eq.1).or.(i_pseudo.eq.1)) then
                
                    x1r(bin_index,bin_count(bin_index)) = x0
                    y1r(bin_index,bin_count(bin_index)) = y0
                    x2r(bin_index,bin_count(bin_index)) = x0
                    y2r(bin_index,bin_count(bin_index)) = y0
                    
                    x1r_coord(bin_index,bin_count(bin_index)) = 0
                    y1r_coord(bin_index,bin_count(bin_index)) = 0
                    x2r_coord(bin_index,bin_count(bin_index)) = 0
                    y2r_coord(bin_index,bin_count(bin_index)) = 0
                
                    if (i_pseudo.eq.1) then
                    x1r_pseudo(bin_index,bin_count(bin_index)) = x0
                    y1r_pseudo(bin_index,bin_count(bin_index)) = y0
                    x2r_pseudo(bin_index,bin_count(bin_index)) = x0
                    y2r_pseudo(bin_index,bin_count(bin_index)) = y0
                    
                    x1r_pseudo_coord(bin_index,bin_count(bin_index)) = 0
                    y1r_pseudo_coord(bin_index,bin_count(bin_index)) = 0
                    x2r_pseudo_coord(bin_index,bin_count(bin_index)) = 0
                    y2r_pseudo_coord(bin_index,bin_count(bin_index)) = 0
                        
                    endif
                    
                    elseif (i_eddy.eq.1) then
                    x1r_eddy(bin_index,bin_count(bin_index)) = x0
                    y1r_eddy(bin_index,bin_count(bin_index)) = y0
                    x2r_eddy(bin_index,bin_count(bin_index)) = x0
                    y2r_eddy(bin_index,bin_count(bin_index)) = y0
                    
                    x1r_eddy_coord(bin_index,bin_count(bin_index)) = 0
                    y1r_eddy_coord(bin_index,bin_count(bin_index)) = 0
                    x2r_eddy_coord(bin_index,bin_count(bin_index)) = 0
                    y2r_eddy_coord(bin_index,bin_count(bin_index)) = 0
                    endif
                else
                !print*,'bin_count(bin_index) =',bin_count(bin_index)
                bin_stop(bin_index) = 1
                endif
        

        
        index = 0
        do p = 1,nbins
         if (bin_stop(p) .eq. 1) then
          index = index +1
          endif
        enddo  
        !print*,'index=',index
        if (index.ne.nbins) then
            points_count = points_count + 1
        goto 10
        endif
        

        if(i_full.eq.1) then
        
        call write_width(full_name,bin_count,nbins,k)
        call write_bin_boundaries(full_name,k,nbins
     & ,Y_bin)

        endif
	
	     if(i_eddy.eq.1) then
        
         call write_width(eddy_name,bin_count,nbins,k)
         call write_bin_boundaries(eddy_name,k,nbins
     & ,Y_bin)

         endif
	
	     if(i_pseudo.eq.1) then
        
         call write_width(pseudo_name,bin_count,nbins,k)
         call write_bin_boundaries(pseudo_name,k,nbins
     & ,Y_bin)

         endif
         

	
        do p = 1,nbins
            
            if (i_full.eq.1) then
                        
            
            x1 = x1r(p,:)
            x2 = x2r(p,:)
            y1 = y1r(p,:)
            y2 = y2r(p,:)
            
            x1_coord = x1r_coord(p,:)
            x2_coord = x2r_coord(p,:)
            y1_coord = y1r_coord(p,:)
            y2_coord = y2r_coord(p,:)
                        

            call write_binned_file(full_name,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord
     &       ,release_time(k),nrel(k))
     
            endif
     
            if (i_eddy.eq.1) then
            
            x1 = x1r_eddy(p,:)
            x2 = x2r_eddy(p,:)
            y1 = y1r_eddy(p,:)
            y2 = y2r_eddy(p,:)
            
            x1_coord = x1r_eddy_coord(p,:)
            x2_coord = x2r_eddy_coord(p,:)
            y1_coord = y1r_eddy_coord(p,:)
            y2_coord = y2r_eddy_coord(p,:)
            
            call write_binned_file(eddy_name,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord,y2_coord
     &       ,release_time(k),nrel(k))
            endif
            
            if (i_pseudo.eq.1) then
            
            x1 = x1r_pseudo(p,:)
            x2 = x2r_pseudo(p,:)
            y1 = y1r_pseudo(p,:)
            y2 = y2r_pseudo(p,:)
            
            x1_coord = x1r_pseudo_coord(p,:)
            x2_coord = x2r_pseudo_coord(p,:)
            y1_coord = y1r_pseudo_coord(p,:)
            y2_coord = y2r_pseudo_coord(p,:)
            
            call write_binned_file(pseudo_name,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord
     &       ,y2_coord,release_time(k),nrel(k))
            endif
  
        enddo 
        nrel(k) = nrel(k) + 1
        
        write(*,*)'Lagrangian particles randomly generated'

        else
        
c ---------------------------------------------------------------------
c ---------------------- UNIFORMLY BIN THE DOMAIN AND SEED PARTICLES -----------

c bin domain

      write(*,*) 'Binning domain uniformally'

      d_bin = dfloat(jj)/nbins
      
      do p = 1,nbins
        y_bin(p) = p*d_bin
      enddo
      
      do p = 1,nbins
      do n = 1,npoints
      
        x0 = ran1(iseed)*dfloat(ii)
        if (p.eq.1) then
        y0 = ran1(iseed)*d_bin
        else
        y0 = (ran1(iseed)*d_bin)+y_bin(p-1)
        endif
        
        if (i_full.eq.1. .or. i_pseudo.eq.1) then
            x1r(p,n) = x0
            y1r(p,n) = y0
            x2r(p,n) = x0
            y2r(p,n) = y0
            
            x1r_coord(p,n) = 0
            y1r_coord(p,n) = 0
            x2r_coord(p,n) = 0
            y2r_coord(p,n) = 0
        endif
        if (i_eddy.eq.1) then
            x1r_eddy(p,n) = x0
            y1r_eddy(p,n) = y0
            x2r_eddy(p,n) = x0
            y2r_eddy(p,n) = y0
            
            x1r_eddy_coord(p,n) = 0
            y1r_eddy_coord(p,n) = 0
            x2r_eddy_coord(p,n) = 0
            y2r_eddy_coord(p,n) = 0
        endif
        if (i_pseudo .eq. 1) then
            x1r_pseudo(p,n) = x0
            y1r_pseudo(p,n) = y0
            x2r_pseudo(p,n) = x0
            y2r_pseudo(p,n) = y0
            
            x1r_pseudo_coord(p,n) = 0
            y1r_pseudo_coord(p,n) = 0
            x2r_pseudo_coord(p,n) = 0
            y2r_pseudo_coord(p,n) = 0
        endif
      
      enddo
      enddo
    
      
      do p=1,nbins
      
                        if (i_full.eq.1) then
                        
            
            x1 = x1r(p,:)
            x2 = x2r(p,:)
            y1 = y1r(p,:)
            y2 = y2r(p,:)
            
            x1_coord = x1r_coord(p,:)
            x2_coord = x2r_coord(p,:)
            y1_coord = y1r_coord(p,:)
            y2_coord = y2r_coord(p,:)
                        

            call write_binned_file(full_name,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord
     &       ,release_time(k),nrel(k))
     
            endif
     
            if (i_eddy.eq.1) then
            
            x1 = x1r_eddy(p,:)
            x2 = x2r_eddy(p,:)
            y1 = y1r_eddy(p,:)
            y2 = y2r_eddy(p,:)
            
            x1_coord = x1r_eddy_coord(p,:)
            x2_coord = x2r_eddy_coord(p,:)
            y1_coord = y1r_eddy_coord(p,:)
            y2_coord = y2r_eddy_coord(p,:)
            
            call write_binned_file(eddy_name,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord,y2_coord
     &       ,release_time(k),nrel(k))
            endif
            
            if (i_pseudo.eq.1) then
            
            x1 = x1r_pseudo(p,:)
            x2 = x2r_pseudo(p,:)
            y1 = y1r_pseudo(p,:)
            y2 = y2r_pseudo(p,:)
            
            x1_coord = x1r_pseudo_coord(p,:)
            x2_coord = x2r_pseudo_coord(p,:)
            y1_coord = y1r_pseudo_coord(p,:)
            y2_coord = y2r_pseudo_coord(p,:)
            
            call write_binned_file(pseudo_name,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord
     &       ,y2_coord,release_time(k),nrel(k))
            endif
  
        enddo 
        nrel(k) = nrel(k) + 1
        
        
        endif
        
        else


            time_day = time_o(t)*tscale/86400.
            k_s(k) = k_s(k) + 1

        

         
        
c ----- find entry m in time such that time(m) < time_o(k) <= time(m+1)  
        do m = 1,t_len
            if (time(m) >= time_dim(t)) then
            k_new = m
                if(k_new .eq. 2)then
                k_new = 3
                else if (k_new.eq.t_len)then
                k_new = k_new-1
                endif
            exit
            endif
        enddo
c ----- do same for k-1         
        do m = 1,t_len
            if (time(m) >= time_dim(t-1)) then
            k_old = m
                if(k_old .eq. 2)then
                k_old = 3
                elseif(k_old.eq.1)then
                k_old = 3
                else if (k_old.eq.t_len)then
                k_old= k_old-1
                endif
            exit
            endif
        enddo
        
        
        
        time_half = (time_dim(t) + time_dim(t-1))/2
        
        do m = 1,t_len
            if (time(m) >= time_half) then
            k_half = m
                if(k_half .eq. 2)then
                k_half = 3
                else if (k_half.eq.t_len)then
                k_half= k_half-1
                endif
            exit
            endif
        enddo
        

            
            
c ----- cubic interpolate psi in time to find psi
c ----- at time_o(k), time_o(k-1)

c ---- time_o(k)
        
        call read_psi(file_name,ii,jj,k_new-2,psi1,psi2)
  
        
        time_cubic = time(k_new-2:k_new+1)
     
        call interp_time(ii,jj,time_cubic,time_dim(t),psi1,psi1_new)
        call interp_time(ii,jj,time_cubic,time_dim(t),psi2,psi2_new)
        
c ----- time_o(k-1)

        call read_psi(file_name,ii,jj,k_old-2,psi1,psi2)

        
        time_cubic = time(k_old-2:k_old+1)
     
        call interp_time(ii,jj,time_cubic,time_dim(t-1),psi1,psi1_old)
        call interp_time(ii,jj,time_cubic,time_dim(t-1),psi2,psi2_old)
        
c ----- time_half

        call read_psi(file_name,ii,jj,k_half-2,psi1,psi2)
        
        time_cubic = time(k_half-2:k_half+1)
     
        call interp_time(ii,jj,time_cubic,time_half,psi1,psi1_half)
        call interp_time(ii,jj,time_cubic,time_half,psi2,psi2_half)

     
            if (i_eddy.eq.1) then
            
            ! determine eddy streamfunction
            
            do i = 1,ii
            do j = 1,jj
                psi1_eddy_old(i,j) = psi1_old(i,j) - psi1_av(i,j)
                psi2_eddy_old(i,j) = psi2_old(i,j) - psi2_av(i,j)

                
                psi1_eddy_half(i,j) = psi1_half(i,j) - psi1_av(i,j)
                psi2_eddy_half(i,j) = psi2_half(i,j) - psi2_av(i,j)

                
                psi1_eddy_new(i,j) = psi1_new(i,j) - psi1_av(i,j)
                psi2_eddy_new(i,j) = psi2_new(i,j) - psi2_av(i,j)

            enddo
            enddo
            
            endif
            
            !print*,'eddy created'

c --------------- CALCULATE COEFFICIENTS FOR SPATIAL INTERPOLATION -----------  

c --------  BICUBIC INTERPOLATION -----------    
      if (isolve.eq.0) then
      
      if(t.eq.2) then
      
      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then
        
      call A_matrix(ii,jj,psi1_old,M1_old)
      call A_matrix(ii,jj,psi2_old,M2_old)
      
      endif
      
      if (i_eddy.eq.1) then
      
      call A_matrix(ii,jj,psi1_eddy_old,M1_eddy_old)
      call A_matrix(ii,jj,psi2_eddy_old,M2_eddy_old)
      
      endif
      
      else
      
      M1_old = M1_new
      M2_old = M2_new
      
      if (i_eddy.eq.1) then
      
      M1_eddy_old = M1_eddy_new
      M2_eddy_old = M2_eddy_new
      
      endif
  
      endif

      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then
        
      
      call A_matrix(ii,jj,psi1_half,M1_half)
      call A_matrix(ii,jj,psi2_half,M2_half)
      call A_matrix(ii,jj,psi1_new,M1_new)
      call A_matrix(ii,jj,psi2_new,M2_new)
      
      endif
      
      
      if (i_eddy.eq.1) then
      call A_matrix(ii,jj,psi1_eddy_half,M1_eddy_half)
      call A_matrix(ii,jj,psi2_eddy_half,M2_eddy_half)
      call A_matrix(ii,jj,psi1_eddy_new,M1_eddy_new)
      call A_matrix(ii,jj,psi2_eddy_new,M2_eddy_new)
      endif
      
c --------- 2D CUBIC INTERPOLATION 
      
      
      else if (isolve.eq.1) then
      
      if (t.eq.2) then
      
      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then
      
      call cubic_coeff_x(ii,jj,psi1_old
     & ,a1_old,b1_old,c1_old,d1_old)
      call cubic_coeff_x(ii,jj,psi2_old
     & ,a2_old,b2_old,c2_old,d2_old)
     
      endif
      
      if (i_eddy.eq.1) then
      call cubic_coeff_x(ii,jj,psi1_eddy_old
     & ,a1_eddy_old,b1_eddy_old,c1_eddy_old,d1_eddy_old)
      call cubic_coeff_x(ii,jj,psi2_eddy_old
     & ,a2_eddy_old,b2_eddy_old,c2_eddy_old,d2_eddy_old)
      endif
     
      else
      
      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then
      
      a1_old = a1_new
      a2_old = a2_new
      b1_old = b1_new
      b2_old = b2_new
      c1_old = c1_new
      c2_old = c2_new
      d1_old = d1_new
      d2_old = d2_new
      
      endif
      
      if(i_eddy.eq.1) then
      
      a1_eddy_old = a1_eddy_new
      a2_eddy_old = a2_eddy_new
      b1_eddy_old = b1_eddy_new
      b2_eddy_old = b2_eddy_new
      c1_eddy_old = c1_eddy_new
      c2_eddy_old = c2_eddy_new
      d1_eddy_old = d1_eddy_new
      d2_eddy_old = d2_eddy_new
      endif
      
      end if
      
      if ((i_pseudo .eq. 1) .or. (i_full.eq.1)) then
     
    
      call cubic_coeff_x(ii,jj,psi1_half
     &,a1_half,b1_half,c1_half,d1_half)
      call cubic_coeff_x(ii,jj,psi2_half
     &,a2_half,b2_half,c2_half,d2_half)
      call cubic_coeff_x(ii,jj,psi1_new
     &,a1_new,b1_new,c1_new,d1_new)
      call cubic_coeff_x(ii,jj,psi2_new
     &,a2_new,b2_new,c2_new,d2_new)
     
      endif
     
      if(i_eddy.eq.1) then
      call cubic_coeff_x(ii,jj,psi1_eddy_half
     &,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half)
      call cubic_coeff_x(ii,jj,psi2_eddy_half
     &,a2_eddy_half,b2_eddy_half,c2_eddy_half,d2_eddy_half)
      call cubic_coeff_x(ii,jj,psi1_eddy_new
     &,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new)
      call cubic_coeff_x(ii,jj,psi2_eddy_new
     &,a2_eddy_new,b2_eddy_new,c2_eddy_new,d2_eddy_new)
      endif
     
      endif
      
            
   
c ------------------ RK4 - FULL ADVECTION ---------------------------
            
            
c FULL ADVECTION
            
            if ((i_full.eq.1).or.(i_pseudo.eq.1)) then
            do p = 1,nbins
            do n = 1,npoints
            
        if (isolve.eq.0) then
        
        call rk4_bicubic(ii,jj,x1r(p,n),y1r(p,n)
     &   ,dt_nondim,U_0,M1_old,M1_half,M1_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2r(p,n),y2r(p,n)
     &   ,dt_nondim,dfloat(0),M2_old,M2_half,M2_new
     & ,x_diff2,y_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1r(p,n),y1r(p,n)
     &   ,dt_nondim,U_0,a1_old,b1_old,c1_old
     & ,d1_old,a1_half,b1_half,c1_half,d1_half
     & ,a1_new,b1_new,c1_new,d1_new
     & ,x_diff1,y_diff1)
        call rk4_2dcubic(ii,jj,x2r(p,n),y2r(p,n)
     &   ,dt_nondim,dfloat(0),a2_old,b2_old,c2_old
     & ,d2_old,a2_half,b2_half,c2_half,d2_half
     & ,a2_new,b2_new,c2_new,d2_new
     & ,x_diff2,y_diff2)
        
        endif

      if (y_diff1 .gt. 20) then

        print*,'y_diff1 = ',y_diff1
        print*,'psi_old = ',psi1_old
       stop
      endif
        
c MEAN CONTRIBUTION FOR THE PSEUDO TRAJECTORIES

        if(i_pseudo.eq.1) then
        
                 if (isolve.eq.0) then
         
        call rk4_bicubic(ii,jj,x1r(p,n),y1r(p,n),dt_nondim,U_0
     & ,M1_av,M1_av,M1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_bicubic(ii,jj,x2r(p,n),y2r(p,n),dt_nondim,dfloat(0)
     & ,M2_av,M2_av,M2_av
     & ,x_av_diff2,y_av_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1r(p,n),y1r(p,n),dt_nondim,U_0
     & ,a1_av,b1_av,c1_av
     & ,d1_av,a1_av,b1_av,c1_av,d1_av
     & ,a1_av,b1_av,c1_av,d1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_2dcubic(ii,jj,x2r(p,n),y2r(p,n),dt_nondim,dfloat(0)
     & ,a2_av,b2_av,c2_av
     & ,d2_av,a2_av,b2_av,c2_av,d2_av
     & ,a2_av,b2_av,c2_av,d2_av
     & ,x_av_diff2,y_av_diff2)
        
        endif 
        
        ENDIF
        
        !print*, 'full location'
        !print*, x1r(n,t), y1r(n,t)
        
        
        x1r(p,n) = x1r(p,n) + x_diff1
        y1r(p,n) = y1r(p,n) + y_diff1
        x2r(p,n) = x2r(p,n) + x_diff2
        y2r(p,n) = y2r(p,n) + y_diff2
        
        
        
            if(x1r(p,n).lt.0)then
                x1r(p,n) = dfloat(ii)+x1r(p,n)
                x1r_coord(p,n)=x1r_coord(p,n)-1
            end if
            if(x1r(p,n).gt.dfloat(ii))then
                x1r(p,n) = x1r(p,n) - dfloat(ii)
                x1r_coord(p,n) = x1r_coord(p,n)+1
            endif
            if(y1r(p,n).lt.0.)then
              y1r(p,n)=dfloat(jj)+y1r(p,n)
              y1r_coord(p,n) = y1r_coord(p,n) -1
            endif
            if(y1r(p,n).gt.dfloat(jj))then
              y1r(p,n)=y1r(p,n)-dfloat(jj)
              y1r_coord(p,n)= y1r_coord(p,n) +1
            endif
            
            if(x2r(p,n).lt.0)then
                x2r(p,n) = dfloat(ii)+x2r(p,n)
                x2r_coord(p,n)=x2r_coord(p,n)-1
            end if
            if(x2r(p,n).gt.dfloat(ii))then
                x2r(p,n) = x2r(p,n) - dfloat(ii)
                x2r_coord(p,n) = x2r_coord(p,n)+1
            endif
            if(y2r(p,n).lt.0.)then
              y2r(p,n)=dfloat(jj)+y2r(p,n)
              y2r_coord(p,n) = y2r_coord(p,n) -1
            endif
            if(y2r(p,n).gt.dfloat(jj))then
              y2r(p,n)=y2r(p,n)-dfloat(jj)
              y2r_coord(p,n)= y2r_coord(p,n) +1
            endif
            
            
            
            if (i_pseudo.eq.1) then
            
        x1r_pseudo(p,n) = x1r_pseudo(p,n) + x_diff1 - x_av_diff1
        x2r_pseudo(p,n) = x2r_pseudo(p,n) + x_diff2 - x_av_diff2
        y1r_pseudo(p,n) = y1r_pseudo(p,n) + y_diff1 - y_av_diff1
        y2r_pseudo(p,n) = y2r_pseudo(p,n) + y_diff2 - y_av_diff2
        
                    if(x1r_pseudo(p,n).lt.0)then
                x1r_pseudo(p,n) = dfloat(ii)+x1r_pseudo(p,n)
                x1r_pseudo_coord(p,n)=x1r_pseudo_coord(p,n)-1
            end if
            if(x1r_pseudo(p,n).gt.dfloat(ii))then
                x1r_pseudo(p,n) = x1r_pseudo(p,n) - dfloat(ii)
                x1r_pseudo_coord(p,n) = x1r_pseudo_coord(p,n)+1
            endif
            if(y1r_pseudo(p,n).lt.0.)then
              y1r_pseudo(p,n)=dfloat(jj)+y1r_pseudo(p,n)
              y1r_pseudo_coord(p,n) = y1r_pseudo_coord(p,n) -1
            endif
            if(y1r_pseudo(p,n).gt.dfloat(jj))then
              y1r_pseudo(p,n)=y1r_pseudo(p,n)-dfloat(jj)
              y1r_pseudo_coord(p,n)= y1r_pseudo_coord(p,n) +1
            endif
            
            if(x2r_pseudo(p,n).lt.0)then
                x2r_pseudo(p,n) = dfloat(ii)+x2r_pseudo(p,n)
                x2r_pseudo_coord(p,n)=x2r_pseudo_coord(p,n)-1
            end if
            if(x2r_pseudo(p,n).gt.dfloat(ii))then
                x2r_pseudo(p,n) = x2r_pseudo(p,n) - dfloat(ii)
                x2r_pseudo_coord(p,n) = x2r_pseudo_coord(p,n)+1
            endif
            if(y2r_pseudo(p,n).lt.0.)then
              y2r_pseudo(p,n)=dfloat(jj)+y2r_pseudo(p,n)
              y2r_pseudo_coord(p,n) = y2r_pseudo_coord(p,n) -1
            endif
            if(y2r_pseudo(p,n).gt.dfloat(jj))then
              y2r_pseudo(p,n)=y2r_pseudo(p,n)-dfloat(jj)
              y2r_pseudo_coord(p,n)= y2r_pseudo_coord(p,n) +1
            endif
        
        endif
            
            
            enddo
            enddo
            endif
            
            
c ----------------- EDDY INDUCED PARTICLE ADVECTION ---------------
        if (i_eddy.eq.1) then

            do p = 1,nbins
            do n = 1,npoints
            
        if (isolve.eq.0) then
        
        call rk4_bicubic(ii,jj,x1r_eddy(p,n),y1r_eddy(p,n),dt_nondim
     &   ,dfloat(0)
     & ,M1_eddy_old,M1_eddy_half,M1_eddy_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2r_eddy(p,n),y2r_eddy(p,n),dt_nondim
     &   ,dfloat(0)
     & ,M2_eddy_old,M2_eddy_half,M2_eddy_new
     & ,x_diff2,y_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1r_eddy(p,n),y1r_eddy(p,n)
     &   ,dt_nondim,
     &   dfloat(0)
     & ,a1_eddy_old,b1_eddy_old,c1_eddy_old
     & ,d1_eddy_old,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half
     & ,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new
     & ,x_diff1,y_diff1)
     
        !print*, 'x1r_eddy,y1r_eddy'
        !print*, x1r_eddy(p,n), y1r_eddy(p,n)
        call rk4_2dcubic(ii,jj,x2r_eddy(p,n),y2r_eddy(p,n),dt_nondim
     &   ,dfloat(0)
     & ,a2_eddy_old,b2_eddy_old,c2_eddy_old
     & ,d2_eddy_old,a2_eddy_half,b2_eddy_half,c2_eddy_half,d2_eddy_half
     & ,a2_eddy_new,b2_eddy_new,c2_eddy_new,d2_eddy_new
     & ,x_diff2,y_diff2)
        
        endif
            
        x1r_eddy(p,n) = x1r_eddy(p,n) + x_diff1
        y1r_eddy(p,n) = y1r_eddy(p,n) + y_diff1
        x2r_eddy(p,n) = x2r_eddy(p,n) + x_diff2
        y2r_eddy(p,n) = y2r_eddy(p,n) + y_diff2
        
            if(x1r_eddy(p,n).lt.0)then
                x1r_eddy(p,n) = dfloat(ii)+x1r_eddy(p,n)
                x1r_eddy_coord(p,n)=x1r_eddy_coord(p,n)-1
            end if
            if(x1r_eddy(p,n).gt.dfloat(ii))then
                x1r_eddy(p,n) = x1r_eddy(p,n) - dfloat(ii)
                x1r_eddy_coord(p,n) = x1r_eddy_coord(p,n)+1
            endif
            if(y1r_eddy(p,n).lt.0.)then
              y1r_eddy(p,n)=dfloat(jj)+y1r_eddy(p,n)
              y1r_eddy_coord(p,n) = y1r_eddy_coord(p,n) -1
            endif
            if(y1r_eddy(p,n).gt.dfloat(jj))then
              y1r_eddy(p,n)=y1r_eddy(p,n)-dfloat(jj)
              y1r_eddy_coord(p,n)= y1r_eddy_coord(p,n) +1
            endif
            
            if(x2r_eddy(p,n).lt.0)then
                x2r_eddy(p,n) = dfloat(ii)+x2r_eddy(p,n)
                x2r_eddy_coord(p,n)=x2r_eddy_coord(p,n)-1
            end if
            if(x2r_eddy(p,n).gt.dfloat(ii))then
                x2r_eddy(p,n) = x2r_eddy(p,n) - dfloat(ii)
                x2r_eddy_coord(p,n) = x2r_eddy_coord(p,n)+1
            endif
            if(y2r_eddy(p,n).lt.0.)then
              y2r_eddy(p,n)=dfloat(jj)+y2r_eddy(p,n)
              y2r_eddy_coord(p,n) = y2r_eddy_coord(p,n) -1
            endif
            if(y2r_eddy(p,n).gt.dfloat(jj))then
              y2r_eddy(p,n)=y2r_eddy(p,n)-dfloat(jj)
              y2r_eddy_coord(p,n)= y2r_eddy_coord(p,n) +1
            endif
            
            
            enddo
            enddo
            endif

            

             
            if (k_s(k) .eq. l_tot_day*k_save) then
            
            k_s(k) = 0
            
            print *, 'Writing to NETCDF files at time',time_day
            
            
     
            if (i_full.eq.1) then
            
            do p = 1,nbins
            
            x1 = x1r(p,:)
            x2 = x2r(p,:)
            y1 = y1r(p,:)
            y2 = y2r(p,:)
            
            x1_coord = x1r_coord(p,:)
            x2_coord = x2r_coord(p,:)
            y1_coord = y1r_coord(p,:)
            y2_coord = y2r_coord(p,:)

            call write_binned_file(full_name,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord
     &       ,time_day,nrel(k))
     
            enddo
     
            endif
     
            if (i_eddy.eq.1) then
            
            
            do p = 1,nbins
            
            
            x1 = x1r_eddy(p,:)
            x2 = x2r_eddy(p,:)
            y1 = y1r_eddy(p,:)
            y2 = y2r_eddy(p,:)
            
            x1_coord = x1r_eddy_coord(p,:)
            x2_coord = x2r_eddy_coord(p,:)
            y1_coord = y1r_eddy_coord(p,:)
            y2_coord = y2r_eddy_coord(p,:)
            
            call write_binned_file(eddy_name,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord,y2_coord
     &       ,time_day,nrel(k))
            enddo
            endif
            
            
            if (i_pseudo.eq.1) then
            do p = 1,nbins
            
            x1 = x1r_pseudo(p,:)
            x2 = x2r_pseudo(p,:)
            y1 = y1r_pseudo(p,:)
            y2 = y2r_pseudo(p,:)
            
            x1_coord = x1r_pseudo_coord(p,:)
            x2_coord = x2r_pseudo_coord(p,:)
            y1_coord = y1r_pseudo_coord(p,:)
            y2_coord = y2r_pseudo_coord(p,:)
            
            call write_binned_file(pseudo_name,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord
     &       ,y2_coord,time_day,nrel(k))
            enddo
            endif
        
            nrel(k) = nrel(k) + 1
            
            endif
            
            
            endif
            enddo
          
            
        enddo
      
      end program offline_transport
