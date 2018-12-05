C CODE THAT CALCULATES BOUNDARIES OF PV BINS
C AND ALSO CALCULATES THE WIDTH BY UNIFORMLY SEEDING PARTICLES AND DIVIDING 
C THEM INTO THEIR RESPECTIVE BINS

      program create_bins
      
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
        eddy_name = 'eddy_PV_bins_trajectories.nc'
        call create_binned_file(eddy_name,nbins,npoints,release_no)
      endif
      if(i_pseudo .eq. 1) then
         pseudo_name = 'pseudo_PV_bins_trajectories.nc'
         call create_binned_file(pseudo_name,nbins,npoints,release_no)
      endif
      if (i_full .eq. 1) then
          full_name = 'full_PV_bins_trajectories_boundaries.nc'
          call create_binned_file(full_name,nbins,npoints,release_no)
      endif
      
      else
    
            if (i_eddy .eq. 1) then
        eddy_name = 'eddy_uniform_bins_trajectories.nc'
        call create_binned_file(eddy_name,nbins,npoints,release_no)
      endif
      if(i_pseudo .eq. 1) then
         pseudo_name = 'pseudo_uniform_bins_trajectories.nc'
         call create_binned_file(pseudo_name,nbins,npoints,release_no)
      endif
      if (i_full .eq. 1) then
          full_name = 'full_uniform_bins_trajectories.nc'
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
        print*,'npoints',npoints
        
        nrec = 0
        do k = 1,release_no
        
        print*, 'Starting particle advection for release number',k
        

        do t = 1,t_tot
            time_o(t) = (t-1)*dt_nondim + release_time(k)*86400/tscale
c            if (time_o(k).gt.(max_time/tscale)) then
c            time_o(k) = time_o(k) - max_time/tscale
c            endif 
            time_dim(t) = time_o(t)*tscale/86400
            
        enddo
 
      t = 1
      
c ------- CONSTRUCT BINS UNIFORM IN PV --------------------
      
c Calculate zonally time-averaged full PV across whole domain
c Add extra domain above and blow as the PV snapshot may fall outside of range of PV_bar
      
      call PV_bar_from_psi(ii,jj,jj1,basinscale,beta
     & ,Rd,H1,H2,U_0,psi1_av,psi2_av,PV_bar)
     
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

c ----- cubic interpolate psi in time to find psi
c ----- at time_o(k), time_o(k-1)

c ---- time_o(k)
        
        call read_psi(file_name,ii,jj,k_new-2,psi1,psi2)
        
  
        
        time_cubic = time(k_new-2:k_new+1)
     
        call interp_time(ii,jj,time_cubic,time_dim(t),psi1,psi1_new)
        call interp_time(ii,jj,time_cubic,time_dim(t),psi2,psi2_new) ! psi1_new is the stream function
        
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
                PV(i,j) = zeta1_av(i,j) + beta1_y(j)
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
        
        
c Find the four indices of PV_bar which surround PV_bin

      do j = 1,jj1
      
      if (PV_bar(j).ge.PV_bin(p)) then
        k_index = j
        exit
      endif

      enddo
      
c treat PV at the end of the domain in order to perform interpolation
      
        if(k_index .eq. 1) then
        do j =1,4
            PV_x(j) = PV_bar(j)
            y_index(j) = j_bin(j)
        enddo
        elseif(k_index .eq. 2) then
        do j = 1,4
            PV_x(j) = PV_bar(j)
            y_index(j) = j_bin(j)
        enddo
        elseif(k_index.eq.jj1) then
        do j = 1,4
            PV_x(j) = PV_bar(jj1-4+j)
            y_index(j) = j_bin(jj1-4+j)
        enddo
        else
        do j = 1,4
            PV_x(j) = PV_bar(k_index+j-3)
            y_index(j)= j_bin(k_index+j-3)
        enddo
        endif
    
      
        
c Map PV_bin to Y_bin
        call interp_1d(PV_x,4,y_index,PV_bin(p),Y_bin(p))
        

      enddo
      
      print*,'bin boundaries =',Y_bin
      print*,'PV boundaries = ',PV_bin
      
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
     
        call cubic_interp(ii,jj,PV_x,y0,PV_0) 
        
        ! find the four indices in PV_bar which surround PV_t

        do j = 1,jj1
       
            if(PV_bar(j).ge.PV_0) then
                k_index = j
                exit
            endif
        enddo
        if (PV_0.gt.maxval(PV_bar)) then
        k_index = jj
        endif
        
        
        if(k_index .eq. 1) then
        do j =1,4
            PV_x(j) = PV_bar(j)
            y_index(j) = y_c(j)
        enddo
        elseif(k_index .eq. 2) then
        do j = 1,4
            PV_x(j) = PV_bar(j)
            y_index(j) = y_c(j)
        enddo
        elseif(k_index .eq . jj1) then
        do j =1,4
            PV_x(j) = PV_bar(k_index+j-4)
            y_index(j)= y_c(k_index+j-4)
        enddo
        else
        do j = 1,4
            PV_x(j) = PV_bar(k_index+j-3)
            y_index(j)= y_c(k_index+j-3)
        enddo
        endif
        
C MAP TO Y
        
        ! INTERPOLATE PV_x WHICH CONTAINS PV_bar AT THE 4 POINTS SURROUNDING PV_t
        ! TO FIND Y CORRESPONDING TO PV_t.
        
        call interp_1d(PV_x,4,y_index,PV_0,y_map)
        ! DEBUG
        if(abs(y_map) .gt. 10000) then
        print*,'error,y_map = ',y_map
c        if (PV_0.gt.maxval(PV_bar)) then
c            bin_index = nbins
c            go to 11
c        elseif (PV_0.gt.minval(PV_bar)) then
c            bin_index = 1
c            go to 11
c        else
            print*,'stop code as PV interpolation failed'
            print*,'min,max PV = ',minval(PV),maxval(PV)
            print*,'PV_0 = ',PV_0
            stop
        endif
        !print*,'y_map = ',y_map
        
        print*,'particle location = ',x0,y0
        print*,'y_map = ',y_map
        print*,'PV_0 = ',PV_0
        
      
C BIN

        do p = 1,nbins
            if (y_map.lt.Y_bin(p)) then
                bin_index = p
                go to 11
            endif
        enddo
        
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
         enddo
      
      end program create_bins
