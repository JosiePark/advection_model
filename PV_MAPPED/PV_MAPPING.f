C CODE THAT CALCULATES THE PV MAPPED DISPERSION

      program PV_MAPPING
      
      use mod_PVmap_constants
      use mod_pvmap_input
      use mod_qg2_netcdf
      use mod_variables
      use mod_traj_netcdf
      use mod_time_interp
      use mod_2dcubic
      use mod_1dinterp
      use mod_PVdisp_netcdf
      
      implicit none
      
c set non-dimensional variables

      uscale=1.
      scale=basinscale/dfloat(ii) 
      tscale=scale/uscale
      SS=(scale/Rd)**2

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      beta_nondim=beta*scale*scale/uscale
      
      BETA_NONDIM_U1=BETA_NONDIM+U_0*S1
      BETA_NONDIM_U2=BETA_NONDIM-U_0*S2
      
c grid on which PV is saved
      
      allocate(x_c(ii),y_c(jj1))

      do i = 1,ii
        x_c(i) = dfloat(i - 1)
      enddo
      
      do j = 1,jj1
        y_c(j) = dfloat(j - 1)
      enddo
      
     
c reading data

      call read_ave_file(ave_file,ii,jj,psi1_av,psi2_av,time_av) !  READ TIME-AVERAGED STREAM FUNCTION
      call rel_from_psi(ii,jj,psi1_av,psi2_av,rel1_av,rel2_av) ! calculate time-averaged relative vorticity
      call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1_av,psi2_av,rel1_av
     & ,rel2_av,zeta1_av,zeta2_av) ! time averaged PV anomaly
     
      print*,'time-averaged stream function read and zeta calculated'
c CALCULATE ZONAL FULL TIME-AVERAGED PV

      do j = 1,int((coord_range+1)*jj)
            beta1_y(j) = beta_nondim_u1*(dfloat(j-1)) ! only do for the top layer
      enddo
      
      do j=1,jj
      do kk = 1,coord_range+1
            idx = j+(kk-1)*jj
            PV_bar(idx) = 0.
            do i = 1,ii
                PV_av(i,idx) = zeta1_av(i,j) + beta1_y(idx)
                PV_bar(idx) = PV_bar(idx) + PV_av(i,idx)
            enddo
            PV_bar(idx) = PV_bar(idx)/dfloat(ii)
      enddo
      enddo
      
      
c READ TRAJECTORY DATA ARRAY SIZES
      
      call read_time(qg_file,qg_time,qg_t_len)  ! read time array from dynamical model data file  
      call read_binned_file_dimensions(traj_file,nbins,nrel
     & ,npoints,traj_t_len)                                ! reads properties of trajectoriy data file
      print*,'nbins = ',nbins
      print*,'nrel=',nrel
      allocate(x1(npoints),y1(npoints)
     & ,x1_coord(npoints),y1_coord(npoints),y_loc(npoints))
     
      allocate(traj_time(nrel,traj_t_len))
      
      allocate(bin_width(nbins,nrel),bin_boundaries(nbins,nrel))
      
      call read_binned_file_time(traj_file,nrel,nbins
     & ,traj_t_len,traj_time)
       print*,'time read'
      call read_release_bin_width(bin_file,nbins,nrel,bin_width)
       print*,'bin_width read'
      do k = 1,nrel
      
      bin_boundaries(1,k) = bin_width(1,k)
      
      do p = 2,nbins
      
        bin_boundaries(p,k) = bin_boundaries(p-1,k) + bin_width(p,k)
      
      enddo
      enddo
      
      call create_binned_PVdisp_file(disp_file
     & ,npoints,nbins,nrel,bin_boundaries)
      
      allocate(y0_map(npoints),y_map(npoints)
     & ,PV_0(npoints),PV_t(npoints))
      
      do p=1,nbins
      print*,'starting bin =',p
      do k = 1,nrel
      print*,'release = ',k
      do t= 1,traj_t_len
      !do t=1,50
      
      print*,'time =',t
    
      
c read trajectory data from release file 
        call read_PV_binned_trajectories(traj_file,npoints,k,p
     & ,t
     & ,x1,y1
     & ,x1_coord,y1_coord )
c remove non-entries (set to zero)

         do n =1,npoints
            if (x1(n).lt.0) then
            x1(n) = 0
            endif
            if (y1(n).lt.0) then
            y1(n) = 0
            endif

            
            if (x1(n).gt.ii) then
            x1(n) = 0
            endif
            if (y1(n).gt.jj) then
            y1(n) = 0
            endif

            
            if (x1_coord(n).lt.-100) then
            x1_coord(n) = 0
            endif
            if (y1_coord(n).lt.-100) then
            y1_coord(n) = 0
            endif

            
            if (x1_coord(n).gt.100) then
            x1_coord(n) = 0
            endif
            if (y1_coord(n).gt.100) then
            y1_coord(n) = 0
            endif

            
         enddo
         

c calculate real y location - across the range of domain coordinates
        do n =1,npoints
         y_loc(n) = y1(n) +dfloat(jj)*dfloat(y1_coord(n)-coord_min)
        enddo         

c ----- find entry m in time such that time_qg(m) < time_traj(k_index) <= time_qg(m+1) 
            exact = 0
        do m = 1,qg_t_len
        
            if(qg_time(m).eq.traj_time(k,t)) then

            exact = 1
            k_index = m
     
            exit
     
            elseif (qg_time(m) > traj_time(k,t)) then
            k_index = m
                if(k_index .eq. 2)then
                k_index = 3
                else if (k_index.eq.qg_t_len)then
                k_index = k_index-1
                endif
            exit
            endif
        enddo
        
c read the 4 stream functions around the trajectory time
c and interpolate in time if necessary

        if (exact .eq. 0) then
        call read_psi(qg_file,ii,jj,k_index-2,psi1,psi2)
        
        
        
        time_cubic = qg_time(k_index-2:k_index+1)
        
        call interp_time(ii,jj,time_cubic,traj_time(k,t),psi1
     & ,psi1_interp)
        call interp_time(ii,jj,time_cubic,traj_time(k,t),psi2
     & ,psi2_interp)
 
        elseif(exact.eq.1) then
        
        call read_psi(qg_file,ii,jj,k_index,psi1,psi2)
        
        psi1_interp = psi1(:,:,1)
        psi2_interp = psi2(:,:,1)

        
        endif
        
        if ((sim .eq. 2) .or. (sim .eq. 3)) then
        
            do i = 1,ii
            do j = 1,jj
        
                psi1_interp(i,j) = psi1_interp(i,j) - psi1_av(i,j) ! take eddying streamfunction
                psi2_interp(i,j) = psi2_interp(i,j) - psi2_av(i,j)
            
            enddo
            enddo
            print*,'Eddying streamfunction taken'
        
        endif
c now calculate the full PV at that time

        call rel_from_psi(ii,jj,psi1_interp,psi2_interp,rel1,rel2)
c PV anomaly
        call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1_interp,psi2_interp
     &   ,rel1,rel2,zeta1,zeta2)

c Calculate instantaneous full PV

        do j=1,jj
        do kk = 1,coord_range+1
        do i = 1,ii
            idx = j+(kk-1)*jj
        
            PV(i,idx) = zeta1(i,j) + beta1_y(idx)
        
        enddo
        enddo
        enddo 
        
        call cubic_coeff_x(ii,jj1
     &    ,PV
     &    ,a,b,c,d)      ! calculate coefficients required for interpolation   
      
      do n = 1,npoints
      
        call cubic_poly_x(ii,jj1,x1(n) ! construct polynomial PV_x
     &    ,y_loc(n),a,b,c,d,PV_x)
     
        call cubic_interp(ii,jj1,PV_x,y_loc(n),PV_t(n)) 
        
        
        ! find the four indices in PV_bar which surround PV_t

        do j = 1,jj1
       
            if(PV_bar(j).ge.PV_t(n)) then
                k_index = j
                exit
            endif
        enddo
        
        
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
        else
        do j = 1,4
            PV_x(j) = PV_bar(k_index+j-3)
            y_index(j)= y_c(k_index+j-3)
        enddo
        endif
        
        ! INTERPOLATE PV_x WHICH CONTAINS PV_bar AT THE 4 POINTS SURROUNDING PV_t
        ! TO FIND Y CORRESPONDING TO PV_t.
        
        call interp_1d(PV_x,4,y_index,PV_t(n),y_map(n))
        ! DEBUG
        if(abs(y_map(n)) .gt. 10000) then
        stop
        endif
        
      
      enddo
      
      call write_binned_PVdisp_file(disp_file,npoints,k,p,y_map,PV_t,t)

      enddo
      enddo
      enddo

      end program PV_MAPPING
