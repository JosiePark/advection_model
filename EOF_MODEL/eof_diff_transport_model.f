c PROGRAM THAT PERFORMS OFFLINE ADVECTION
c RELEASES PARTICLES UNIFORMLY AT DIFFERENT TIMES

c USES EOFS AS THE INPUT VELOCITY FIELD
c BUT THIS TIME DEDUCTS THE EOFS FROM THE TOTAL STREAM FUNCTION

c NEED TO CALCULATE FFE OR PSEUDO TRAJECTORIES

      program eof_transport_model
      
      use mod_advection_input
      use mod_PVbin_advection_constants
      use mod_qg2_netcdf
      use mod_traj_netcdf
      use mod_2dcubic
      use mod_bicubic
      use mod_random
      use mod_time_interp
      use mod_laplace
      use mod_rk4
      use mod_eof_netcdf
      
      
      
      implicit none 
      
      
c ---- DEFINE DIMENSIONAL VARIABLES
      
      uscale=1.
      scale=basinscale/dfloat(ii) 
      tscale=scale/uscale
      SS=(scale/Rd)**2

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      beta_nondim=beta*scale*scale/uscale
        
c ------ CREATE TRAJECTORY DATA FILE ----------------
        
c ----- Initialise trajectory arrays
        
        allocate(x1r(nbins,npoints), 
     &   x2r(nbins,npoints)
     & ,y1r(nbins,npoints),y2r(nbins,npoints)
     & ,x1r_coord(nbins,npoints)
     & ,x2r_coord(nbins,npoints)
     & ,y1r_coord(nbins,npoints)
     & ,y2r_coord(nbins,npoints))
     
             allocate(x1r_pseudo(nbins,npoints)
     & ,x2r_pseudo(nbins,npoints)
     & ,y1r_pseudo(nbins,npoints),y2r_pseudo(nbins,npoints)
     & ,x1r_pseudo_coord(nbins,npoints)
     & ,x2r_pseudo_coord(nbins,npoints)
     & ,y1r_pseudo_coord(nbins,npoints)
     & ,y2r_pseudo_coord(nbins,npoints)) 
     
      if (i_eddy.eq.1) then 
     
        allocate(x1r_eddy(nbins,npoints)
     & ,x2r_eddy(nbins,npoints)
     & ,y1r_eddy(nbins,npoints),y2r_eddy(nbins,npoints)
     & ,x1r_eddy_coord(nbins,npoints)
     & ,x2r_eddy_coord(nbins,npoints)
     & ,y1r_eddy_coord(nbins,npoints)
     & ,y2r_eddy_coord(nbins,npoints))
        endif
        
        print*,'traj_file', traj_file
        print*,'full_traj_file', full_traj_file
        print*,'eddy_traj_file',eddy_traj_file
        print*,'file_name', file_name
        print*,'ave_name', ave_name
        if (i_start .eq. 0) then
        print*,'i_start = 0'
        call create_binned_file(traj_file,nbins,npoints,release_no)
        call create_binned_file(full_traj_file,nbins,npoints,release_no)
        if (i_eddy .eq. 1) then
        call create_binned_file(eddy_traj_file,nbins,npoints,release_no)
        endif
        release_start = 1
        else
        !!!! EDIT TO READ BOTH LAYERS !!!! 
        !!!! ALSO EDIT TO READ FROM FULL TRAJECTORY FILE !!!!
        call read_binned_file(traj_file,npoints,nbins,x1r_pseudo
     &   ,y1r_pseudo,x2r_pseudo,y2r_pseudo
     &   ,x1r_pseudo_coord,y1r_pseudo_coord
     &   ,x2r_pseudo_coord,y2r_pseudo_coord,start_info)
        call read_binned_file(full_traj_file,npoints,nbins,x1r
     &   ,y1r,x2r,y2r
     &   ,x1r_coord,y1r_coord
     &   ,x2r_coord,y2r_coord,start_info)
        release_start = start_info(1)
        endif

c --- READ TIME AND TIME AVERAGE STREAM FUNCTION ----

	     print*,'file_name=',file_name
      
        call read_time(file_name,time,t_len)
	     print*,'time read'
	     print*,'ave_name =',ave_name
        call read_ave_file(ave_name,ii,jj,psi1_av,psi2_av,time_av)
	     print*,'ave_read'
         
c ---- Calculate coefficients for time-averaged stream function
c ---- Needed for spatial interpolation

      if (isolve.eq.0) then
      
        
      call A_matrix(ii,jj,psi1_av,M1_av)
      call A_matrix(ii,jj,psi2_av,M2_av)
    
      
      else if (isolve.eq.1) then
    
      
      call cubic_coeff_x(ii,jj,psi1_av
     & ,a1_av,b1_av,c1_av,d1_av)
      call cubic_coeff_x(ii,jj,psi2_av
     & ,a2_av,b2_av,c2_av,d2_av)
     
      endif

        
c ----- CALCULATE THE TRUE TIME ARRAY (WHICH TAKES INTO ACCOUNT THE CHOSEN OFFLINE TIME STEP)

        release_length_secs = release_length*86400 ! time recording in seconds for a release
        max_time_lag = time(t_len)*86400
        
        t_tot = int(release_length_secs/dt) ! total number of time steps per release
        l_tot_day = int(86400/dt) ! total number of time steps per day

        dt_nondim = dt/tscale
            
        allocate(time_o(t_tot))
        allocate(time_dim(t_tot))

        allocate(release_time(release_no) ,nrel(release_no))
        allocate(k_s(release_no))
        
        
        if (i_start .eq. 0) then
        do k = 1,release_no    
            release_time(k) = (k-1)*release_interval
            nrel(k) = 1
            k_s(k) = 0
        enddo
        else
        do k = release_start,release_no    
            release_time(k) = (k-1)*release_interval
            nrel(k) = start_info(2)
            k_s(k) = 0
        enddo
        endif
 

cc  LOOP THROUGH THE DIFFERENT TEMPORAL RELEASES 
      do k = release_start,release_no
      
        print*, 'release number =',k
        do t = 1,t_tot
            time_o(t) = (t-1)*dt_nondim + release_time(k)*86400/tscale
            time_dim(t) = time_o(t)*tscale/86400
            
        enddo
        
        if ((i_start .eq. 0) .or. ((i_start .eq. 1) .and. 
     &        (k .ne. release_start))) then
        
c BIN PARTICLES UNIFORMLY


        iseed = 102
        
        nrec = 0
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
        
            x1r(p,n) = x0
            y1r(p,n) = y0
            x2r(p,n) = x0
            y2r(p,n) = y0
            
            x1r_coord(p,n) = 0
            y1r_coord(p,n) = 0
            x2r_coord(p,n) = 0
            y2r_coord(p,n) = 0
            
            x1r_pseudo(p,n) = x0
            y1r_pseudo(p,n) = y0
            x2r_pseudo(p,n) = x0
            y2r_pseudo(p,n) = y0
            
            x1r_pseudo_coord(p,n) = 0
            y1r_pseudo_coord(p,n) = 0
            x2r_pseudo_coord(p,n) = 0
            y2r_pseudo_coord(p,n) = 0
            
            if (i_eddy .eq. 1) then
            
            x1r_eddy(p,n) = x0
            y1r_eddy(p,n) = y0
            x2r_eddy(p,n) = x0
            y2r_eddy(p,n) = y0
            
            x1r_eddy_coord(p,n) = 0
            y1r_eddy_coord(p,n) = 0
            x2r_eddy_coord(p,n) = 0
            y2r_eddy_coord(p,n) = 0
            
            endif

      
      enddo
      enddo
    
      
      do p=1,nbins
      
            x1 = x1r(p,:)
            x2 = x2r(p,:)
            y1 = y1r(p,:)
            y2 = y2r(p,:)
            
            x1_coord = x1r_coord(p,:)
            x2_coord = x2r_coord(p,:)
            y1_coord = y1r_coord(p,:)
            y2_coord = y2r_coord(p,:)
                        

            call write_binned_file(traj_file,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord
     &       ,release_time(k),nrel(k))
     
            call write_binned_file(full_traj_file,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord
     &       ,release_time(k),nrel(k))
     
            if (i_eddy .eq. 1) then
     
            call write_binned_file(eddy_traj_file,npoints,p,k
     &       ,x1,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord,x2_coord
     &       ,y2_coord
     &       ,release_time(k),nrel(k))
     
            endif
     
     
      enddo
      
      write(*,*)'Lagrangian particles randomly generated'
      
      endif
      
      if (i_start .eq. 0) then
        t_start = 1
      else
        if (k .eq. release_start) then
            t_start = start_info(2)
        else
            t_start = 1
        endif
      endif
            
        do t = t_start+1,t_tot
        !print*,'t=',t
        
c ------------ MAIN CYCLE -----------------------------

        k_s(k) = k_s(k) + 1


        time_day = time_o(t)*tscale/86400.
        !k_s = k_s +dt_nondim*tscale/86400.
        !print*, 'k_s=',k_s

         
        
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
        
        call read_eof(eof_file,ii,jj,k_new-2,psi1_eof,psi2_eof)
        !call read_psi(file_name,ii,jj,k_new-2,psi1,psi2)
  
        
        time_cubic = time(k_new-2:k_new+1)
        
        call interp_time(ii,jj,time_cubic,time_dim(t)
     &        ,psi1_eof,psi1_new)
        call interp_time(ii,jj,time_cubic,time_dim(t)
     &        ,psi2_eof,psi2_new)
     
     
        if (i_eddy .eq. 1) then
     
        psi1_eddy_new = psi1_new - psi1_av
        psi2_eddy_new = psi2_new - psi2_av
        
        endif
        !call interp_time(ii,jj,time_cubic,time_dim(t),psi1,psi1_new)
        !call interp_time(ii,jj,time_cubic,time_dim(t),psi2,psi2_new)
        
c ----- time_o(k-1)

        call read_eof(eof_file,ii,jj,k_old-2,psi1_eof,psi2_eof)
        !call read_psi(file_name,ii,jj,k_old-2,psi1,psi2)

        
        time_cubic = time(k_old-2:k_old+1)
        
        call interp_time(ii,jj,time_cubic,time_dim(t-1)
     &        ,psi1_eof,psi1_old)
        call interp_time(ii,jj,time_cubic,time_dim(t-1)
     &        ,psi2_eof,psi2_old)
     
        if (i_eddy .eq. 1) then
     
        psi1_eddy_old = psi1_old - psi1_av
        psi2_eddy_old = psi2_old - psi2_av
        
        endif
     
        !call interp_time(ii,jj,time_cubic,time_dim(t-1),psi1,psi1_old)
        !call interp_time(ii,jj,time_cubic,time_dim(t-1),psi2,psi2_old)
        
c ----- time_half

        call read_eof(eof_file,ii,jj,k_half-2,psi1_eof,psi2_eof)
        !call read_psi(file_name,ii,jj,k_half-2,psi1,psi2)
        
        time_cubic = time(k_half-2:k_half+1)
        
        call interp_time(ii,jj,time_cubic,time_half
     &        ,psi1_eof,psi1_half)
        call interp_time(ii,jj,time_cubic,time_half
     &        ,psi2_eof,psi2_half)
     
        if (i_eddy .eq. 1) then
     
        psi1_eddy_half = psi1_half - psi1_av
        psi2_eddy_half = psi2_half - psi2_av
        
        endif
     
        !call interp_time(ii,jj,time_cubic,time_half,psi1,psi1_half)
        !call interp_time(ii,jj,time_cubic,time_half,psi2,psi2_half)
        


c --------------- CALCULATE COEFFICIENTS FOR SPATIAL INTERPOLATION ----------- 

      if (isolve.eq.0) then
      
      if(k.eq.2) then
        
      call A_matrix(ii,jj,psi1_old,M1_old)
      call A_matrix(ii,jj,psi2_old,M2_old)
      
      if (i_eddy .eq. 1) then
      
      call A_matrix(ii,jj,psi1_eddy_old,M1_eddy_old)
      call A_matrix(ii,jj,psi2_eddy_old,M2_eddy_old)
      
      endif
      
      
      else
      
      M1_old = M1_new
      M2_old = M2_new
      
      if (i_eddy .eq. 1) then
      
      M1_eddy_old = M1_eddy_new
      M2_eddy_old = M2_eddy_new
      
      endif
      
      
      
      endif
      
      call A_matrix(ii,jj,psi1_half,M1_half)
      call A_matrix(ii,jj,psi2_half,M2_half)
      call A_matrix(ii,jj,psi1_new,M1_new)
      call A_matrix(ii,jj,psi2_new,M2_new)
      
      if (i_eddy .eq. 1) then
      
      call A_matrix(ii,jj,psi1_eddy_half,M1_eddy_half)
      call A_matrix(ii,jj,psi2_eddy_half,M2_eddy_half)
      call A_matrix(ii,jj,psi1_eddy_new,M1_eddy_new)
      call A_matrix(ii,jj,psi2_eddy_new,M2_eddy_new)
      
      endif
 
      else if (isolve.eq.1) then
      
      if (k.eq.2) then
      
      call cubic_coeff_x(ii,jj,psi1_old
     & ,a1_old,b1_old,c1_old,d1_old)
      call cubic_coeff_x(ii,jj,psi2_old
     & ,a2_old,b2_old,c2_old,d2_old)
      if (i_eddy.eq.1) then
      call cubic_coeff_x(ii,jj,psi1_eddy_old
     & ,a1_eddy_old,b1_eddy_old,c1_eddy_old,d1_eddy_old)
      call cubic_coeff_x(ii,jj,psi2_eddy_old
     & ,a2_eddy_old,b2_eddy_old,c2_eddy_old,d2_eddy_old)
      endif
     
      else
      
      a1_old = a1_new
      a2_old = a2_new
      b1_old = b1_new
      b2_old = b2_new
      c1_old = c1_new
      c2_old = c2_new
      d1_old = d1_new
      d2_old = d2_new
      
      if (i_eddy .eq. 1) then
      
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
     
    
      call cubic_coeff_x(ii,jj,psi1_half
     &,a1_half,b1_half,c1_half,d1_half)
      call cubic_coeff_x(ii,jj,psi2_half
     &,a2_half,b2_half,c2_half,d2_half)
      call cubic_coeff_x(ii,jj,psi1_new
     &,a1_new,b1_new,c1_new,d1_new)
      call cubic_coeff_x(ii,jj,psi2_new
     &,a2_new,b2_new,c2_new,d2_new)
     
      if (i_eddy .eq. 1) then
     
      call cubic_coeff_x(ii,jj,psi1_eddy_half
     &,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half)
      call cubic_coeff_x(ii,jj,psi2_eddy_half
     &,a2_eddy_half,b2_eddy_half,c2_eddy_half,d2_eddy_half)
      call cubic_coeff_x(ii,jj,psi1_eddy_new
     &,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new)
      call cubic_coeff_x(ii,jj,psi2_eddy_new
     &,a2_new,b2_eddy_new,c2_eddy_new,d2_eddy_new)
     
      endif
     
      endif
      
      !print*,'coefficients created'
            
   
c ------------------ RK4 - FULL ADVECTION ---------------------------
       
c FULL ADVECTION
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
        
c MEAN CONTRIBUTION FOR THE PSEUDO TRAJECTORIES

        
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

c EDDY ONLY TRAJECTORIES

        if (i_eddy .eq. 1) then

        if (isolve.eq.0) then
        
        call rk4_bicubic(ii,jj,x1r_eddy(p,n),y1r_eddy(p,n)
     &   ,dt_nondim,U_0,M1_eddy_old,M1_eddy_half,M1_eddy_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2r_eddy(p,n),y2r_eddy(p,n)
     &   ,dt_nondim,dfloat(0),M2_eddy_old,M2_eddy_half,M2_eddy_new
     & ,x_diff2,y_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1r_eddy(p,n),y1r_eddy(p,n)
     &   ,dt_nondim,U_0,a1_eddy_old,b1_eddy_old,c1_eddy_old
     & ,d1_eddy_old,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half
     & ,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new
     & ,x_diff1,y_diff1)
        call rk4_2dcubic(ii,jj,x2r_eddy(p,n),y2r_eddy(p,n)
     &   ,dt_nondim,dfloat(0),a2_eddy_old,b2_eddy_old,c2_eddy_old
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
            
            endif
        
        

        
            
            
            enddo
            enddo
            

            
             
             if (k_s(k) .eq. l_tot_day*k_save) then
            !nrec= nrec + 1
c            do t=1,release_no
c            print*, 't=',t

c            if((time_dim(k).gt.release_time(t)) .and.
c     &       ((time_dim(k).lt.release_time(t) + release_length)))then 
            print*, 'saving release time', release_time(k)
            
            k_s(k) = 0
            
            print *, 'Writing to NETCDF files at time',time_day
            !print *, 'Writing release number',t
            
     
            !if (i_full.eq.1) then
            
            do p = 1,nbins
            
            x1 = x1r_pseudo(p,:)
            x2 = x2r_pseudo(p,:)
            y1 = y1r_pseudo(p,:)
            y2 = y2r_pseudo(p,:)
            
            x1_coord = x1r_pseudo_coord(p,:)
            x2_coord = x2r_pseudo_coord(p,:)
            y1_coord = y1r_pseudo_coord(p,:)
            y2_coord = y2r_pseudo_coord(p,:)

            call write_binned_file(traj_file,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord,y2_coord
     &       ,time_day,nrel(k))
     
            x1 = x1r(p,:)
            x2 = x2r(p,:)
            y1 = y1r(p,:)
            y2 = y2r(p,:)
            
            x1_coord = x1r_coord(p,:)
            x2_coord = x2r_coord(p,:)
            y1_coord = y1r_coord(p,:)
            y2_coord = y2r_coord(p,:)
            
            call write_binned_file(full_traj_file,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord,y2_coord
     &       ,time_day,nrel(k))
     
            if (i_eddy .eq. 1) then
     
            x1 = x1r_eddy(p,:)
            x2 = x2r_eddy(p,:)
            y1 = y1r_eddy(p,:)
            y2 = y2r_eddy(p,:)
            
            x1_coord = x1r_eddy_coord(p,:)
            x2_coord = x2r_eddy_coord(p,:)
            y1_coord = y1r_eddy_coord(p,:)
            y2_coord = y2r_eddy_coord(p,:)
            
            call write_binned_file(eddy_traj_file,npoints,p,k
     &       ,x1
     &       ,y1
     &       ,x2,y2
     &       ,x1_coord,y1_coord
     &       ,x2_coord,y2_coord
     &       ,time_day,nrel(k))
     
            endif
            
            enddo
            nrel(k) = nrel(k) + 1
            
            endif
            
            
    
            enddo
          
        
        enddo
      
      end program eof_transport_model
