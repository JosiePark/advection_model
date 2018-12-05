c PROGRAM THAT PERFORMS OFFLINE ADVECTION
c RELEASES PARTICLES UNIFORMLY AT DIFFERENT TIMES

      program offline_transport
      
      use mod_advection_input
      use mod_advection_constants
      use mod_variables
      use mod_qg2_netcdf
      use mod_traj_netcdf
      use mod_2dcubic
      use mod_bicubic
      use mod_random
      use mod_time_interp
      use mod_laplace
      use mod_rk4
      
      implicit none 
      
        scale = basinscale/dfloat(ii)
        uscale = 1
        tscale = scale/uscale
        
c ------ CREATE TRAJECTORY DATA FILE ----------------

        if(i_full.eq.1.) then
      call name_files(home_dir,regime,1,file_name,full_name
     &       ,disp_file,ave_name)
      call create_release_file(full_name,npoints,release_no)
      endif
      if (i_pseudo.eq.1) then
      call name_files(home_dir,regime,3,file_name,pseudo_name
     &       ,disp_file,ave_name)
      call create_release_file(pseudo_name,npoints,release_no)
      endif
      if (i_eddy.eq.1) then
      call name_files(home_dir,regime,2,file_name,eddy_name
     &       ,disp_file,ave_name)
       call create_release_file(eddy_name,npoints,release_no)
      endif

 
c --- READ TIME AND TIME AVERAGE STREAM FUNCTION ----
      
      call read_time(file_name,time,t_len)
      call read_ave_file(ave_name,ii,jj,psi1_av,psi2_av,time_av)


c ----- CALCULATE THE TRUE TIME ARRAY (WHICH TAKES INTO ACCOUNT THE CHOSEN OFFLINE TIME STEP)

        release_length_secs = release_length*86400 ! time recording in seconds for a release
        max_time_lag = time(t_len)*86400
        
        t_tot = int(release_length_secs/dt) ! total number of time steps per release
        l_tot_day = int(86400/dt) ! total number of time steps per day
        print*, 't_tot=',t_tot
        dt_nondim = dt/tscale
            
        allocate(time_o(t_tot))
        allocate(time_dim(t_tot))

        
c ----- Initialise trajectory arrays
        
        if((i_full.eq.1).or.(i_pseudo.eq.1)) then
        allocate(x1r(npoints,release_no), x2r(npoints,release_no)
     & ,y1r(npoints,release_no),y2r(npoints,release_no)
     & ,x1r_coord(npoints,release_no),x2r_coord(npoints,release_no)
     & ,y1r_coord(npoints,release_no),y2r_coord(npoints,release_no))
        endif
        
        if (i_eddy.eq.1) then
        allocate(x1r_eddy(npoints,release_no)
     & ,x2r_eddy(npoints,release_no)
     & ,y1r_eddy(npoints,release_no),y2r_eddy(npoints,release_no)
     & ,x1r_eddy_coord(npoints,release_no)
     & ,x2r_eddy_coord(npoints,release_no)
     & ,y1r_eddy_coord(npoints,release_no)
     & ,y2r_eddy_coord(npoints,release_no))
        endif
        
        if (i_pseudo.eq.1) then        
        allocate(x1r_pseudo(npoints,release_no)
     & ,x2r_pseudo(npoints,release_no)
     & ,y1r_pseudo(npoints,release_no),y2r_pseudo(npoints,release_no)
     & ,x1r_pseudo_coord(npoints,release_no)
     & ,x2r_pseudo_coord(npoints,release_no)
     & ,y1r_pseudo_coord(npoints,release_no)
     & ,y2r_pseudo_coord(npoints,release_no)
     & ,x1r_disp(npoints,release_no),y1r_disp(npoints,release_no)
     & ,x2r_disp(npoints,release_no),y2r_disp(npoints,release_no)
     & ,x1r_mean(npoints,release_no),y1r_mean(npoints,release_no)
     & ,x2r_mean(npoints,release_no),y2r_mean(npoints,release_no))
        endif    
        
        allocate(release_time(release_no) ,nrel(release_no))
        allocate(k_s(release_no))
        
        
        do k = 1,release_no    
            release_time(k) = (k-1)*release_interval
            nrel(k) = 1
            k_s(k) = 0
        enddo
        
        
c ------ Retrieve time averaged strem function
        
      if ((i_pseudo.eq.1).or.(i_eddy.eq.1)) then
      

      
c ---- Calculate coefficients for time-averaged stream function
c ---- Needed for spatial interpolation


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
      
c ------- CONSTRUCT BINS UNIFORM IN PV --------------------

      if(i_bin .eq. 1) then
      
c Calculate zonally time-averaged full PV across whole domain
      
      call PV_bar_from_psi(ii,jj,basinscale,beta
     & ,Rd,H1,H2,U_0,psi1_av,psi2_av,PV_bar)
     
c Bin the PV
     
      d_bin = (max(PV_bar)-min(PV_bar))/dfloat(nbins)
      
      do j = 1,jj
      j_bin(j) = j
      enddo
      
     
      do k = 1,nbins

        PV_bin(k) = min(PV_bar) + dfloat(k)*d_bin
c Map PV_bin to Y_bin
        call interp_1d(j_bin,jj,PV_bar,PV_bin(k),Y_bin(k))
        

      enddo
      
      print*,'PV_bin = ',PV_bin
      print*,'Y_bin = ',Y_bin
       
      
      endif
      

      
c ------------ LAGRANGIAN PARTICLES ------------------


        iseed = 102
        
        nrec = 0
        do k = 1,release_no
            n = 0

            do j = 1,npoints_sqrt
            do i = 1,npoints_sqrt
                n = n+1
                x0(n) = ran1(iseed)*dfloat(ii)
                y0(n) = ran1(iseed)*dfloat(jj)
                !y0(n) = ran1(iseed)*dfloat(300)+100
            
            enddo
            enddo
            
            
        
            do n = 1,npoints
                if((i_full.eq.1).or.(i_pseudo.eq.1)) then
        x1r(n,k) = x0(n)
        y1r(n,k) = y0(n)
        x2r(n,k) = x0(n)
        y2r(n,k) = y0(n)
        x1r_coord(n,k) = 0
        x2r_coord(n,k) = 0
        y1r_coord(n,k) = 0
        y2r_coord(n,k) = 0 
                endif
                if (i_eddy.eq.1) then
        x1r_eddy(n,k) = x0(n)
        y1r_eddy(n,k) = y0(n)
        x2r_eddy(n,k) = x0(n)
        y2r_eddy(n,k) = y0(n)
        x1r_eddy_coord(n,k) = 0
        y1r_eddy_coord(n,k) = 0
        x2r_eddy_coord(n,k) = 0
        y2r_eddy_coord(n,k) = 0
                endif
                if (i_pseudo.eq.1) then
        x1r_pseudo(n,k) = x0(n)
        y1r_pseudo(n,k) = y0(n)
        x2r_pseudo(n,k) = x0(n)
        y2r_pseudo(n,k) = y0(n)
        x1r_pseudo_coord(n,k) = 0
        y1r_pseudo_coord(n,k) = 0
        x2r_pseudo_coord(n,k) = 0
        y2r_pseudo_coord(n,k) = 0
                endif
            enddo 
                        if (i_full.eq.1) then

            call write_release_file(full_name,npoints,release_no
     &       ,x1r(:,k),y1r(:,k)
     &       ,x2r(:,k),y2r(:,k)
     &       ,x1r_coord(:,k),y1r_coord(:,k),x2r_coord(:,k)
     &       ,y2r_coord(:,k)
     &       ,k,release_time(k),nrel(k))
     
            endif
     
            if (i_eddy.eq.1) then
            call write_release_file(eddy_name,npoints,release_no
     &       ,x1r_eddy(:,k)
     &       ,y1r_eddy(:,k)
     &       ,x2r_eddy(:,k),y2r_eddy(:,k)
     &       ,x1r_eddy_coord(:,k),y1r_eddy_coord(:,k)
     &       ,x2r_eddy_coord(:,k),y2r_eddy_coord(:,k)
     &       ,k,release_time(k),nrel(k))
            endif
            
            if (i_pseudo.eq.1) then
            call write_release_file(pseudo_name,npoints,release_no
     &       ,x1r_pseudo(:,k)
     &       ,y1r_pseudo(:,k)
     &       ,x2r_pseudo(:,k),y2r_pseudo(:,k)
     &       ,x1r_pseudo_coord(:,k),y1r_pseudo_coord(:,k)
     &       ,x2r_pseudo_coord(:,k)
     &       ,y2r_pseudo_coord(:,k),k,release_time(k),nrel(k))
            endif
            
            nrel(k) = nrel(k) + 1
            
            
        enddo 
        
        write(*,*)'Lagrangian particles randomly generated'

        
c ------------ MAIN CYCLE -----------------------------
        do t=1,release_no
        print*, 'release number',t
        do k = 1,t_tot
            time_o(k) = (k-1)*dt_nondim + release_time(t)*86400/tscale
            if (time_o(k).gt.(max_time/tscale)) then
            time_o(k) = time_o(k) - max_time/tscale
            endif 
            time_dim(k) = time_o(k)*tscale/86400
            
        enddo
 
        do k = 2,t_tot
        k_s(t) = k_s(t) + 1


        time_day = time_o(k)*tscale/86400.
        !k_s = k_s +dt_nondim*tscale/86400.
        !print*, 'k_s=',k_s

         
        
c ----- find entry m in time such that time(m) < time_o(k) <= time(m+1)  
        do m = 1,t_len
            if (time(m) >= time_dim(k)) then
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
            if (time(m) >= time_dim(k-1)) then
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
        
        
        
        time_half = (time_dim(k) + time_dim(k-1))/2
        
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
     
        call interp_time(ii,jj,time_cubic,time_dim(k),psi1,psi1_new)
        call interp_time(ii,jj,time_cubic,time_dim(k),psi2,psi2_new)
        
c ----- time_o(k-1)

        call read_psi(file_name,ii,jj,k_old-2,psi1,psi2)

        
        time_cubic = time(k_old-2:k_old+1)
     
        call interp_time(ii,jj,time_cubic,time_dim(k-1),psi1,psi1_old)
        call interp_time(ii,jj,time_cubic,time_dim(k-1),psi2,psi2_old)
        
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
      if (isolve.eq.0) then
      
      if(k.eq.2) then
        
      call A_matrix(ii,jj,psi1_old,M1_old)
      call A_matrix(ii,jj,psi2_old,M2_old)
      
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
      
      call A_matrix(ii,jj,psi1_half,M1_half)
      call A_matrix(ii,jj,psi2_half,M2_half)
      call A_matrix(ii,jj,psi1_new,M1_new)
      call A_matrix(ii,jj,psi2_new,M2_new)
      
      
      if (i_eddy.eq.1) then
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
     
    
      call cubic_coeff_x(ii,jj,psi1_half
     &,a1_half,b1_half,c1_half,d1_half)
      call cubic_coeff_x(ii,jj,psi2_half
     &,a2_half,b2_half,c2_half,d2_half)
      call cubic_coeff_x(ii,jj,psi1_new
     &,a1_new,b1_new,c1_new,d1_new)
      call cubic_coeff_x(ii,jj,psi2_new
     &,a2_new,b2_new,c2_new,d2_new)
     
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
      
      !print*,'coefficients created'
            
   
c ------------------ RK4 - FULL ADVECTION ---------------------------

c        do t=1,release_no

c            if((time_dim(k).gt.release_time(t)) .and.
c     &       ((time_dim(k).lt.release_time(t) + release_length)))then 
            !print*, 'advecting release time', release_time(t)
            
            
c FULL ADVECTION
            
            if ((i_full.eq.1).or.(i_pseudo.eq.1)) then
            do n = 1,npoints
            
        if (isolve.eq.0) then
        
        call rk4_bicubic(ii,jj,x1r(n,t),y1r(n,t)
     &   ,dt_nondim,U_0,M1_old,M1_half,M1_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2r(n,t),y2r(n,t)
     &   ,dt_nondim,U_0,M2_old,M2_half,M2_new
     & ,x_diff2,y_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1r(n,t),y1r(n,t)
     &   ,dt_nondim,U_0,a1_old,b1_old,c1_old
     & ,d1_old,a1_half,b1_half,c1_half,d1_half
     & ,a1_new,b1_new,c1_new,d1_new
     & ,x_diff1,y_diff1)
        call rk4_2dcubic(ii,jj,x2r(n,t),y2r(n,t)
     &   ,dt_nondim,U_0,a2_old,b2_old,c2_old
     & ,d2_old,a2_half,b2_half,c2_half,d2_half
     & ,a2_new,b2_new,c2_new,d2_new
     & ,x_diff2,y_diff2)
        
        endif
        
c MEAN CONTRIBUTION FOR THE PSEUDO TRAJECTORIES

        if(i_pseudo.eq.1) then
        
                 if (isolve.eq.0) then
         
        call rk4_bicubic(ii,jj,x1r(n,t),y1r(n,t),dt_nondim,U_0
     & ,M1_av,M1_av,M1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_bicubic(ii,jj,x2r(n,t),y2r(n,t),dt_nondim,U_0
     & ,M2_av,M2_av,M2_av
     & ,x_av_diff2,y_av_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1r(n,t),y1r(n,t),dt_nondim,U_0
     & ,a1_av,b1_av,c1_av
     & ,d1_av,a1_av,b1_av,c1_av,d1_av
     & ,a1_av,b1_av,c1_av,d1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_2dcubic(ii,jj,x2r(n,t),y2r(n,t),dt_nondim,U_0
     & ,a2_av,b2_av,c2_av
     & ,d2_av,a2_av,b2_av,c2_av,d2_av
     & ,a2_av,b2_av,c2_av,d2_av
     & ,x_av_diff2,y_av_diff2)
        
        endif 
        
        ENDIF
        
        !print*, 'full location'
        !print*, x1r(n,t), y1r(n,t)
        
        
        x1r(n,t) = x1r(n,t) + x_diff1
        y1r(n,t) = y1r(n,t) + y_diff1
        x2r(n,t) = x2r(n,t) + x_diff2
        y2r(n,t) = y2r(n,t) + y_diff2
        
        
            if(x1r(n,t).lt.0)then
                x1r(n,t) = dfloat(ii)+x1r(n,t)
                x1r_coord(n,t)=x1r_coord(n,t)-1
            end if
            if(x1r(n,t).gt.dfloat(ii))then
                x1r(n,t) = x1r(n,t) - dfloat(ii)
                x1r_coord(n,t) = x1r_coord(n,t)+1
            endif
            if(y1r(n,t).lt.0.)then
              y1r(n,t)=dfloat(jj)+y1r(n,t)
              y1r_coord(n,t) = y1r_coord(n,t) -1
            endif
            if(y1r(n,t).gt.dfloat(jj))then
              y1r(n,t)=y1r(n,t)-dfloat(jj)
              y1r_coord(n,t)= y1r_coord(n,t) +1
            endif
            
            if(x2r(n,t).lt.0)then
                x2r(n,t) = dfloat(ii)+x2r(n,t)
                x2r_coord(n,t)=x2r_coord(n,t)-1
            end if
            if(x2r(n,t).gt.dfloat(ii))then
                x2r(n,t) = x2r(n,t) - dfloat(ii)
                x2r_coord(n,t) = x2r_coord(n,t)+1
            endif
            if(y2r(n,t).lt.0.)then
              y2r(n,t)=dfloat(jj)+y2r(n,t)
              y2r_coord(n,t) = y2r_coord(n,t) -1
            endif
            if(y2r(n,t).gt.dfloat(jj))then
              y2r(n,t)=y2r(n,t)-dfloat(jj)
              y2r_coord(n,t)= y2r_coord(n,t) +1
            endif
            
            
            
            if (i_pseudo.eq.1) then
            
        x1r_pseudo(n,t) = x1r_pseudo(n,t) + x_diff1 - x_av_diff1
        x2r_pseudo(n,t) = x2r_pseudo(n,t) + x_diff2 - x_av_diff2
        y1r_pseudo(n,t) = y1r_pseudo(n,t) + y_diff1 - y_av_diff1
        y2r_pseudo(n,t) = y2r_pseudo(n,t) + y_diff2 - y_av_diff2
        
                    if(x1r_pseudo(n,t).lt.0)then
                x1r_pseudo(n,t) = dfloat(ii)+x1r_pseudo(n,t)
                x1r_pseudo_coord(n,t)=x1r_pseudo_coord(n,t)-1
            end if
            if(x1r_pseudo(n,t).gt.dfloat(ii))then
                x1r_pseudo(n,t) = x1r_pseudo(n,t) - dfloat(ii)
                x1r_pseudo_coord(n,t) = x1r_pseudo_coord(n,t)+1
            endif
            if(y1r_pseudo(n,t).lt.0.)then
              y1r_pseudo(n,t)=dfloat(jj)+y1r_pseudo(n,t)
              y1r_pseudo_coord(n,t) = y1r_pseudo_coord(n,t) -1
            endif
            if(y1r_pseudo(n,t).gt.dfloat(jj))then
              y1r_pseudo(n,t)=y1r_pseudo(n,t)-dfloat(jj)
              y1r_pseudo_coord(n,t)= y1r_pseudo_coord(n,t) +1
            endif
            
            if(x2r_pseudo(n,t).lt.0)then
                x2r_pseudo(n,t) = dfloat(ii)+x2r_pseudo(n,t)
                x2r_pseudo_coord(n,t)=x2r_pseudo_coord(n,t)-1
            end if
            if(x2r_pseudo(n,t).gt.dfloat(ii))then
                x2r_pseudo(n,t) = x2r_pseudo(n,t) - dfloat(ii)
                x2r_pseudo_coord(n,t) = x2r_pseudo_coord(n,t)+1
            endif
            if(y2r_pseudo(n,t).lt.0.)then
              y2r_pseudo(n,t)=dfloat(jj)+y2r_pseudo(n,t)
              y2r_pseudo_coord(n,t) = y2r_pseudo_coord(n,t) -1
            endif
            if(y2r_pseudo(n,t).gt.dfloat(jj))then
              y2r_pseudo(n,t)=y2r_pseudo(n,t)-dfloat(jj)
              y2r_pseudo_coord(n,t)= y2r_pseudo_coord(n,t) +1
            endif
        
        endif
            
            
            enddo
            endif
            
            
c ----------------- EDDY INDUCED PARTICLE ADVECTION ---------------
        if (i_eddy.eq.1) then

            do n = 1,npoints
            
        if (isolve.eq.0) then
        
        call rk4_bicubic(ii,jj,x1r_eddy(n,t),y1r_eddy(n,t),dt_nondim
     &   ,dfloat(0)
     & ,M1_eddy_old,M1_eddy_half,M1_eddy_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2r_eddy(n,t),y2r_eddy(n,t),dt_nondim
     &   ,dfloat(0)
     & ,M2_eddy_old,M2_eddy_half,M2_eddy_new
     & ,x_diff2,y_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1r_eddy(n,t),y1r_eddy(n,t),dt_nondim,
     &   dfloat(0)
     & ,a1_eddy_old,b1_eddy_old,c1_eddy_old
     & ,d1_eddy_old,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half
     & ,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new
     & ,x_diff1,y_diff1)
     
        !print*, 'x1r_eddy,y1r_eddy'
        !print*, x1r_eddy(n,t), y1r_eddy(n,t)
        call rk4_2dcubic(ii,jj,x2r_eddy(n,t),y2r_eddy(n,t),dt_nondim
     &   ,dfloat(0)
     & ,a2_eddy_old,b2_eddy_old,c2_eddy_old
     & ,d2_eddy_old,a2_eddy_half,b2_eddy_half,c2_eddy_half,d2_eddy_half
     & ,a2_eddy_new,b2_eddy_new,c2_eddy_new,d2_eddy_new
     & ,x_diff2,y_diff2)
        
        endif
            
        x1r_eddy(n,t) = x1r_eddy(n,t) + x_diff1
        y1r_eddy(n,t) = y1r_eddy(n,t) + y_diff1
        x2r_eddy(n,t) = x2r_eddy(n,t) + x_diff2
        y2r_eddy(n,t) = y2r_eddy(n,t) + y_diff2
        
            if(x1r_eddy(n,t).lt.0)then
                x1r_eddy(n,t) = dfloat(ii)+x1r_eddy(n,t)
                x1r_eddy_coord(n,t)=x1r_eddy_coord(n,t)-1
            end if
            if(x1r_eddy(n,t).gt.dfloat(ii))then
                x1r_eddy(n,t) = x1r_eddy(n,t) - dfloat(ii)
                x1r_eddy_coord(n,t) = x1r_eddy_coord(n,t)+1
            endif
            if(y1r_eddy(n,t).lt.0.)then
              y1r_eddy(n,t)=dfloat(jj)+y1r_eddy(n,t)
              y1r_eddy_coord(n,t) = y1r_eddy_coord(n,t) -1
            endif
            if(y1r_eddy(n,t).gt.dfloat(jj))then
              y1r_eddy(n,t)=y1r_eddy(n,t)-dfloat(jj)
              y1r_eddy_coord(n,t)= y1r_eddy_coord(n,t) +1
            endif
            
            if(x2r_eddy(n,t).lt.0)then
                x2r_eddy(n,t) = dfloat(ii)+x2r_eddy(n,t)
                x2r_eddy_coord(n,t)=x2r_eddy_coord(n,t)-1
            end if
            if(x2r_eddy(n,t).gt.dfloat(ii))then
                x2r_eddy(n,t) = x2r_eddy(n,t) - dfloat(ii)
                x2r_eddy_coord(n,t) = x2r_eddy_coord(n,t)+1
            endif
            if(y2r_eddy(n,t).lt.0.)then
              y2r_eddy(n,t)=dfloat(jj)+y2r_eddy(n,t)
              y2r_eddy_coord(n,t) = y2r_eddy_coord(n,t) -1
            endif
            if(y2r_eddy(n,t).gt.dfloat(jj))then
              y2r_eddy(n,t)=y2r_eddy(n,t)-dfloat(jj)
              y2r_eddy_coord(n,t)= y2r_eddy_coord(n,t) +1
            endif
            
           
            
            enddo
            endif
            

             
             if (k_s(t) .eq. l_tot_day*k_save) then
            print*, 'saving release number', t
            
            k_s(t) = 0
            
            print *, 'Writing to NETCDF files at time',time_day
            
     
            if (i_full.eq.1) then

            call write_release_file(full_name,npoints,release_no
     &       ,x1r(:,t),y1r(:,t)
     &       ,x2r(:,t),y2r(:,t)
     &       ,x1r_coord(:,t),y1r_coord(:,t),x2r_coord(:,t)
     &       ,y2r_coord(:,t)
     &       ,t,time_day,nrel(t))
     
            endif
     
            if (i_eddy.eq.1) then
            call write_release_file(eddy_name,npoints,release_no
     &       ,x1r_eddy(:,t)
     &       ,y1r_eddy(:,t)
     &       ,x2r_eddy(:,t),y2r_eddy(:,t)
     &       ,x1r_eddy_coord(:,t),y1r_eddy_coord(:,t)
     &       ,x2r_eddy_coord(:,t),y2r_eddy_coord(:,t)
     &       ,t,time_day,nrel(t))
            endif
            
            if (i_pseudo.eq.1) then
            call write_release_file(pseudo_name,npoints,release_no
     &       ,x1r_pseudo(:,t)
     &       ,y1r_pseudo(:,t)
     &       ,x2r_pseudo(:,t),y2r_pseudo(:,t)
     &       ,x1r_pseudo_coord(:,t),y1r_pseudo_coord(:,t)
     &       ,x2r_pseudo_coord(:,t)
     &       ,y2r_pseudo_coord(:,t),t,time_day,nrel(t))
            endif
            nrel(t) = nrel(t) + 1
            
            endif
            
            
            
            enddo
          
            
        enddo
      
      end program offline_transport
