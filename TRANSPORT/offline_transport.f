c PROGRAM THAT PERFORMS OFFLINE ADVECTION

      program offline_transport
      
      use mod_lag_input
      use mod_qg_input
      use mod_lagr_constants
      use mod_time_array
      use mod_qg2_netcdf
      use mod_read_psi_av
      use mod_traj_netcdf
      use mod_2dcubic
      use mod_bicubic
      use mod_random
      use mod_time_interp
      use mod_laplace
      use mod_rk4
      use mod_name_files
      
      implicit none 
      

      k_t = 1
      
      
c ------ CREATE TRAJECTORY DATA FILE ----------------

        

      if(i_full.eq.1.) then
      call name_files(home_dir,regime,1,file_name,full_name
     &       ,disp_file,ave_name)
      call create_traj_file(full_name,npoints)
      endif
      if (i_pseudo.eq.1) then
      call name_files(home_dir,regime,2,file_name,pseudo_name
     &       ,disp_file,ave_name)
      call create_traj_file(pseudo_name,npoints)
      endif
      if (i_eddy.eq.1) then
      call name_files(home_dir,regime,3,file_name,eddy_name
     &       ,disp_file,ave_name)
      call create_traj_file(eddy_name,npoints)
      endif
      
      
c --- READ TIME AND TIME AVERAGE STREAM FUNCTION ----
      
      call read_time(file_name,time,t_len)

        scale = basinscale/dfloat(ii)
        uscale = 1
        tscale = scale/uscale
        
      if (i_pseudo.eq.1) then
      
      call psi_av(ave_name,ii,jj,psi1_av,psi2_av)
      
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
      
      elseif (i_eddy.eq.1) then
      
      call psi_av(ave_name,ii,jj,psi1_av,psi2_av)
      
      endif
      

      
c ----- CALCULATE THE TRUE TIME ARRAY (WHICH TAKES INTO ACCOUNT THE CHOSEN OFFLINE TIME STEP)

        max_time_lag = time(t_len)*86400 ! largest time recording in seconds
        
        print*, 'max_time =', max_time_lag
        print*, 'time(t_len)= ',time(t_len)
    
        
        t_tot = int(dfloat(max_time_lag)/dt) ! total number of time steps
        l_tot_day = int(86400./dt) ! total number of time steps in a day
        
        print*, 't_tot=',t_tot
        
        !stop
        
        dt_nondim = dt/tscale
        
        allocate(time_o(t_tot))
        allocate(time_dim(t_tot))
        
        do k = 1,t_tot
            time_o(k) = k*dt_nondim ! non dimensional time
            time_dim(k) = time_o(k)*tscale/86400 ! dimensional time in days
        enddo
            
        dt05 = dt_nondim/2
        dt6 = dt_nondim/6

c ------------ LAGRANGIAN PARTICLES ------------------


        iseed = 102
        n=0
        nrec = 1
        do j = 1,npoints_sqrt
          do i = 1,npoints_sqrt
            n = n+1
            x0(n) = ran1(iseed)*dfloat(ii)
            y0(n) = ran1(iseed)*dfloat(jj)
            
          enddo
        enddo
        
        do n = 1,npoints
        if((i_full.eq.1).or.(i_pseudo.eq.1)) then
        x1(n) = x0(n)
        y1(n) = y0(n)
        x2(n) = x0(n)
        y2(n) = y0(n)
        x1_coord(n) = 0
        x2_coord(n) = 0
        y1_coord(n) = 0
        y2_coord(n) = 0 
        endif
        if (i_eddy.eq.1) then
        x1_eddy(n) = x0(n)
        y1_eddy(n) = y0(n)
        x2_eddy(n) = x0(n)
        y2_eddy(n) = y0(n)
        x1_eddy_coord(n) = 0
        y1_eddy_coord(n) = 0
        x2_eddy_coord(n) = 0
        y2_eddy_coord(n) = 0
        endif
        if (i_pseudo.eq.1) then
        x1_pseudo(n) = x0(n)
        y1_pseudo(n) = y0(n)
        x2_pseudo(n) = x0(n)
        y2_pseudo(n) = y0(n)
        x1_pseudo_coord(n) = 0
        y1_pseudo_coord(n) = 0
        x2_pseudo_coord(n) = 0
        y2_pseudo_coord(n) = 0

        endif
        enddo  
        write(*,*)'Lagrangian particles randomly generated'
        
        time_day = 0.
        
        if (i_full.eq.1) then

            call write_traj_file(full_name,npoints,x1,y1,x2,y2
     &       ,x1_coord,y1_coord,x2_coord,y2_coord,time_day,nrec)
     
            endif
     
            if (i_eddy.eq.1) then
            call write_traj_file(eddy_name,npoints,x1_eddy,y1_eddy
     &       ,x2_eddy,y2_eddy
     &       ,x1_eddy_coord,y1_eddy_coord,x2_eddy_coord,y2_eddy_coord
     &       ,time_day,nrec)
            endif
            
            if (i_pseudo.eq.1) then
            call write_traj_file(pseudo_name,npoints,x1_pseudo,y1_pseudo
     &       ,x2_pseudo,y2_pseudo
     &       ,x1_pseudo_coord,y1_pseudo_coord,x2_pseudo_coord
     &       ,y2_pseudo_coord,time_day,nrec)
            endif
            
            !endif
        

        
        
c ------------ MAIN CYCLE -----------------------------
        do k = 2,t_tot
        !print*, 'time = ',time(k)
        time_day = time_o(k)*tscale/86400.
        print*,'time_o(k)=',time_o(k)
        
        k_t = k_t + 1
        
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

     

        call read_netcdf(FILE_NAME,psi1_old,psi2_old,ii,jj,time(k_old),
     +  new_tim, step_tim)
            
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
            

            
            
c ------------------ RK4 - FULL ADVECTION ---------------------------            
            if ((i_full.eq.1).or.(i_pseudo.eq.1)) then
            do n = 1,npoints
            
        if (isolve.eq.0) then
        
        call rk4_bicubic(ii,jj,x1(n),y1(n),dt_nondim,U_0,M1_old,M1_half
     &   ,M1_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2(n),y2(n),dt_nondim,U_0,M2_old,M2_half
     &   ,M2_new
     & ,x_diff2,y_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1(n),y1(n),dt_nondim,U_0,a1_old,b1_old,
     &   c1_old
     & ,d1_old,a1_half,b1_half,c1_half,d1_half
     & ,a1_new,b1_new,c1_new,d1_new
     & ,x_diff1,y_diff1)
        call rk4_2dcubic(ii,jj,x2(n),y2(n),dt_nondim,U_0,a2_old,b2_old
     &   ,c2_old
     & ,d2_old,a2_half,b2_half,c2_half,d2_half
     & ,a2_new,b2_new,c2_new,d2_new
     & ,x_diff2,y_diff2)
        
        endif
        
c MEAN CONTRIBUTION FOR THE PSEUDO TRAJECTORIES

        if(i_pseudo.eq.1) then
        
                 if (isolve.eq.0) then
         
        call rk4_bicubic(ii,jj,x1(n),y1(n),dt_nondim,U_0
     & ,M1_av,M1_av,M1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_bicubic(ii,jj,x2(n),y2(n),dt_nondim,U_0
     & ,M2_av,M2_av,M2_av
     & ,x_av_diff2,y_av_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1(n),y1(n),dt_nondim,U_0
     & ,a1_av,b1_av,c1_av
     & ,d1_av,a1_av,b1_av,c1_av,d1_av
     & ,a1_av,b1_av,c1_av,d1_av
     & ,x_av_diff1,y_av_diff1)
        call rk4_2dcubic(ii,jj,x2(n),y2(n),dt_nondim,U_0
     & ,a2_av,b2_av,c2_av
     & ,d2_av,a2_av,b2_av,c2_av,d2_av
     & ,a2_av,b2_av,c2_av,d2_av
     & ,x_av_diff2,y_av_diff2)
        
        endif 
        
        ENDIF
     
        
        x1(n) = x1(n) + x_diff1
        y1(n) = y1(n) + y_diff1
        x2(n) = x2(n) + x_diff2
        y2(n) = y2(n) + y_diff2
        
            if(x1(n).lt.0)then
                x1(n) = dfloat(ii)+x1(n)
                x1_coord(n)=x1_coord(n)-1
            end if
            if(x1(n).gt.dfloat(ii))then
                x1(n) = x1(n) - dfloat(ii)
                x1_coord(n) = x1_coord(n)+1
            endif
            if(y1(n).lt.0.)then
              y1(n)=dfloat(jj)+y1(n)
              y1_coord(n) = y1_coord(n) -1
            endif
            if(y1(n).gt.dfloat(jj))then
              y1(n)=y1(n)-dfloat(jj)
              y1_coord(n)= y1_coord(n) +1
            endif
            
            if(x2(n).lt.0)then
                x2(n) = dfloat(ii)+x2(n)
                x2_coord(n)=x2_coord(n)-1
            end if
            if(x2(n).gt.dfloat(ii))then
                x2(n) = x2(n) - dfloat(ii)
                x2_coord(n) = x2_coord(n)+1
            endif
            if(y2(n).lt.0.)then
              y2(n)=dfloat(jj)+y2(n)
              y2_coord(n) = y2_coord(n) -1
            endif
            if(y2(n).gt.dfloat(jj))then
              y2(n)=y2(n)-dfloat(jj)
              y2_coord(n)= y2_coord(n) +1
            endif
            
            if(i_pseudo.eq.1) then
            
        x1_pseudo(n) = x1_pseudo(n) + x_diff1 - x_av_diff1
        x2_pseudo(n) = x2_pseudo(n) + x_diff2 - x_av_diff2
        y1_pseudo(n) = y1_pseudo(n) + y_diff1 - y_av_diff1
        y2_pseudo(n) = y2_pseudo(n) + y_diff2 - y_av_diff2
        
            if(x1_pseudo(n).lt.0)then
                x1_pseudo(n) = dfloat(ii)+x1_pseudo(n)
                x1_pseudo_coord(n)=x1_pseudo_coord(n)-1
            end if
            if(x1_pseudo(n).gt.dfloat(ii))then
                x1_pseudo(n) = x1_pseudo(n) - dfloat(ii)
                x1_pseudo_coord(n) = x1_pseudo_coord(n)+1
            endif
            if(y1_pseudo(n).lt.0.)then
              y1_pseudo(n)=dfloat(jj)+y1_pseudo(n)
              y1_pseudo_coord(n) = y1_pseudo_coord(n) -1
            endif
            if(y1_pseudo(n).gt.dfloat(jj))then
              y1_pseudo(n)=y1_pseudo(n)-dfloat(jj)
              y1_pseudo_coord(n)= y1_pseudo_coord(n) +1
            endif
            
            if(x2_pseudo(n).lt.0)then
                x2_pseudo(n) = dfloat(ii)+x2_pseudo(n)
                x2_pseudo_coord(n)=x2_pseudo_coord(n)-1
            end if
            if(x2_pseudo(n).gt.dfloat(ii))then
                x2_pseudo(n) = x2_pseudo(n) - dfloat(ii)
                x2_pseudo_coord(n) = x2_pseudo_coord(n)+1
            endif
            if(y2_pseudo(n).lt.0.)then
              y2_pseudo(n)=dfloat(jj)+y2_pseudo(n)
              y2_pseudo_coord(n) = y2_pseudo_coord(n) -1
            endif
            if(y2_pseudo(n).gt.dfloat(jj))then
              y2_pseudo(n)=y2_pseudo(n)-dfloat(jj)
              y2_pseudo_coord(n)= y2_pseudo_coord(n) +1
            endif
        
           endif
            
            
            enddo
            endif
            
            
c ----------------- EDDY INDUCED PARTICLE ADVECTION ---------------
        if (i_eddy.eq.1) then
            do n = 1,npoints
            
        if (isolve.eq.0) then
        
        call rk4_bicubic(ii,jj,x1_eddy(n),y1_eddy(n),dt_nondim,dfloat(0)
     & ,M1_eddy_old,M1_eddy_half,M1_eddy_new
     & ,x_diff1,y_diff1)
        call rk4_bicubic(ii,jj,x2_eddy(n),y2_eddy(n),dt_nondim,dfloat(0)
     & ,M2_eddy_old,M2_eddy_half,M2_eddy_new
     & ,x_diff2,y_diff2)
        
        elseif (isolve.eq.1) then
        
        
        call rk4_2dcubic(ii,jj,x1_eddy(n),y1_eddy(n),dt_nondim,dfloat(0)
     & ,a1_eddy_old,b1_eddy_old,c1_eddy_old
     & ,d1_eddy_old,a1_eddy_half,b1_eddy_half,c1_eddy_half,d1_eddy_half
     & ,a1_eddy_new,b1_eddy_new,c1_eddy_new,d1_eddy_new
     & ,x_diff1,y_diff1)
        call rk4_2dcubic(ii,jj,x2_eddy(n),y2_eddy(n),dt_nondim,dfloat(0)
     & ,a2_eddy_old,b2_eddy_old,c2_eddy_old
     & ,d2_eddy_old,a2_eddy_half,b2_eddy_half,c2_eddy_half,d2_eddy_half
     & ,a2_eddy_new,b2_eddy_new,c2_eddy_new,d2_eddy_new
     & ,x_diff2,y_diff2)
        
        endif
            
        x1_eddy(n) = x1_eddy(n) + x_diff1
        y1_eddy(n) = y1_eddy(n) + y_diff1
        x2_eddy(n) = x2_eddy(n) + x_diff2
        y2_eddy(n) = y2_eddy(n) + y_diff2
        
            if(x1_eddy(n).lt.0)then
                x1_eddy(n) = dfloat(ii)+x1_eddy(n)
                x1_eddy_coord(n)=x1_eddy_coord(n)-1
            end if
            if(x1_eddy(n).gt.dfloat(ii))then
                x1_eddy(n) = x1_eddy(n) - dfloat(ii)
                x1_eddy_coord(n) = x1_eddy_coord(n)+1
            endif
            if(y1_eddy(n).lt.0.)then
              y1_eddy(n)=dfloat(jj)+y1_eddy(n)
              y1_eddy_coord(n) = y1_eddy_coord(n) -1
            endif
            if(y1_eddy(n).gt.dfloat(jj))then
              y1_eddy(n)=y1_eddy(n)-dfloat(jj)
              y1_eddy_coord(n)= y1_eddy_coord(n) +1
            endif
            
            if(x2_eddy(n).lt.0)then
                x2_eddy(n) = dfloat(ii)+x2_eddy(n)
                x2_eddy_coord(n)=x2_eddy_coord(n)-1
            end if
            if(x2_eddy(n).gt.dfloat(ii))then
                x2_eddy(n) = x2_eddy(n) - dfloat(ii)
                x2_eddy_coord(n) = x2_eddy_coord(n)+1
            endif
            if(y2_eddy(n).lt.0.)then
              y2_eddy(n)=dfloat(jj)+y2_eddy(n)
              y2_eddy_coord(n) = y2_eddy_coord(n) -1
            endif
            if(y2_eddy(n).gt.dfloat(jj))then
              y2_eddy(n)=y2_eddy(n)-dfloat(jj)
              y2_eddy_coord(n)= y2_eddy_coord(n) +1
            endif
            
            enddo
            endif
            
            
             if (k_t .eq. int(l_tot_day*k_save)) then
            nrec= nrec + 1
            k_t = 0

           print*, 'time=',time_day,'days'
           
     
            if (i_full.eq.1) then

            call write_traj_file(full_name,npoints,x1,y1,x2,y2
     &       ,x1_coord,y1_coord,x2_coord,y2_coord,time_day,nrec)
     
            endif
     
            if (i_eddy.eq.1) then
            call write_traj_file(eddy_name,npoints,x1_eddy,y1_eddy
     &       ,x2_eddy,y2_eddy
     &       ,x1_eddy_coord,y1_eddy_coord,x2_eddy_coord,y2_eddy_coord
     &       ,time_day,nrec)
            endif
            
            if (i_pseudo.eq.1) then
            call write_traj_file(pseudo_name,npoints,x1_pseudo,y1_pseudo
     &       ,x2_pseudo,y2_pseudo
     &       ,x1_pseudo_coord,y1_pseudo_coord,x2_pseudo_coord
     &       ,y2_pseudo_coord,time_day,nrec)
            endif
            
            endif
          
            
        enddo
      
      end program offline_transport
