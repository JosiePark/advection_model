	  �  8   k820309    4          12.1        �d\                                                                                                           
       mod_qg2_netcdf.f MOD_QG2_NETCDF                                                    
       #         @                                                    #ERRCODE                                                           #         @                                                    #CREATE_NETCDF_FILE%LEN    #CREATE_NETCDF_FILE%DFLOAT    #FILE_NAME    #NUM_X    #NUM_Y 	   #BASINSCALE 
                                                   LEN                                                 DFLOAT           D @                                                    1           D @                                                     D @                               	                                                      
     
       #         @                                                     #FILE_NAME    #PSI1    #PSI2    #NUM_X    #NUM_Y    #READ_TIM    #NEW_TIM    #STEP_TIM              D @                                                    1          D @                                                  
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               D @                                                  
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                                                                                                                                                                   
                 D                                     
                 D                                             #         @                                                  
   #FILE_NAME    #PSI1    #PSI2    #EPOT    #EKIN1    #EKIN2    #NUM_X    #NUM_Y    #SAVE_TIM    #STEP_TIM              D @                                                    1                                                              
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                                                                   
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                                                     
                                                      
                                                      
                                                                                                                               D @                                   
                 D @                                           #         @                                                    #CREATE_AVE_FILE%LEN     #CREATE_AVE_FILE%DFLOAT !   #FILE_NAME "   #II #   #JJ $   #BASINSCALE %                                                    LEN                                            !     DFLOAT           D @                               "                     1           D @                               #                      D @                               $                                                      %     
       #         @                                  &                   #FILE_NAME '   #II (   #JJ )   #PSI1_AV *   #PSI2_AV +   #TIME ,             D @                               '                     1                                            (                                                       )                                                     *                    
       p        5 � p        r (   p          5 � p        r (     5 � p        r )       5 � p        r (     5 � p        r )                                                              +                    
       p        5 � p        r (   p          5 � p        r (     5 � p        r )       5 � p        r (     5 � p        r )                               D @                              ,     
       #         @                                  -                   #FILE_NAME .   #II /   #JJ 0   #PSI1_AV 1   #PSI2_AV 2   #TIME_AV 3             D @                               .                     1                                            /                                                       0                     D @                              1                    
 #      p        5 � p        r /   p          5 � p        r /     5 � p        r 0       5 � p        r /     5 � p        r 0                              D @                              2                    
 $      p        5 � p        r /   p          5 � p        r /     5 � p        r 0       5 � p        r /     5 � p        r 0                               D @                              3     
       #         @                                  4                   #FILE_NAME 5   #TIME 6   #T_LEN 7             D @                               5                     1         D @                              6                   
 '              &                                                     D @                               7               �   (      fn#fn !   �   @   J   MOD_NETCDF_ERROR ,     U       HANDLE_ERR+MOD_NETCDF_ERROR 4   ]  @   a   HANDLE_ERR%ERRCODE+MOD_NETCDF_ERROR #   �  �       CREATE_NETCDF_FILE '   U  <      CREATE_NETCDF_FILE%LEN *   �  ?      CREATE_NETCDF_FILE%DFLOAT -   �  L   a   CREATE_NETCDF_FILE%FILE_NAME )     @   a   CREATE_NETCDF_FILE%NUM_X )   \  @   a   CREATE_NETCDF_FILE%NUM_Y .   �  @   a   CREATE_NETCDF_FILE%BASINSCALE    �  �       READ_NETCDF &   �  L   a   READ_NETCDF%FILE_NAME !   �  $  a   READ_NETCDF%PSI1 !   �  $  a   READ_NETCDF%PSI2 "     @   a   READ_NETCDF%NUM_X "   Z  @   a   READ_NETCDF%NUM_Y %   �  @   a   READ_NETCDF%READ_TIM $   �  @   a   READ_NETCDF%NEW_TIM %     @   a   READ_NETCDF%STEP_TIM    Z  �       WRITE_NETCDF '   	  L   a   WRITE_NETCDF%FILE_NAME "   c	  $  a   WRITE_NETCDF%PSI1 "   �
  $  a   WRITE_NETCDF%PSI2 "   �  @   a   WRITE_NETCDF%EPOT #   �  @   a   WRITE_NETCDF%EKIN1 #   +  @   a   WRITE_NETCDF%EKIN2 #   k  @   a   WRITE_NETCDF%NUM_X #   �  @   a   WRITE_NETCDF%NUM_Y &   �  @   a   WRITE_NETCDF%SAVE_TIM &   +  @   a   WRITE_NETCDF%STEP_TIM     k  �       CREATE_AVE_FILE $     <      CREATE_AVE_FILE%LEN '   S  ?      CREATE_AVE_FILE%DFLOAT *   �  L   a   CREATE_AVE_FILE%FILE_NAME #   �  @   a   CREATE_AVE_FILE%II #     @   a   CREATE_AVE_FILE%JJ +   ^  @   a   CREATE_AVE_FILE%BASINSCALE    �  �       WRITE_AVE_FILE )   )  L   a   WRITE_AVE_FILE%FILE_NAME "   u  @   a   WRITE_AVE_FILE%II "   �  @   a   WRITE_AVE_FILE%JJ '   �  $  a   WRITE_AVE_FILE%PSI1_AV '     $  a   WRITE_AVE_FILE%PSI2_AV $   =  @   a   WRITE_AVE_FILE%TIME    }  �       READ_AVE_FILE (     L   a   READ_AVE_FILE%FILE_NAME !   W  @   a   READ_AVE_FILE%II !   �  @   a   READ_AVE_FILE%JJ &   �  $  a   READ_AVE_FILE%PSI1_AV &   �  $  a   READ_AVE_FILE%PSI2_AV &     @   a   READ_AVE_FILE%TIME_AV    _  l       READ_TIME $   �  L   a   READ_TIME%FILE_NAME      �   a   READ_TIME%TIME     �  @   a   READ_TIME%T_LEN 