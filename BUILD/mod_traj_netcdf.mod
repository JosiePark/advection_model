	  XF     k820309    4          12.1        ùd\                                                                                                           
       mod_traj_netcdf.f MOD_TRAJ_NETCDF                                                    
       #         @                                                    #ERRCODE                                                           #         @                                                    #CREATE_TRAJ_FILE%LEN    #FILE_NAME    #NPOINTS                                                    LEN           D @                                                    1           D @                                           #         @                                                     #FILE_NAME 	   #NPOINTS 
   #X1    #Y1    #X2    #Y2    #X1_COORD    #Y1_COORD    #X2_COORD    #Y2_COORD    #SAVE_TIME    #STEP_TIME              D @                               	                     1                                            
                                                                         
     p          5  p        r 
       5  p        r 
                                                                                  
     p          5  p        r 
       5  p        r 
                                                                                  
     p          5  p        r 
       5  p        r 
                                                                                  
 	    p          5  p        r 
       5  p        r 
                                                                                        p          5  p        r 
       5  p        r 
                                                                                        p          5  p        r 
       5  p        r 
                                                                                        p          5  p        r 
       5  p        r 
                                                                                        p          5  p        r 
       5  p        r 
                               D @                                   
                 D @                                           #         @                                                     #FILE_NAME    #NPOINTS    #X1    #Y1    #X2    #Y2    #X1_COORD    #Y1_COORD    #X2_COORD    #Y2_COORD    #READ_TIM     #NEW_TIM !   #STEP_TIM "             D @                                                    1                                                                 D @                                                  
     p          5  p        r        5  p        r                               D @                                                  
     p          5  p        r        5  p        r                               D @                                                  
     p          5  p        r        5  p        r                               D @                                                  
     p          5  p        r        5  p        r                               D @                                                        p          5  p        r        5  p        r                               D @                                                        p          5  p        r        5  p        r                               D @                                                        p          5  p        r        5  p        r                               D @                                                        p          5  p        r        5  p        r                                                                      
                 D                                !     
                 D                                 "            #         @                                  #                  #CREATE_RELEASE_FILE%LEN $   #FILE_NAME %   #NPOINTS &   #RELEASE_NO '                                              $     LEN           D @                               %                     1           D @                               &                      D @                               '            #         @                                  (                   #FILE_NAME )   #NPOINTS *   #RELEASE_NO +   #X1 ,   #Y1 -   #X2 .   #Y2 /   #X1_COORD 0   #Y1_COORD 1   #X2_COORD 2   #Y2_COORD 3   #RELEASE 4   #SAVE_TIME 5   #STEP_TIME 6             D @                               )                     1                                            *                                                       +                                                     ,                    
 !    p          5  p        r *       5  p        r *                                                              -                    
 "    p          5  p        r *       5  p        r *                                                              .                    
 #    p          5  p        r *       5  p        r *                                                              /                    
 $    p          5  p        r *       5  p        r *                                                               0                         p          5  p        r *       5  p        r *                                                               1                         p          5  p        r *       5  p        r *                                                               2                         p          5  p        r *       5  p        r *                                                               3                          p          5  p        r *       5  p        r *                                                                4                      D @                              5     
                                                  6            #         @                                  7                   #FILE_NAME 8   #TIME 9   #T_LEN ;   #NREL :             D @                               8                     1          D @                              9                    
 -      p        5  p        r :   p          5  p        r :     5  p        r ;       5  p        r :     5  p        r ;                                                                ;                                                       :            #         @                                  <                   #FILE_NAME =   #T_LEN >   #NREL ?   #NPOINTS @             D @                               =                     1           D @                               >                      D @                               ?                      D @                               @            #         @                                  A                   #FILE_NAME B   #NPOINTS C   #X1 D   #Y1 E   #X1_COORD F   #Y1_COORD G   #REL_NO H   #T_INDEX I             
@ @                               B                    1           
                        @         C                    D @                              D                    
 0    p          5  p        r C       5  p        r C                              D @                              E                    
 1    p          5  p        r C       5  p        r C                              D @                               F                     .    p          5  p        r C       5  p        r C                              D @                               G                     /    p          5  p        r C       5  p        r C                               
                                  H                     
                                  I           #         @                                  J                  #CREATE_BINNED_FILE%LEN K   #FILE_NAME L   #NBINS M   #NPOINTS N   #RELEASE_NO O                                              K     LEN           
@ @                               L                    1           
@ @                               M                     
@ @                               N                     
@ @                               O           #         @                                  P                   #FILE_NAME Q   #NPOINTS R   #BIN_NO S   #RELEASE T   #X1 U   #Y1 V   #X2 W   #Y2 X   #X1_COORD Y   #Y1_COORD Z   #X2_COORD [   #Y2_COORD \   #SAVE_TIME ]   #STEP_TIME ^             
@ @                               Q                    1           
                                  R                     
                                  S                     
                                  T                    
                                 U                    
 :   p          5  p        r R       5  p        r R                              
                                 V                    
 ;   p          5  p        r R       5  p        r R                              
                                 W                    
 <   p          5  p        r R       5  p        r R                              
                                 X                    
 =   p          5  p        r R       5  p        r R                              
                                  Y                     >   p          5  p        r R       5  p        r R                              
                                  Z                     ?   p          5  p        r R       5  p        r R                              
                                  [                     @   p          5  p        r R       5  p        r R                              
                                  \                     A   p          5  p        r R       5  p        r R                               
@ @                              ]     
                
                                  ^           #         @                                  _                   #FILE_NAME `   #BIN_COUNT a   #NBINS b   #K c             D @                               `                     1          D @                               a                     K    p          5  p        r b       5  p        r b                                                                b                                                       c            #         @                                  d                   #FILE_NAME e   #REL_NO f   #NBINS g   #BIN_BOUNDARIES h             D @                               e                     1                                            f                                                       g                     D @                              h                    
 N    p          5  p        r g       5  p        r g                     #         @                                  i                   #FILE_NAME j   #NPOINTS k   #NBINS l   #X1 m   #Y1 n   #X2 o   #Y2 p   #X1_COORD q   #Y1_COORD r   #X2_COORD s   #Y2_COORD t   #START_INFO u             D @                               j                     1                                            k                                                       l                     D @                              m                    
 R      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                              D @                              n                    
 S      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                              D @                              o                    
 V      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                              D @                              p                    
 W      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                              D @                               q                     T      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                              D @                               r                     U      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                              D @                               s                     X      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                              D @                               t                     Y      p        5  p        r l   p          5  p        r l     5  p        r k       5  p        r l     5  p        r k                               D @                               u                    Q    p          p            p                          #         @                                  v                   #FILE_NAME w   #NREL x   #NBINS y   #T_LEN z   #TIME {             D @                               w                     1                                            x                                                       y                                                       z                     D @                              {                    
 \      p        5  p        r x   p          5  p        r x     5  p        r z       5  p        r x     5  p        r z                     #         @                                  |                   #FILE_NAME }   #NBINS ~   #NREL    #NPOINTS    #T_LEN              D @                               }                     1           D @                               ~                      D @                                                     D @                                                     D @                                           #         @                                                     #FILE_NAME    #NBINS    #BIN_WIDTH              D @                                                    1           D @                                                    D @                                                  
 _    p          5  p        r        5  p        r                      #         @                                                     #FILE_NAME    #NBINS    #NREL    #BIN_WIDTH              D @                                                    1                                                                                                                        D @                                                  
 `      p        5  p        r    p          5  p        r      5  p        r        5  p        r      5  p        r                             *      fn#fn !   Ê   @   J   MOD_NETCDF_ERROR ,   
  U       HANDLE_ERR+MOD_NETCDF_ERROR 4   _  @   a   HANDLE_ERR%ERRCODE+MOD_NETCDF_ERROR !     ~       CREATE_TRAJ_FILE %     <      CREATE_TRAJ_FILE%LEN +   Y  L   a   CREATE_TRAJ_FILE%FILE_NAME )   ¥  @   a   CREATE_TRAJ_FILE%NPOINTS     å  Ú       WRITE_TRAJ_FILE *   ¿  L   a   WRITE_TRAJ_FILE%FILE_NAME (     @   a   WRITE_TRAJ_FILE%NPOINTS #   K  ´   a   WRITE_TRAJ_FILE%X1 #   ÿ  ´   a   WRITE_TRAJ_FILE%Y1 #   ³  ´   a   WRITE_TRAJ_FILE%X2 #   g  ´   a   WRITE_TRAJ_FILE%Y2 )     ´   a   WRITE_TRAJ_FILE%X1_COORD )   Ï  ´   a   WRITE_TRAJ_FILE%Y1_COORD )     ´   a   WRITE_TRAJ_FILE%X2_COORD )   7	  ´   a   WRITE_TRAJ_FILE%Y2_COORD *   ë	  @   a   WRITE_TRAJ_FILE%SAVE_TIME *   +
  @   a   WRITE_TRAJ_FILE%STEP_TIME    k
  å       READ_TRAJ_FILE )   P  L   a   READ_TRAJ_FILE%FILE_NAME '     @   a   READ_TRAJ_FILE%NPOINTS "   Ü  ´   a   READ_TRAJ_FILE%X1 "     ´   a   READ_TRAJ_FILE%Y1 "   D  ´   a   READ_TRAJ_FILE%X2 "   ø  ´   a   READ_TRAJ_FILE%Y2 (   ¬  ´   a   READ_TRAJ_FILE%X1_COORD (   `  ´   a   READ_TRAJ_FILE%Y1_COORD (     ´   a   READ_TRAJ_FILE%X2_COORD (   È  ´   a   READ_TRAJ_FILE%Y2_COORD (   |  @   a   READ_TRAJ_FILE%READ_TIM '   ¼  @   a   READ_TRAJ_FILE%NEW_TIM (   ü  @   a   READ_TRAJ_FILE%STEP_TIM $   <         CREATE_RELEASE_FILE (   Í  <      CREATE_RELEASE_FILE%LEN .   	  L   a   CREATE_RELEASE_FILE%FILE_NAME ,   U  @   a   CREATE_RELEASE_FILE%NPOINTS /     @   a   CREATE_RELEASE_FILE%RELEASE_NO #   Õ  ÷       WRITE_RELEASE_FILE -   Ì  L   a   WRITE_RELEASE_FILE%FILE_NAME +     @   a   WRITE_RELEASE_FILE%NPOINTS .   X  @   a   WRITE_RELEASE_FILE%RELEASE_NO &     ´   a   WRITE_RELEASE_FILE%X1 &   L  ´   a   WRITE_RELEASE_FILE%Y1 &      ´   a   WRITE_RELEASE_FILE%X2 &   ´  ´   a   WRITE_RELEASE_FILE%Y2 ,   h  ´   a   WRITE_RELEASE_FILE%X1_COORD ,     ´   a   WRITE_RELEASE_FILE%Y1_COORD ,   Ð  ´   a   WRITE_RELEASE_FILE%X2_COORD ,     ´   a   WRITE_RELEASE_FILE%Y2_COORD +   8  @   a   WRITE_RELEASE_FILE%RELEASE -   x  @   a   WRITE_RELEASE_FILE%SAVE_TIME -   ¸  @   a   WRITE_RELEASE_FILE%STEP_TIME "   ø  v       READ_RELEASE_TIME ,   n  L   a   READ_RELEASE_TIME%FILE_NAME '   º  $  a   READ_RELEASE_TIME%TIME (   Þ  @   a   READ_RELEASE_TIME%T_LEN '     @   a   READ_RELEASE_TIME%NREL )   ^  y       READ_RELEASE_TIME_LENGTH 3   ×  L   a   READ_RELEASE_TIME_LENGTH%FILE_NAME /   #  @   a   READ_RELEASE_TIME_LENGTH%T_LEN .   c  @   a   READ_RELEASE_TIME_LENGTH%NREL 1   £  @   a   READ_RELEASE_TIME_LENGTH%NPOINTS "   ã  ©       READ_RELEASE_FILE ,      L   a   READ_RELEASE_FILE%FILE_NAME *   Ø   @   a   READ_RELEASE_FILE%NPOINTS %   !  ´   a   READ_RELEASE_FILE%X1 %   Ì!  ´   a   READ_RELEASE_FILE%Y1 +   "  ´   a   READ_RELEASE_FILE%X1_COORD +   4#  ´   a   READ_RELEASE_FILE%Y1_COORD )   è#  @   a   READ_RELEASE_FILE%REL_NO *   ($  @   a   READ_RELEASE_FILE%T_INDEX #   h$         CREATE_BINNED_FILE '   %  <      CREATE_BINNED_FILE%LEN -   ?%  L   a   CREATE_BINNED_FILE%FILE_NAME )   %  @   a   CREATE_BINNED_FILE%NBINS +   Ë%  @   a   CREATE_BINNED_FILE%NPOINTS .   &  @   a   CREATE_BINNED_FILE%RELEASE_NO "   K&  ó       WRITE_BINNED_FILE ,   >'  L   a   WRITE_BINNED_FILE%FILE_NAME *   '  @   a   WRITE_BINNED_FILE%NPOINTS )   Ê'  @   a   WRITE_BINNED_FILE%BIN_NO *   
(  @   a   WRITE_BINNED_FILE%RELEASE %   J(  ´   a   WRITE_BINNED_FILE%X1 %   þ(  ´   a   WRITE_BINNED_FILE%Y1 %   ²)  ´   a   WRITE_BINNED_FILE%X2 %   f*  ´   a   WRITE_BINNED_FILE%Y2 +   +  ´   a   WRITE_BINNED_FILE%X1_COORD +   Î+  ´   a   WRITE_BINNED_FILE%Y1_COORD +   ,  ´   a   WRITE_BINNED_FILE%X2_COORD +   6-  ´   a   WRITE_BINNED_FILE%Y2_COORD ,   ê-  @   a   WRITE_BINNED_FILE%SAVE_TIME ,   *.  @   a   WRITE_BINNED_FILE%STEP_TIME    j.  x       WRITE_WIDTH &   â.  L   a   WRITE_WIDTH%FILE_NAME &   ./  ´   a   WRITE_WIDTH%BIN_COUNT "   â/  @   a   WRITE_WIDTH%NBINS    "0  @   a   WRITE_WIDTH%K %   b0         WRITE_BIN_BOUNDARIES /   ä0  L   a   WRITE_BIN_BOUNDARIES%FILE_NAME ,   01  @   a   WRITE_BIN_BOUNDARIES%REL_NO +   p1  @   a   WRITE_BIN_BOUNDARIES%NBINS 4   °1  ´   a   WRITE_BIN_BOUNDARIES%BIN_BOUNDARIES !   d2  ×       READ_BINNED_FILE +   ;3  L   a   READ_BINNED_FILE%FILE_NAME )   3  @   a   READ_BINNED_FILE%NPOINTS '   Ç3  @   a   READ_BINNED_FILE%NBINS $   4  $  a   READ_BINNED_FILE%X1 $   +5  $  a   READ_BINNED_FILE%Y1 $   O6  $  a   READ_BINNED_FILE%X2 $   s7  $  a   READ_BINNED_FILE%Y2 *   8  $  a   READ_BINNED_FILE%X1_COORD *   »9  $  a   READ_BINNED_FILE%Y1_COORD *   ß:  $  a   READ_BINNED_FILE%X2_COORD *   <  $  a   READ_BINNED_FILE%Y2_COORD ,   '=     a   READ_BINNED_FILE%START_INFO &   »=         READ_BINNED_FILE_TIME 0   <>  L   a   READ_BINNED_FILE_TIME%FILE_NAME +   >  @   a   READ_BINNED_FILE_TIME%NREL ,   È>  @   a   READ_BINNED_FILE_TIME%NBINS ,   ?  @   a   READ_BINNED_FILE_TIME%T_LEN +   H?  $  a   READ_BINNED_FILE_TIME%TIME ,   l@         READ_BINNED_FILE_DIMENSIONS 6   ð@  L   a   READ_BINNED_FILE_DIMENSIONS%FILE_NAME 2   <A  @   a   READ_BINNED_FILE_DIMENSIONS%NBINS 1   |A  @   a   READ_BINNED_FILE_DIMENSIONS%NREL 4   ¼A  @   a   READ_BINNED_FILE_DIMENSIONS%NPOINTS 2   üA  @   a   READ_BINNED_FILE_DIMENSIONS%T_LEN    <B  q       READ_BIN_WIDTH )   ­B  L   a   READ_BIN_WIDTH%FILE_NAME %   ùB  @   a   READ_BIN_WIDTH%NBINS )   9C  ´   a   READ_BIN_WIDTH%BIN_WIDTH '   íC  {       READ_RELEASE_BIN_WIDTH 1   hD  L   a   READ_RELEASE_BIN_WIDTH%FILE_NAME -   ´D  @   a   READ_RELEASE_BIN_WIDTH%NBINS ,   ôD  @   a   READ_RELEASE_BIN_WIDTH%NREL 1   4E  $  a   READ_RELEASE_BIN_WIDTH%BIN_WIDTH 