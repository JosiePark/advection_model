	  ź	     k820309    4          12.1        ső[                                                                                                           
       mod_diffusivity_functions.f MOD_DIFFUSIVITY_FUNCTIONS                                                    
       #         @                                                    #X_C    #NX    #F_C    #X    #F                                                                 
     p          5 O p            5 O p                                                                                                                                              
     p          5 O p            5 O p                                                                         
                                                      
       #         @                                                   #DIFFUSIVITY_INTERP_1D%DFLOAT 	   #XC 
   #NBINS    #KC    #II    #X    #K                                               	     DFLOAT                                          
                    
     p          5  p        r        5  p        r                                                                                                                                          
     p          5  p        r        5  p        r                                 @                                                     D @                                   
                 D @                                   
       #         @                                                     #XC    #NBINS    #KC    #II    #X    #DKDX             D @                                                  
     p          5  p        r        5  p        r                                D @                                                    D @                                                  
     p          5  p        r        5  p        r                                D @                                                     D @                                   
                 D                                     
              >      fn#fn    Ţ   @   J   MOD_1DINTERP '     p       INTERP_1D+MOD_1DINTERP +     ¤   a   INTERP_1D%X_C+MOD_1DINTERP *   2  @   a   INTERP_1D%NX+MOD_1DINTERP +   r  ¤   a   INTERP_1D%F_C+MOD_1DINTERP )     @   a   INTERP_1D%X+MOD_1DINTERP )   V  @   a   INTERP_1D%F+MOD_1DINTERP &            DIFFUSIVITY_INTERP_1D -   1  ?      DIFFUSIVITY_INTERP_1D%DFLOAT )   p  ´   a   DIFFUSIVITY_INTERP_1D%XC ,   $  @   a   DIFFUSIVITY_INTERP_1D%NBINS )   d  ´   a   DIFFUSIVITY_INTERP_1D%KC )     @   a   DIFFUSIVITY_INTERP_1D%II (   X  @   a   DIFFUSIVITY_INTERP_1D%X (     @   a   DIFFUSIVITY_INTERP_1D%K '   Ř  |       DIFFUSIVITY_DERIVATIVE *   T  ´   a   DIFFUSIVITY_DERIVATIVE%XC -     @   a   DIFFUSIVITY_DERIVATIVE%NBINS *   H  ´   a   DIFFUSIVITY_DERIVATIVE%KC *   ü  @   a   DIFFUSIVITY_DERIVATIVE%II )   <	  @   a   DIFFUSIVITY_DERIVATIVE%X ,   |	  @   a   DIFFUSIVITY_DERIVATIVE%DKDX 