  �  ,   k820309    ?          19.1        ,�f                                                                                                          
       D:\Projects\Models\swatp_paddy\src\ch_pathogen_module.f90 CH_PATHOGEN_MODULE                                                        #CHPATH_ADD                                                           #CHPATH_DIV                                                           #CHPATH_MULT                  @                               '            #PTH_IN    #PTH_OUT    #DRE_OFF    #REGROW    #WATER 	   #BENTHIC 
           �                                          	                             	                         0.        �                                         	                             	                         0.        �                                         	                             	                         0.        �                                         	                             	                         0.        �                               	          	                             	                         0.        �                               
          	                             	                         0.                                                          &                       #CH_PATHOGEN_OUTPUT                                                              &                       #CH_PATHOGEN_OUTPUT                                                              &                       #CH_PATHOGEN_OUTPUT                                                              &                       #CH_PATHOGEN_OUTPUT                                                    #CH_PATHOGEN_OUTPUT                                                    #CH_PATHOGEN_OUTPUT                                                    #CH_PATHOGEN_OUTPUT                                                    #CH_PATHOGEN_OUTPUT                                                    #CH_PATHOGEN_OUTPUT                  @                               '�            #DAY    #MO    #DAY_MO    #YRC    #ISD    #ID    #NAME    #PTH_IN    #PTH_OUT    #DRE_OFF    #REGROW    #WATER     #BENTHIC !           �                                                                                                C  jday                        �                                                                                               C   mon                        �                                                                                               C   day                        �                                                                                               C    yr                        �                                                                            	                   C   unit                         �                                                                             	                   C gis_id                         �                                      (                                                         C name                                   �                                      8                                                         C path_in_#cfu/10                        �                                      H   	                                                      C path_out_#cfu/                        �                                      W   
                                                      C  dre_off_#cfu/                        �                                      f                                                         C    regrow_#cfu                        �                                       u                                                         C water_stor_#cf                        �                              !        �                                                         C   benthic_#cfu                                                        "     �   #CH_PATHOGEN_HEADER    &     @     X                                                   #CHO1 #   #CHO2 $   #CH_PATHOGEN_OUTPUT          
                                  #           #CH_PATHOGEN_OUTPUT          
                                  $           #CH_PATHOGEN_OUTPUT    &     @     X                                                   #CH1 %   #CONST &   #CH_PATHOGEN_OUTPUT          
                                  %           #CH_PATHOGEN_OUTPUT          
                                  &     	  &     @     X                                                   #CONST '   #CHN1 (   #CH_PATHOGEN_OUTPUT          
                                  '     	        
                                  (           #CH_PATHOGEN_OUTPUT       �   U      fn#fn    �   L      i@    A  L      i@    �  M      i@ #   �  �       CH_PATHOGEN_OUTPUT *   h  ~   a   CH_PATHOGEN_OUTPUT%PTH_IN +   �  ~   a   CH_PATHOGEN_OUTPUT%PTH_OUT +   d  ~   a   CH_PATHOGEN_OUTPUT%DRE_OFF *   �  ~   a   CH_PATHOGEN_OUTPUT%REGROW )   `  ~   a   CH_PATHOGEN_OUTPUT%WATER +   �  ~   a   CH_PATHOGEN_OUTPUT%BENTHIC    \  |       CHPTH_D    �  |       CHPTH_M    T  |       CHPTH_Y    �  |       CHPTH_A    L  P       BCHPTH_D    �  P       BCHPTH_M    �  P       BCHPTH_Y    <  P       BCHPTH_A    �  P       CHPTHZ #   �  �       CH_PATHOGEN_HEADER '   �	  �   a   CH_PATHOGEN_HEADER%DAY &   >
  �   a   CH_PATHOGEN_HEADER%MO *   �
  �   a   CH_PATHOGEN_HEADER%DAY_MO '   d  �   a   CH_PATHOGEN_HEADER%YRC '   �  �   a   CH_PATHOGEN_HEADER%ISD &   �  �   a   CH_PATHOGEN_HEADER%ID (   !  �   a   CH_PATHOGEN_HEADER%NAME *   �  �   a   CH_PATHOGEN_HEADER%PTH_IN +   [  �   a   CH_PATHOGEN_HEADER%PTH_OUT +   �  �   a   CH_PATHOGEN_HEADER%DRE_OFF *   �  �   a   CH_PATHOGEN_HEADER%REGROW )   /  �   a   CH_PATHOGEN_HEADER%WATER +   �  �   a   CH_PATHOGEN_HEADER%BENTHIC    g  P       CHPATH_HDR    �  t       CHPATH_ADD     +  T   a   CHPATH_ADD%CHO1       T   a   CHPATH_ADD%CHO2    �  t       CHPATH_DIV    G  T   a   CHPATH_DIV%CH1 !   �  8   a   CHPATH_DIV%CONST    �  u       CHPATH_MULT "   H  8   a   CHPATH_MULT%CONST !   �  T   a   CHPATH_MULT%CHN1 