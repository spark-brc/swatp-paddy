  �;  R   k820309    ?          19.1        ���f                                                                                                          
       D:\Projects\Models\swatp_paddy\src\res_cs_module.f90 RES_CS_MODULE                   @                              '0                    #INFLOW    #OUTFLOW    #SEEP    #SETTLE    #RCTN    #PROD    #FERT    #IRRIG 	   #DIV 
   #MASS    #CONC    #VOLM                �                                              	                                                	                                 0.                �                                             	                                                	                                 0.                �                                             	                                                	                                 0.                �                                             	                                                	                                 0.                �                                             	                                                	                                 0.                �                                             	                                                	                                 0.                �                                             	                                                	                                 0.                �                               	              	                                                	                                 0.                �                               
             	  	                                                	                                 0.                �                                    $       
  	                                                	                                 0.                �                                    (         	                                                	                                 0.                �                                    ,         	                                                	                                 0.                      @               @                'H                    #CS               �                                                   0             #RES_CS_BALANCE              &                                                                                             0       #RES_CS_BALANCE                                                          H                        &                                           #RES_CS_OUTPUT                                                          H                        &                                           #RES_CS_OUTPUT                                                          H                        &                                           #RES_CS_OUTPUT                                                          H                        &                                           #RES_CS_OUTPUT                                                          H                        &                                           #RES_CS_OUTPUT                                                          H                        &                                           #RES_CS_OUTPUT                                                          H                        &                                           #RES_CS_OUTPUT                                                          H                        &                                           #RES_CS_OUTPUT                      @                               'L                    #NAME    #V_SEO4    #V_SEO3    #V_BORN    #K_SEO4    #K_SEO3    #K_BORN     #THETA_SEO4 !   #THETA_SEO3 "   #THETA_BORN #   #C_SEO4 $   #C_SEO3 %   #C_BORN &                �                                                                     �                                             	                                                	                 o�:            0.001                �                                              	                                                	                 o�:            0.001                �                                    $         	                                                	                 o�:            0.001                �                                    (         	                                                	                 ��L=            0.05                �                                    ,         	                                                	                 ���<            0.03                �                                     0         	                                                	                                 0.00                �                               !     4         	                                                	                 q=�?            1.08                �                               "     8       	  	                                                	                 q=�?            1.08                �                               #     <       
  	                                                	                 q=�?            1.08                �                               $     @         	                                                	                                 0.                �                               %     D         	                                                	                                 0.                �                               &     H         	                                                	                                 0.                                               '            L                        &                                           #RESERVOIR_CS_DATA                      @                          (     '*             (      #DAY )   #MO *   #DAY_MO +   #YRC ,   #ISD -   #ID .   #SEO4IN /   #SEO3IN 0   #BORNIN 1   #SEO4OUT 2   #SEO3OUT 3   #BORNOUT 4   #SEO4SEEP 5   #SEO3SEEP 6   #BORNSEEP 7   #SEO4SETL 8   #SEO3SETL 9   #BORNSETL :   #SEO4RCTN ;   #SEO3RCTN <   #BORNRCTN =   #SEO4PROD >   #SEO3PROD ?   #BORNPROD @   #SEO4FERT A   #SEO3FERT B   #BORNFERT C   #SEO4IRR D   #SEO3IRR E   #BORNIRR F   #SEO4DIV G   #SEO3DIV H   #BORNDIV I   #SEO4 J   #SEO3 K   #BORN L   #SEO4C M   #SEO3C N   #BORNC O   #VOLM P               �                              )                                                                                                          C  jday                                �                              *                                                                                                         C   mon                                �                              +                                                                                                         C   day                                �                              ,                                                                                                         C    yr                                �                              -                                                                                  	                       C   unit                                 �                              .                                                                                                          C gis_id                                     �                              /            ,                                                                                             Cseo4_in                                        �                              0            ;                                                                                             Cseo3_in                                        �                              1            J       	                                                                                      Cborn_in                                        �                              2            Y       
                                                                                      Cseo4_out                                       �                              3            h                                                                                             Cseo3_out                                       �                              4            w                                                                                             Cborn_out                                       �                              5            �                                                                                             Cseo4_seep                                      �                              6            �                                                                                             Cseo3_seep                                      �                              7            �                                                                                             Cborn_seep                                      �                              8            �                                                                                             Cseo4_setl                                      �                              9            �                                                                                             Cseo3_setl                                      �                              :            �                                                                                             Cborn_setl                                      �                              ;            �                                                                                             Cseo4_rctn                                      �                              <            �                                                                                             Cseo3_rctn                                      �                              =            �                                                                                             Cborn_rctn                                      �                              >                                                                                                        Cseo4_prod                                      �                              ?                                                                                                        Cseo3_prod                                      �                              @            +                                                                                            Cborn_prod                                      �                              A            :                                                                                            Cseo4_fert                                      �                              B            I                                                                                            Cseo3_fert                                      �                              C            X                                                                                            Cborn_fert                                      �                              D            g                                                                                            Cseo4_irrg                                      �                              E            v                                                                                            Cseo3_irrg                                      �                              F            �                                                                                            Cborn_irrg                                      �                              G            �                                                                                            Cseo4_div                                       �                              H            �                                                                                             Cseo3_div                                       �                              I            �      !                                                                                      Cborn_div                                       �                              J            �      "                                                                                      Cseo4_mass                                      �                              K            �      #                                                                                      Cseo3_mass                                      �                              L            �      $                                                                                      Cborn_mass                                      �                              M            �      %                                                                                      Cseo4_conc                                      �                              N            �      &                                                                                      Cseo3_conc                                      �                              O                  '                                                                                      Cborn_conc                                      �                              P                  (                                                                                      Cvol_m3                                                                          Q     *      #RES_CS_HEADER (      �   K      fn#fn    �   �       RES_CS_BALANCE &   �  �   a   RES_CS_BALANCE%INFLOW '   `  �   a   RES_CS_BALANCE%OUTFLOW $     �   a   RES_CS_BALANCE%SEEP &   �  �   a   RES_CS_BALANCE%SETTLE $   R  �   a   RES_CS_BALANCE%RCTN $   �  �   a   RES_CS_BALANCE%PROD $   �  �   a   RES_CS_BALANCE%FERT %   D  �   a   RES_CS_BALANCE%IRRIG #   �  �   a   RES_CS_BALANCE%DIV $   �  �   a   RES_CS_BALANCE%MASS $   6  �   a   RES_CS_BALANCE%CONC $   �  �   a   RES_CS_BALANCE%VOLM    �	  X       RES_CS_OUTPUT !   �	  �   a   RES_CS_OUTPUT%CS    �
  T       RES_CSBZ    �
  �       RESCS_D    u  �       RESCS_M      �       RESCS_Y    �  �       RESCS_A    R  �       WETCS_D    �  �       WETCS_M    �  �       WETCS_Y    /  �       WETCS_A "   �  �       RESERVOIR_CS_DATA '   �  P   a   RESERVOIR_CS_DATA%NAME )     �   a   RESERVOIR_CS_DATA%V_SEO4 )   �  �   a   RESERVOIR_CS_DATA%V_SEO3 )   f  �   a   RESERVOIR_CS_DATA%V_BORN )     �   a   RESERVOIR_CS_DATA%K_SEO4 )   �  �   a   RESERVOIR_CS_DATA%K_SEO3 )   _  �   a   RESERVOIR_CS_DATA%K_BORN -     �   a   RESERVOIR_CS_DATA%THETA_SEO4 -   �  �   a   RESERVOIR_CS_DATA%THETA_SEO3 -   W  �   a   RESERVOIR_CS_DATA%THETA_BORN )   �  �   a   RESERVOIR_CS_DATA%C_SEO4 )   �  �   a   RESERVOIR_CS_DATA%C_SEO3 )   K  �   a   RESERVOIR_CS_DATA%C_BORN    �  �       RES_CS_DATA    �  ;      RES_CS_HEADER "   �  �   a   RES_CS_HEADER%DAY !   �  �   a   RES_CS_HEADER%MO %   U  �   a   RES_CS_HEADER%DAY_MO "     �   a   RES_CS_HEADER%YRC "   �  �   a   RES_CS_HEADER%ISD !   �  �   a   RES_CS_HEADER%ID %   i   �   a   RES_CS_HEADER%SEO4IN %   5!  �   a   RES_CS_HEADER%SEO3IN %   "  �   a   RES_CS_HEADER%BORNIN &   �"  �   a   RES_CS_HEADER%SEO4OUT &   �#  �   a   RES_CS_HEADER%SEO3OUT &   e$  �   a   RES_CS_HEADER%BORNOUT '   1%  �   a   RES_CS_HEADER%SEO4SEEP '   �%  �   a   RES_CS_HEADER%SEO3SEEP '   �&  �   a   RES_CS_HEADER%BORNSEEP '   �'  �   a   RES_CS_HEADER%SEO4SETL '   a(  �   a   RES_CS_HEADER%SEO3SETL '   -)  �   a   RES_CS_HEADER%BORNSETL '   �)  �   a   RES_CS_HEADER%SEO4RCTN '   �*  �   a   RES_CS_HEADER%SEO3RCTN '   �+  �   a   RES_CS_HEADER%BORNRCTN '   ],  �   a   RES_CS_HEADER%SEO4PROD '   )-  �   a   RES_CS_HEADER%SEO3PROD '   �-  �   a   RES_CS_HEADER%BORNPROD '   �.  �   a   RES_CS_HEADER%SEO4FERT '   �/  �   a   RES_CS_HEADER%SEO3FERT '   Y0  �   a   RES_CS_HEADER%BORNFERT &   %1  �   a   RES_CS_HEADER%SEO4IRR &   �1  �   a   RES_CS_HEADER%SEO3IRR &   �2  �   a   RES_CS_HEADER%BORNIRR &   �3  �   a   RES_CS_HEADER%SEO4DIV &   U4  �   a   RES_CS_HEADER%SEO3DIV &   !5  �   a   RES_CS_HEADER%BORNDIV #   �5  �   a   RES_CS_HEADER%SEO4 #   �6  �   a   RES_CS_HEADER%SEO3 #   �7  �   a   RES_CS_HEADER%BORN $   Q8  �   a   RES_CS_HEADER%SEO4C $   9  �   a   RES_CS_HEADER%SEO3C $   �9  �   a   RES_CS_HEADER%BORNC #   �:  �   a   RES_CS_HEADER%VOLM    �;  S       RESCS_HDR 