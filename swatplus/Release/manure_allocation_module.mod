  Û7  q   k820309    ?          19.1        ¶'óf                                                                                                          
       D:\Projects\Models\swatp_paddy\src\manure_allocation_module.f90 MANURE_ALLOCATION_MODULE                                                 
                                                              #MALLOUT_ADD                                                           #MALLO_DIV_CONST                  @                               '$            #FERTNM    #FMINN    #FMINP    #FORGN    #FORGP 	   #FNH3N 
                                                                                                         C                                                                                 	                             	                         0.                                                 	                             	                         0.                                                 	                             	                         0.                                       	          	                             	                         0.                                       
           	                             	                         0.              @                              '            #MALLO_OBJ    #SRC_OBJ    #APP_T_HA    #APP_METHOD                                                                                                            0                                                                                                       0                                                 	                             	                         0.                                                                                                       0                                                #MANURE_DEMAND_AMOUNT                  @                              '            #STOR    #PROD    #WITHDR                                                      	                             	                         0.                                                 	                             	                         0.                                                 	                             	                         0.                                                #SOURCE_MANURE_OUTPUT                  @                              '            #NUM    #MOIS_TYP    #MANURE_TYP    #LAT    #LONG    #STOR_INIT    #STOR_MAX    #PROD_MON    #FERTDB    #BAL_D     #BAL_M !   #BAL_Y "   #BAL_A #                                                                                                                                                                                                                        	                                                $      	                                                (      	                                                ,      	                                                   0         	  p      p        p                                                           `   	                                                       d   
   #SOURCE_MANURE_OUTPUT                                            !        p      #SOURCE_MANURE_OUTPUT                                            "        |      #SOURCE_MANURE_OUTPUT                                            #              #SOURCE_MANURE_OUTPUT                  @              A           $     'Ô            #NUM %   #OB_TYP &   #OB_NUM '   #DTBL (   #RIGHT )   #DTBL_NUM *   #MANURE_AMT +   #WITHDR ,   #WITHDR_M -   #WITHDR_Y .   #WITHDR_A /                                           %                                                      &     
                                                 '                                                     (                                                     )        -                                              *     0                                                 +        4      #MANURE_DEMAND_AMOUNT                                          ,        D         	        &                                                             -        h      	   	        &                                                             .              
   	        &                                                             /        °         	        &                                     @               Á           0     '            #NAME 1   #RULE_TYP 2   #SRC_OBS 3   #DMD_OBS 4   #TOT 5   #SRC 6   #DMD 7                                          1                                                      2                                                      3     4                                                 4     8                                                 5        <      #SOURCE_MANURE_OUTPUT                                          6        H            #MANURE_SOURCE_OBJECTS          &                                                             7        l   Ô         #MANURE_DEMAND_OBJECTS $         &                                                              8                    &                       #MANURE_ALLOCATION 0                 @                          9     '           #DAY :   #MO ;   #DAY_MO <   #YRC =   #IDMD >   #DMD_TYP ?   #DMD_NUM @   #SRC1_OBJ A   #SRC1_TYP B   #SRC1_NUM C   #DMD1 D   #S1OUT E   #S1UN F   #SRC2_TYP G   #SRC2_NUM H   #DMD2 I   #S2OUT J   #S2UN K   #SRC3_TYP L   #SRC3_NUM M   #DMD3 N   #S3OUT O   #S3UN P                                         :                                                                  C  jday                                                      ;                                                                 C	 mon                                                       <                                                                 C day                                                        =                                                                 C yr                                                         >                                              	                   C unit                                                         ?                                                                  Cdmd_typ                                                               @        0                                                         C    dmd_num                                                           A        @                                                         C   src1_obj                                                       B        L   	                                                      C src1_typ                                                         C        X   
                                                      C src1_num	                                                        D        d                                                         C	  demand                                                            E        s                                                         Csrc1_withdraw                                                        F                                                                 C  src1_unmet                                                      G                                                                 C src2_typ                                                         H                                                                 C src2_num	                                                        I        ¦                                                         C	  demand                                                            J        µ                                                         Csrc2_withdraw                                                        K        Ä                                                         C  src2_unmet                                                      L        Ð                                                         C src3_typ                                                         M        Ü                                                         C src3_num	                                                        N        è                                                         C	  demand                                                            O        ÷                                                         Csrc3_withdraw                                                        P                                                                C  src3_unmet                                                        Q       #MALLO_HEADER 9                 @                          R     '           #DAY S   #MO T   #DAY_MO U   #YRC V   #IDMD W   #DMD_TYP X   #DMD_NUM Y   #SRC1_OBJ Z   #SRC1_TYP [   #SRC1_NUM \   #DMD1 ]   #S1OUT ^   #S1UN _   #SRC2_TYP `   #SRC2_NUM a   #DMD2 b   #S2OUT c   #S2UN d   #SRC3_TYP e   #SRC3_NUM f   #DMD3 g   #S3OUT h   #S3UN i                                         S                                               	                   C	                                                             T                                              	                   C	                                                             U                                              	                   C	                                                             V                                              	                   C	                                                             W                                               	                   C	                                                             X        (                                                         C	                                                                     Y        8                                                         C                                                                      Z        H                                                         C                                                                  [        T   	                                                      C	                                                                 \        `   
                                   	                   C                                                              ]        h                                                         Cm^3                                                                  ^        w                                                         Cm^3                                                                  _     	                                         
                   Cm^3                                                            `                                                                 C                                                                     a                                                                 C                                                                     b        ­                                                         Cm^3                                                                  c        ¼                                                         Cm^3                                                                  d     
   Ë                                                         Cm^3                                                             e        Õ                                                         C                                                                     f        ä                                                         C                                                                     g        ó                                                         Cm^3                                                                  h                                                                Cm^3                                                                  i     
                                                           Cm^3                                                               j       #MALLO_HEADER_UNITS R   &     @     X                                                   #MALLO1 k   #MALLO2 l   #SOURCE_MANURE_OUTPUT          
                                  k           #SOURCE_MANURE_OUTPUT          
                                  l           #SOURCE_MANURE_OUTPUT    &     @     X                                                   #MALLO1 m   #CONST n   #SOURCE_MANURE_OUTPUT          
                                  m           #SOURCE_MANURE_OUTPUT          
                                  n     	         a      fn#fn '     <   J   FERTILIZER_DATA_MODULE    =  M      i@      Q      i@ 5   Û         FERTILIZER_DB+FERTILIZER_DATA_MODULE <   b     a   FERTILIZER_DB%FERTNM+FERTILIZER_DATA_MODULE ;   ÿ  ~   a   FERTILIZER_DB%FMINN+FERTILIZER_DATA_MODULE ;   }  ~   a   FERTILIZER_DB%FMINP+FERTILIZER_DATA_MODULE ;   û  ~   a   FERTILIZER_DB%FORGN+FERTILIZER_DATA_MODULE ;   y  ~   a   FERTILIZER_DB%FORGP+FERTILIZER_DATA_MODULE ;   ÷  ~   a   FERTILIZER_DB%FNH3N+FERTILIZER_DATA_MODULE %   u  ~       MANURE_DEMAND_AMOUNT /   ó  }   a   MANURE_DEMAND_AMOUNT%MALLO_OBJ -   p  }   a   MANURE_DEMAND_AMOUNT%SRC_OBJ .   í  ~   a   MANURE_DEMAND_AMOUNT%APP_T_HA 0   k  }   a   MANURE_DEMAND_AMOUNT%APP_METHOD    è  R       MANURE_AMTZ %   :  d       SOURCE_MANURE_OUTPUT *     ~   a   SOURCE_MANURE_OUTPUT%STOR *   	  ~   a   SOURCE_MANURE_OUTPUT%PROD ,   	  ~   a   SOURCE_MANURE_OUTPUT%WITHDR    
  R       MALLOZ &   j
  á       MANURE_SOURCE_OBJECTS *   K  @   a   MANURE_SOURCE_OBJECTS%NUM /     @   a   MANURE_SOURCE_OBJECTS%MOIS_TYP 1   Ë  @   a   MANURE_SOURCE_OBJECTS%MANURE_TYP *     @   a   MANURE_SOURCE_OBJECTS%LAT +   K  @   a   MANURE_SOURCE_OBJECTS%LONG 0     @   a   MANURE_SOURCE_OBJECTS%STOR_INIT /   Ë  @   a   MANURE_SOURCE_OBJECTS%STOR_MAX /     x   a   MANURE_SOURCE_OBJECTS%PROD_MON -     @   a   MANURE_SOURCE_OBJECTS%FERTDB ,   Ã  Z   a   MANURE_SOURCE_OBJECTS%BAL_D ,     Z   a   MANURE_SOURCE_OBJECTS%BAL_M ,   w  Z   a   MANURE_SOURCE_OBJECTS%BAL_Y ,   Ñ  Z   a   MANURE_SOURCE_OBJECTS%BAL_A &   +  Î       MANURE_DEMAND_OBJECTS *   ù  @   a   MANURE_DEMAND_OBJECTS%NUM -   9  @   a   MANURE_DEMAND_OBJECTS%OB_TYP -   y  @   a   MANURE_DEMAND_OBJECTS%OB_NUM +   ¹  @   a   MANURE_DEMAND_OBJECTS%DTBL ,   ù  @   a   MANURE_DEMAND_OBJECTS%RIGHT /   9  @   a   MANURE_DEMAND_OBJECTS%DTBL_NUM 1   y  Z   a   MANURE_DEMAND_OBJECTS%MANURE_AMT -   Ó  l   a   MANURE_DEMAND_OBJECTS%WITHDR /   ?  l   a   MANURE_DEMAND_OBJECTS%WITHDR_M /   «  l   a   MANURE_DEMAND_OBJECTS%WITHDR_Y /     l   a   MANURE_DEMAND_OBJECTS%WITHDR_A "            MANURE_ALLOCATION '     @   a   MANURE_ALLOCATION%NAME +   T  @   a   MANURE_ALLOCATION%RULE_TYP *     @   a   MANURE_ALLOCATION%SRC_OBS *   Ô  @   a   MANURE_ALLOCATION%DMD_OBS &     Z   a   MANURE_ALLOCATION%TOT &   n     a   MANURE_ALLOCATION%SRC &   õ     a   MANURE_ALLOCATION%DMD    |  {       MALLO    ÷  M      MALLO_HEADER !   D     a   MALLO_HEADER%DAY     ×     a   MALLO_HEADER%MO $   j     a   MALLO_HEADER%DAY_MO !   ý     a   MALLO_HEADER%YRC "        a   MALLO_HEADER%IDMD %   %     a   MALLO_HEADER%DMD_TYP %   Â     a   MALLO_HEADER%DMD_NUM &   _     a   MALLO_HEADER%SRC1_OBJ &   ø     a   MALLO_HEADER%SRC1_TYP &        a   MALLO_HEADER%SRC1_NUM "   *     a   MALLO_HEADER%DMD1 #   Æ     a   MALLO_HEADER%S1OUT "   b     a   MALLO_HEADER%S1UN &   û     a   MALLO_HEADER%SRC2_TYP &         a   MALLO_HEADER%SRC2_NUM "   -!     a   MALLO_HEADER%DMD2 #   É!     a   MALLO_HEADER%S2OUT "   e"     a   MALLO_HEADER%S2UN &   þ"     a   MALLO_HEADER%SRC3_TYP &   #     a   MALLO_HEADER%SRC3_NUM "   0$     a   MALLO_HEADER%DMD3 #   Ì$     a   MALLO_HEADER%S3OUT "   h%     a   MALLO_HEADER%S3UN    &  J       MALLO_HDR #   K&  M      MALLO_HEADER_UNITS '   '     a   MALLO_HEADER_UNITS%DAY &   -(     a   MALLO_HEADER_UNITS%MO *   Â(     a   MALLO_HEADER_UNITS%DAY_MO '   W)     a   MALLO_HEADER_UNITS%YRC (   ì)     a   MALLO_HEADER_UNITS%IDMD +   *     a   MALLO_HEADER_UNITS%DMD_TYP +   +     a   MALLO_HEADER_UNITS%DMD_NUM ,   »+     a   MALLO_HEADER_UNITS%SRC1_OBJ ,   T,     a   MALLO_HEADER_UNITS%SRC1_TYP ,   í,     a   MALLO_HEADER_UNITS%SRC1_NUM (   -     a   MALLO_HEADER_UNITS%DMD1 )   .     a   MALLO_HEADER_UNITS%S1OUT (   º.     a   MALLO_HEADER_UNITS%S1UN ,   P/     a   MALLO_HEADER_UNITS%SRC2_TYP ,   ì/     a   MALLO_HEADER_UNITS%SRC2_NUM (   0     a   MALLO_HEADER_UNITS%DMD2 )   $1     a   MALLO_HEADER_UNITS%S2OUT (   À1     a   MALLO_HEADER_UNITS%S2UN ,   W2     a   MALLO_HEADER_UNITS%SRC3_TYP ,   ó2     a   MALLO_HEADER_UNITS%SRC3_NUM (   3     a   MALLO_HEADER_UNITS%DMD3 )   +4     a   MALLO_HEADER_UNITS%S3OUT (   Ç4     a   MALLO_HEADER_UNITS%S3UN     ^5  P       MALLO_HDR_UNITS    ®5  z       MALLOUT_ADD #   (6  V   a   MALLOUT_ADD%MALLO1 #   ~6  V   a   MALLOUT_ADD%MALLO2     Ô6  y       MALLO_DIV_CONST '   M7  V   a   MALLO_DIV_CONST%MALLO1 &   £7  8   a   MALLO_DIV_CONST%CONST 