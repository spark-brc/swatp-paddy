  ýn    k820309    ?          19.1        µ'óf                                                                                                          
       D:\Projects\Models\swatp_paddy\src\carbon_module.f90 CARBON_MODULE                                                        #CARBON_SOIL_FLUX__ADD    #CARBON_SOIL_GL__ADD    #CARBON_RESIDUE_GL__ADD    #CARBON_PLANT_GL__ADD                                                           #CARBON_SOIL_FLUX_MULT    #CARBON_SOIL_GL_MULT    #CARBON_RESIDUE_GL_MULT    #CARBON_PLANT_GL_MULT                                                           #CARBON_SOIL_FLUX_DIV 	   #CARBON_SOIL_GL_DIV 
   #CARBON_RESIDUE_GL_DIV    #CARBON_PLANT_GL_DIV                  @                               '      %      #ER_POC_PARA    #CFB_PARA    #SF_PARA_SUR    #SF_PARA_SUB    #ABL_PARA    #PEROC_DIC_PARA    #PEROC_DOC_PARA    #PART_DOC_PARA    #HLIFE_DOC_PARA    #ABCO2_PARA_SUR    #ABCO2_PARA_SUB    #ABP_PARA_SUR    #ABP_PARA_SUB    #ALMCO2_PARA_SUR    #ALMCO2_PARA_SUB    #ALSLNCO2_PARA_SUR    #ALSLNCO2_PARA_SUB    #ASP_PARA_SUR    #ASP_PARA_SUB     #ALSLCO2_PARA !   #APCO2_PARA "   #ASCO2_PARA #   #PRMT_51_PARA $   #PRMT_45_PARA %   #BMR_PARA_SUR &   #BMR_PARA_SUB '   #HPR_PARA (   #HSR_PARA )   #LMR_PARA_SUR *   #LMR_PARA_SUB +   #LSR_PARA_SUR ,   #LSR_PARA_SUB -   #XBM_PARA_SUR .   #XBM_PARA_SUB /   #XLSLF_PARA 0   #OX_AA_PARA 1   #OX_BB_PARA 2                                                     	                             	           À?            1.5                                                 	                             	         =
×>            0.42                                                 	                             	         ÍÌL=            0.05                                                 	                             	         ÍÌÌ=            0.10                                                 	                             	                         0.0                                                 	                             	         33s?            0.95                                                 	                             	         333?            0.70                                                 	                             	           zE            4000.                                                	  	                             	           HB            50.                                            $   
  	                             	         ?            0.6                                            (     	                             	                         0.                                            ,     	                             	                         0.0                                            0     	                             	                         0.0                                            4     	                             	         ?            0.6                                            8     	                             	         ÍÌ?            0.55                                            <     	                             	         ?            0.6                                            @     	                             	         ÍÌ?            0.55                                            D     	                             	                         0.0                                             H     	                             	                         0.0                                       !     L     	                             	         >            0.3                                       "     P     	                             	         ÍÌ?            0.55                                       #     T     	                             	         ÍÌ?            0.55                                       $     X     	                             	           ?            1.0                                       %     \     	                             	         ¦D;            0.003                                       &     `     	                             	         KY<            0.0164                                       '     d     	                             	         
×£<            0.02                                       (     h     	                             	         SI7            0.000012                                       )     l     	                             	         ¨§:            0.000548                                       *     p     	                             	         Tã%=            0.0405                                       +     t     	                             	         ÎªO=            0.0507                                       ,     x     	                             	         O/<            0.0107                                       -     |      	                             	         ÐDX<            0.0132                                       .        !  	                             	           ?            1.0                                       /        "  	                             	                         0.0                                       0        #  	                             	                         0.0                                       1        $  	                             	            A            10.0                                       2        %  	                             	         )\=            0.035                                        3        #CARBON_TERRESTRIAL_INPUTS                  @                          4     '(      
      #HP_RATE 5   #HS_RATE 6   #MICROB_RATE 7   #META_RATE 8   #STR_RATE 9   #MICROB_TOP_RATE :   #HS_HP ;   #MICROB_KOC <   #MIN_N_FRAC =   #C_ORG_FRAC >                                          5           	                             	                         0.                                       6          	                             	                         0.                                       7          	                             	                         0.                                       8          	                             	                         0.                                       9          	                             	                         0.                                       :          	                             	                         0.                                       ;          	                             	                         0.                                       <          	                             	                         0.                                       =         	  	                             	                         0.                                       >     $   
  	                             	                         0.                                        ?     (   #CARBON_INPUTS 4                                           @     (   #CARBON_INPUTS 4                 @                           A     '$      	      #ABCO2 B   #ABL C   #ABP D   #ALMCO2 E   #ALSLCO2 F   #ALSLNCO2 G   #APCO2 H   #ASCO2 I   #ASP J                                           B            	                                           C           	                                           D           	                                           E           	                                           F           	                                           G           	                                           H           	                                           I           	                                           J         	   	                                           K     $   #ORGANIC_ALLOCATIONS A                                           L     $   #ORGANIC_ALLOCATIONS A                 @                           M     '            #CDG N   #CS O   #OX P   #SUT Q   #X1 R   #XBMT S   #XLSLF T                                           N            	                                           O           	                                           P           	                                           Q           	                                           R           	                                           S           	                                           T           	                                           U        #ORGANIC_CONTROLS M                 @                           V     '            #LMF W   #LMNF X   #LSF Y   #LSLF Z   #LSNF [                                           W            	                                           X           	                                           Y           	                                           Z           	                                           [           	                                           \        #ORGANIC_FRACTIONS V                 @                           ]     '            #CNR ^   #NCBM _   #NCHP `   #NCHS a                                           ^            	                                           _           	                                           `           	                                           a           	                                           b        #ORGANIC_RATIO ]                 @                           c     '0            #BMCTP d   #BMNTP e   #HSCTP f   #HSNTP g   #HPCTP h   #HPNTP i   #LMCTP j   #LMNTP k   #LSCTP l   #LSLCTP m   #LSLNCTP n   #LSNTP o                                           d            	                                           e           	                                           f           	                                           g           	                                           h           	                                           i           	                                           j           	                                           k           	                                           l         	   	                                           m     $   
   	                                           n     (      	                                           o     ,      	                                           p     0   #ORGANIC_TRANSFORMATIONS c                 @                           q     '      %      #CFMETS1 r   #CFSTRS1 s   #CFSTRS2 t   #EFMETS1 u   #EFSTRS1 v   #EFSTRS2 w   #IMMMETS1 x   #IMMSTRS1 y   #IMMSTRS2 z   #MNRMETS1 {   #MNRSTRS1 |   #MNRSTRS2 }   #CO2FMET ~   #CO2FSTR    #CFS1S2    #CFS1S3    #CFS2S1    #CFS2S3    #CFS3S1    #EFS1S2    #EFS1S3    #EFS2S1    #EFS2S3    #EFS3S1    #IMMS1S2    #IMMS1S3    #IMMS2S1    #IMMS2S3    #IMMS3S1    #MNRS1S2    #MNRS1S3    #MNRS2S1    #MNRS2S3    #MNRS3S1    #CO2FS1    #CO2FS2    #CO2FS3                                            r            	                                           s           	                                           t           	                                           u           	                                           v           	                                           w           	                                           x           	                                           y           	                                           z         	   	                                           {     $   
   	                                           |     (      	                                           }     ,      	                                           ~     0      	                                                4      	                                                8      	                                                <      	                                                @      	                                                D      	                                                H      	                                                L      	                                                P      	                                                T      	                                                X      	                                                \      	                                                `      	                                                d      	                                                h      	                                                l      	                                                p      	                                                t      	                                                x      	                                                |       	                                                   !   	                                                   "   	                                                   #   	                                                   $   	                                                   %   	                                                   #ORGANIC_FLUX q                 @                                '4            #META_MICR    #STR_MICR    #STR_HS    #CO2_META    #CO2_STR    #MICR_HS    #MICR_HP    #HS_MICR     #HS_HP ¡   #HP_MICR ¢   #CO2_MICR £   #CO2_HS ¤   #CO2_HP ¥                                                       	                                                      	                                                      	                                                      	                                                      	                                                      	                                                      	                                                       	                                           ¡         	   	                                           ¢     $   
   	                                           £     (      	                                           ¤     ,      	                                           ¥     0      	                                           ¦     4   #CARBON_SOIL_TRANSFORMATIONS                                           §        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                           ¨        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                           ©        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                           ª        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                           «        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                           ¬        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                           ­        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                           ®        4            &                       #CARBON_SOIL_TRANSFORMATIONS                                            ¯     4   #CARBON_SOIL_TRANSFORMATIONS                                            °     4   #CARBON_SOIL_TRANSFORMATIONS                                            ±     4   #CARBON_SOIL_TRANSFORMATIONS                                            ²     4   #CARBON_SOIL_TRANSFORMATIONS                  @                          ³     '<            #SED_C ´   #SURQ_C µ   #SURQ_DOC ¶   #SURQ_DIC ·   #LATQ_C ¸   #LATQ_DOC ¹   #LATQ_DIC º   #PERC_C »   #PERC_DOC ¼   #PERC_DIC ½   #RES_DECAY_C ¾   #MAN_APP_C ¿   #MAN_GRAZ_C À   #RSP_C Á   #EMIT_C Â                                          ´           	                             	                         0.                                       µ          	                             	                         0.                                       ¶          	                             	                         0.                                       ·          	                             	                         0.                                       ¸          	                             	                         0.                                       ¹          	                             	                         0.                                       º          	                             	                         0.                                       »          	                             	                         0.                                       ¼         	  	                             	                         0.                                       ½     $   
  	                             	                         0.                                       ¾     (     	                             	                         0.                                       ¿     ,     	                             	                         0.                                       À     0     	                             	                         0.                                       Á     4     	                             	                         0.                                       Â     8     	                             	                         0.                                        Ã     <   #CARBON_SOIL_GAIN_LOSSES ³                                          Ä        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                          Å        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                          Æ        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                          Ç        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                          È        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                          É        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                          Ê        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                          Ë        <            &                       #CARBON_SOIL_GAIN_LOSSES ³                                           Ì     <   #CARBON_SOIL_GAIN_LOSSES ³                                           Í     <   #CARBON_SOIL_GAIN_LOSSES ³                                           Î     <   #CARBON_SOIL_GAIN_LOSSES ³                                           Ï     <   #CARBON_SOIL_GAIN_LOSSES ³                 @                          Ð     '            #PLANT_C Ñ   #RES_DECAY_C Ò   #HARV_STOV_C Ó   #EMIT_C Ô                                          Ñ           	                             	                         0.                                       Ò          	                             	                         0.                                       Ó          	                             	                         0.                                       Ô          	                             	                         0.                                        Õ        #CARBON_RESIDUE_GAIN_LOSSES Ð                                          Ö                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                          ×                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                          Ø                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                          Ù                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                          Ú                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                          Û                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                          Ü                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                          Ý                    &                       #CARBON_RESIDUE_GAIN_LOSSES Ð                                           Þ        #CARBON_RESIDUE_GAIN_LOSSES Ð                                           ß        #CARBON_RESIDUE_GAIN_LOSSES Ð                                           à        #CARBON_RESIDUE_GAIN_LOSSES Ð                                           á        #CARBON_RESIDUE_GAIN_LOSSES Ð                 @                          â     '            #NPP_C ã   #HARV_C ä   #DROP_C å   #GRAZEAT_C æ   #EMIT_C ç                                          ã           	                             	                         0.                                       ä          	                             	                         0.                                       å          	                             	                         0.                                       æ          	                             	                         0.                                       ç          	                             	                         0.                                        è        #CARBON_PLANT_GAIN_LOSSES â                                          é                    &                       #CARBON_PLANT_GAIN_LOSSES â                                          ê                    &                       #CARBON_PLANT_GAIN_LOSSES â                                          ë                    &                       #CARBON_PLANT_GAIN_LOSSES â                                          ì                    &                       #CARBON_PLANT_GAIN_LOSSES â                                          í                    &                       #CARBON_PLANT_GAIN_LOSSES â                                          î                    &                       #CARBON_PLANT_GAIN_LOSSES â                                          ï                    &                       #CARBON_PLANT_GAIN_LOSSES â                                          ð                    &                       #CARBON_PLANT_GAIN_LOSSES â                                           ñ        #CARBON_PLANT_GAIN_LOSSES â                                           ò        #CARBON_PLANT_GAIN_LOSSES â                                           ó        #CARBON_PLANT_GAIN_LOSSES â                                           ô        #CARBON_PLANT_GAIN_LOSSES â   &     @     X                                 4                  #HRU1 õ   #HRU2 ö   #CARBON_SOIL_TRANSFORMATIONS          
                                  õ     4      #CARBON_SOIL_TRANSFORMATIONS          
                                  ö     4      #CARBON_SOIL_TRANSFORMATIONS    &     @     X                                 4                  #HRU1 ÷   #CONST ø   #CARBON_SOIL_TRANSFORMATIONS          
                                  ÷     4      #CARBON_SOIL_TRANSFORMATIONS          
                                  ø     	  &     @     X                            	     4                  #HRU1 ù   #CONST ú   #CARBON_SOIL_TRANSFORMATIONS          
                                  ù     4      #CARBON_SOIL_TRANSFORMATIONS          
                                  ú     	  &     @     X                                 <                  #HRU1 û   #HRU2 ü   #CARBON_SOIL_GAIN_LOSSES ³         
                                  û     <      #CARBON_SOIL_GAIN_LOSSES ³         
                                  ü     <      #CARBON_SOIL_GAIN_LOSSES ³   &     @     X                                 <                  #HRU1 ý   #CONST þ   #CARBON_SOIL_GAIN_LOSSES ³         
                                  ý     <      #CARBON_SOIL_GAIN_LOSSES ³         
                                  þ     	  &     @     X                            
     <                  #HRU1 ÿ   #CONST    #CARBON_SOIL_GAIN_LOSSES ³         
                                  ÿ     <      #CARBON_SOIL_GAIN_LOSSES ³         
                                       	  &     @     X                                                   #HRU1   #HRU2   #CARBON_RESIDUE_GAIN_LOSSES Ð         
                                            #CARBON_RESIDUE_GAIN_LOSSES Ð         
                                            #CARBON_RESIDUE_GAIN_LOSSES Ð   &     @     X                                                   #HRU1   #CONST   #CARBON_RESIDUE_GAIN_LOSSES Ð         
                                            #CARBON_RESIDUE_GAIN_LOSSES Ð         
                                      	  &     @     X                                                   #HRU1   #CONST   #CARBON_RESIDUE_GAIN_LOSSES Ð         
                                            #CARBON_RESIDUE_GAIN_LOSSES Ð         
                                      	  &     @     X                                                   #HRU1   #HRU2   #CARBON_PLANT_GAIN_LOSSES â         
                                            #CARBON_PLANT_GAIN_LOSSES â         
                                            #CARBON_PLANT_GAIN_LOSSES â   &     @     X                                                   #HRU1 	  #CONST 
  #CARBON_PLANT_GAIN_LOSSES â         
                                  	          #CARBON_PLANT_GAIN_LOSSES â         
                                  
    	  &     @     X                                                   #HRU1   #CONST   #CARBON_PLANT_GAIN_LOSSES â         
                                            #CARBON_PLANT_GAIN_LOSSES â         
                                      	         K      fn#fn    ë   ¦      i@      ¦      i@    7  ¢      i@ *   Ù  Ü      CARBON_TERRESTRIAL_INPUTS 6   µ     a   CARBON_TERRESTRIAL_INPUTS%ER_POC_PARA 3   4     a   CARBON_TERRESTRIAL_INPUTS%CFB_PARA 6   ´     a   CARBON_TERRESTRIAL_INPUTS%SF_PARA_SUR 6   4     a   CARBON_TERRESTRIAL_INPUTS%SF_PARA_SUB 3   ´     a   CARBON_TERRESTRIAL_INPUTS%ABL_PARA 9   3     a   CARBON_TERRESTRIAL_INPUTS%PEROC_DIC_PARA 9   ³     a   CARBON_TERRESTRIAL_INPUTS%PEROC_DOC_PARA 8   3	     a   CARBON_TERRESTRIAL_INPUTS%PART_DOC_PARA 9   ´	     a   CARBON_TERRESTRIAL_INPUTS%HLIFE_DOC_PARA 9   3
     a   CARBON_TERRESTRIAL_INPUTS%ABCO2_PARA_SUR 9   ²
  ~   a   CARBON_TERRESTRIAL_INPUTS%ABCO2_PARA_SUB 7   0     a   CARBON_TERRESTRIAL_INPUTS%ABP_PARA_SUR 7   ¯     a   CARBON_TERRESTRIAL_INPUTS%ABP_PARA_SUB :   .     a   CARBON_TERRESTRIAL_INPUTS%ALMCO2_PARA_SUR :   ­     a   CARBON_TERRESTRIAL_INPUTS%ALMCO2_PARA_SUB <   -     a   CARBON_TERRESTRIAL_INPUTS%ALSLNCO2_PARA_SUR <   ¬     a   CARBON_TERRESTRIAL_INPUTS%ALSLNCO2_PARA_SUB 7   ,     a   CARBON_TERRESTRIAL_INPUTS%ASP_PARA_SUR 7   «     a   CARBON_TERRESTRIAL_INPUTS%ASP_PARA_SUB 7   *     a   CARBON_TERRESTRIAL_INPUTS%ALSLCO2_PARA 5   ©     a   CARBON_TERRESTRIAL_INPUTS%APCO2_PARA 5   )     a   CARBON_TERRESTRIAL_INPUTS%ASCO2_PARA 7   ©     a   CARBON_TERRESTRIAL_INPUTS%PRMT_51_PARA 7   (     a   CARBON_TERRESTRIAL_INPUTS%PRMT_45_PARA 7   ©     a   CARBON_TERRESTRIAL_INPUTS%BMR_PARA_SUR 7   +     a   CARBON_TERRESTRIAL_INPUTS%BMR_PARA_SUB 3   «     a   CARBON_TERRESTRIAL_INPUTS%HPR_PARA 3   /     a   CARBON_TERRESTRIAL_INPUTS%HSR_PARA 7   ³     a   CARBON_TERRESTRIAL_INPUTS%LMR_PARA_SUR 7   5     a   CARBON_TERRESTRIAL_INPUTS%LMR_PARA_SUB 7   ·     a   CARBON_TERRESTRIAL_INPUTS%LSR_PARA_SUR 7   9     a   CARBON_TERRESTRIAL_INPUTS%LSR_PARA_SUB 7   »     a   CARBON_TERRESTRIAL_INPUTS%XBM_PARA_SUR 7   :     a   CARBON_TERRESTRIAL_INPUTS%XBM_PARA_SUB 5   ¹     a   CARBON_TERRESTRIAL_INPUTS%XLSLF_PARA 5   8     a   CARBON_TERRESTRIAL_INPUTS%OX_AA_PARA 5   ¸     a   CARBON_TERRESTRIAL_INPUTS%OX_BB_PARA    9  W       CBN_TES      Ü       CARBON_INPUTS &   l  ~   a   CARBON_INPUTS%HP_RATE &   ê  ~   a   CARBON_INPUTS%HS_RATE *   h  ~   a   CARBON_INPUTS%MICROB_RATE (   æ  ~   a   CARBON_INPUTS%META_RATE '   d  ~   a   CARBON_INPUTS%STR_RATE .   â  ~   a   CARBON_INPUTS%MICROB_TOP_RATE $   `  ~   a   CARBON_INPUTS%HS_HP )   Þ  ~   a   CARBON_INPUTS%MICROB_KOC )   \  ~   a   CARBON_INPUTS%MIN_N_FRAC )   Ú  ~   a   CARBON_INPUTS%C_ORG_FRAC    X  K       CARBDB    £  K       CARBZ $   î  §       ORGANIC_ALLOCATIONS *     @   a   ORGANIC_ALLOCATIONS%ABCO2 (   Õ  @   a   ORGANIC_ALLOCATIONS%ABL (      @   a   ORGANIC_ALLOCATIONS%ABP +   U   @   a   ORGANIC_ALLOCATIONS%ALMCO2 ,      @   a   ORGANIC_ALLOCATIONS%ALSLCO2 -   Õ   @   a   ORGANIC_ALLOCATIONS%ALSLNCO2 *   !  @   a   ORGANIC_ALLOCATIONS%APCO2 *   U!  @   a   ORGANIC_ALLOCATIONS%ASCO2 (   !  @   a   ORGANIC_ALLOCATIONS%ASP    Õ!  Q       ORG_ALLO    &"  Q       ORG_ALLOZ !   w"         ORGANIC_CONTROLS %   ú"  @   a   ORGANIC_CONTROLS%CDG $   :#  @   a   ORGANIC_CONTROLS%CS $   z#  @   a   ORGANIC_CONTROLS%OX %   º#  @   a   ORGANIC_CONTROLS%SUT $   ú#  @   a   ORGANIC_CONTROLS%X1 &   :$  @   a   ORGANIC_CONTROLS%XBMT '   z$  @   a   ORGANIC_CONTROLS%XLSLF    º$  N       ORG_CON "   %  t       ORGANIC_FRACTIONS &   |%  @   a   ORGANIC_FRACTIONS%LMF '   ¼%  @   a   ORGANIC_FRACTIONS%LMNF &   ü%  @   a   ORGANIC_FRACTIONS%LSF '   <&  @   a   ORGANIC_FRACTIONS%LSLF '   |&  @   a   ORGANIC_FRACTIONS%LSNF    ¼&  O       ORG_FRAC    '  k       ORGANIC_RATIO "   v'  @   a   ORGANIC_RATIO%CNR #   ¶'  @   a   ORGANIC_RATIO%NCBM #   ö'  @   a   ORGANIC_RATIO%NCHP #   6(  @   a   ORGANIC_RATIO%NCHS    v(  K       ORG_RATIO (   Á(  Ë       ORGANIC_TRANSFORMATIONS .   )  @   a   ORGANIC_TRANSFORMATIONS%BMCTP .   Ì)  @   a   ORGANIC_TRANSFORMATIONS%BMNTP .   *  @   a   ORGANIC_TRANSFORMATIONS%HSCTP .   L*  @   a   ORGANIC_TRANSFORMATIONS%HSNTP .   *  @   a   ORGANIC_TRANSFORMATIONS%HPCTP .   Ì*  @   a   ORGANIC_TRANSFORMATIONS%HPNTP .   +  @   a   ORGANIC_TRANSFORMATIONS%LMCTP .   L+  @   a   ORGANIC_TRANSFORMATIONS%LMNTP .   +  @   a   ORGANIC_TRANSFORMATIONS%LSCTP /   Ì+  @   a   ORGANIC_TRANSFORMATIONS%LSLCTP 0   ,  @   a   ORGANIC_TRANSFORMATIONS%LSLNCTP .   L,  @   a   ORGANIC_TRANSFORMATIONS%LSNTP    ,  U       ORG_TRAN    á,        ORGANIC_FLUX %   ÿ.  @   a   ORGANIC_FLUX%CFMETS1 %   ?/  @   a   ORGANIC_FLUX%CFSTRS1 %   /  @   a   ORGANIC_FLUX%CFSTRS2 %   ¿/  @   a   ORGANIC_FLUX%EFMETS1 %   ÿ/  @   a   ORGANIC_FLUX%EFSTRS1 %   ?0  @   a   ORGANIC_FLUX%EFSTRS2 &   0  @   a   ORGANIC_FLUX%IMMMETS1 &   ¿0  @   a   ORGANIC_FLUX%IMMSTRS1 &   ÿ0  @   a   ORGANIC_FLUX%IMMSTRS2 &   ?1  @   a   ORGANIC_FLUX%MNRMETS1 &   1  @   a   ORGANIC_FLUX%MNRSTRS1 &   ¿1  @   a   ORGANIC_FLUX%MNRSTRS2 %   ÿ1  @   a   ORGANIC_FLUX%CO2FMET %   ?2  @   a   ORGANIC_FLUX%CO2FSTR $   2  @   a   ORGANIC_FLUX%CFS1S2 $   ¿2  @   a   ORGANIC_FLUX%CFS1S3 $   ÿ2  @   a   ORGANIC_FLUX%CFS2S1 $   ?3  @   a   ORGANIC_FLUX%CFS2S3 $   3  @   a   ORGANIC_FLUX%CFS3S1 $   ¿3  @   a   ORGANIC_FLUX%EFS1S2 $   ÿ3  @   a   ORGANIC_FLUX%EFS1S3 $   ?4  @   a   ORGANIC_FLUX%EFS2S1 $   4  @   a   ORGANIC_FLUX%EFS2S3 $   ¿4  @   a   ORGANIC_FLUX%EFS3S1 %   ÿ4  @   a   ORGANIC_FLUX%IMMS1S2 %   ?5  @   a   ORGANIC_FLUX%IMMS1S3 %   5  @   a   ORGANIC_FLUX%IMMS2S1 %   ¿5  @   a   ORGANIC_FLUX%IMMS2S3 %   ÿ5  @   a   ORGANIC_FLUX%IMMS3S1 %   ?6  @   a   ORGANIC_FLUX%MNRS1S2 %   6  @   a   ORGANIC_FLUX%MNRS1S3 %   ¿6  @   a   ORGANIC_FLUX%MNRS2S1 %   ÿ6  @   a   ORGANIC_FLUX%MNRS2S3 %   ?7  @   a   ORGANIC_FLUX%MNRS3S1 $   7  @   a   ORGANIC_FLUX%CO2FS1 $   ¿7  @   a   ORGANIC_FLUX%CO2FS2 $   ÿ7  @   a   ORGANIC_FLUX%CO2FS3    ?8  J       ORG_FLUX ,   8  í       CARBON_SOIL_TRANSFORMATIONS 6   v9  @   a   CARBON_SOIL_TRANSFORMATIONS%META_MICR 5   ¶9  @   a   CARBON_SOIL_TRANSFORMATIONS%STR_MICR 3   ö9  @   a   CARBON_SOIL_TRANSFORMATIONS%STR_HS 5   6:  @   a   CARBON_SOIL_TRANSFORMATIONS%CO2_META 4   v:  @   a   CARBON_SOIL_TRANSFORMATIONS%CO2_STR 4   ¶:  @   a   CARBON_SOIL_TRANSFORMATIONS%MICR_HS 4   ö:  @   a   CARBON_SOIL_TRANSFORMATIONS%MICR_HP 4   6;  @   a   CARBON_SOIL_TRANSFORMATIONS%HS_MICR 2   v;  @   a   CARBON_SOIL_TRANSFORMATIONS%HS_HP 4   ¶;  @   a   CARBON_SOIL_TRANSFORMATIONS%HP_MICR 5   ö;  @   a   CARBON_SOIL_TRANSFORMATIONS%CO2_MICR 3   6<  @   a   CARBON_SOIL_TRANSFORMATIONS%CO2_HS 3   v<  @   a   CARBON_SOIL_TRANSFORMATIONS%CO2_HP    ¶<  Y       HSCFZ    =         HSCF_D    =         HSCF_M    >         HSCF_Y    >         HSCF_A    #?         LSCF_D    ¨?         LSCF_M    -@         LSCF_Y    ²@         LCSF_A    7A  Y       BSCF_D    A  Y       BSCF_M    éA  Y       BSCF_Y    BB  Y       BSCF_A (   B        CARBON_SOIL_GAIN_LOSSES .   ©C  ~   a   CARBON_SOIL_GAIN_LOSSES%SED_C /   'D  ~   a   CARBON_SOIL_GAIN_LOSSES%SURQ_C 1   ¥D  ~   a   CARBON_SOIL_GAIN_LOSSES%SURQ_DOC 1   #E  ~   a   CARBON_SOIL_GAIN_LOSSES%SURQ_DIC /   ¡E  ~   a   CARBON_SOIL_GAIN_LOSSES%LATQ_C 1   F  ~   a   CARBON_SOIL_GAIN_LOSSES%LATQ_DOC 1   F  ~   a   CARBON_SOIL_GAIN_LOSSES%LATQ_DIC /   G  ~   a   CARBON_SOIL_GAIN_LOSSES%PERC_C 1   G  ~   a   CARBON_SOIL_GAIN_LOSSES%PERC_DOC 1   H  ~   a   CARBON_SOIL_GAIN_LOSSES%PERC_DIC 4   H  ~   a   CARBON_SOIL_GAIN_LOSSES%RES_DECAY_C 2   I  ~   a   CARBON_SOIL_GAIN_LOSSES%MAN_APP_C 3   I  ~   a   CARBON_SOIL_GAIN_LOSSES%MAN_GRAZ_C .   J  ~   a   CARBON_SOIL_GAIN_LOSSES%RSP_C /   J  ~   a   CARBON_SOIL_GAIN_LOSSES%EMIT_C    K  U       HSCZ    `K         HSC_D    áK         HSC_M    bL         HSC_Y    ãL         HSC_A    dM         LSC_D    åM         LSC_M    fN         LSC_Y    çN         LCS_A    hO  U       BSC_D    ½O  U       BSC_M    P  U       BSC_Y    gP  U       BSC_A +   ¼P         CARBON_RESIDUE_GAIN_LOSSES 3   ;Q  ~   a   CARBON_RESIDUE_GAIN_LOSSES%PLANT_C 7   ¹Q  ~   a   CARBON_RESIDUE_GAIN_LOSSES%RES_DECAY_C 7   7R  ~   a   CARBON_RESIDUE_GAIN_LOSSES%HARV_STOV_C 2   µR  ~   a   CARBON_RESIDUE_GAIN_LOSSES%EMIT_C    3S  X       HRCZ    S         HRC_D    T         HRC_M    T         HRC_Y    U         HRC_A    U         LRC_D    V         LRC_M    £V         LRC_Y    'W         LRS_A    «W  X       BRC_D    X  X       BRC_M    [X  X       BRC_Y    ³X  X       BRC_A )   Y         CARBON_PLANT_GAIN_LOSSES /   Y  ~   a   CARBON_PLANT_GAIN_LOSSES%NPP_C 0   Z  ~   a   CARBON_PLANT_GAIN_LOSSES%HARV_C 0   Z  ~   a   CARBON_PLANT_GAIN_LOSSES%DROP_C 3   [  ~   a   CARBON_PLANT_GAIN_LOSSES%GRAZEAT_C 0   [  ~   a   CARBON_PLANT_GAIN_LOSSES%EMIT_C    \  V       HPCZ    Y\         HPC_D    Û\         HPC_M    ]]         HPC_Y    ß]         HPC_A    a^         LPC_D    ã^         LPC_M    e_         LPC_Y    ç_         LPS_A    i`  V       BPC_D    ¿`  V       BPC_M    a  V       BPC_Y    ka  V       BPC_A &   Áa  }       CARBON_SOIL_FLUX__ADD +   >b  ]   a   CARBON_SOIL_FLUX__ADD%HRU1 +   b  ]   a   CARBON_SOIL_FLUX__ADD%HRU2 &   øb  ~       CARBON_SOIL_FLUX_MULT +   vc  ]   a   CARBON_SOIL_FLUX_MULT%HRU1 ,   Óc  8   a   CARBON_SOIL_FLUX_MULT%CONST %   d  ~       CARBON_SOIL_FLUX_DIV *   d  ]   a   CARBON_SOIL_FLUX_DIV%HRU1 +   æd  8   a   CARBON_SOIL_FLUX_DIV%CONST $   e  y       CARBON_SOIL_GL__ADD )   e  Y   a   CARBON_SOIL_GL__ADD%HRU1 )   ðe  Y   a   CARBON_SOIL_GL__ADD%HRU2 $   If  z       CARBON_SOIL_GL_MULT )   Ãf  Y   a   CARBON_SOIL_GL_MULT%HRU1 *   g  8   a   CARBON_SOIL_GL_MULT%CONST #   Tg  z       CARBON_SOIL_GL_DIV (   Îg  Y   a   CARBON_SOIL_GL_DIV%HRU1 )   'h  8   a   CARBON_SOIL_GL_DIV%CONST '   _h  |       CARBON_RESIDUE_GL__ADD ,   Ûh  \   a   CARBON_RESIDUE_GL__ADD%HRU1 ,   7i  \   a   CARBON_RESIDUE_GL__ADD%HRU2 '   i  }       CARBON_RESIDUE_GL_MULT ,   j  \   a   CARBON_RESIDUE_GL_MULT%HRU1 -   lj  8   a   CARBON_RESIDUE_GL_MULT%CONST &   ¤j  }       CARBON_RESIDUE_GL_DIV +   !k  \   a   CARBON_RESIDUE_GL_DIV%HRU1 ,   }k  8   a   CARBON_RESIDUE_GL_DIV%CONST %   µk  z       CARBON_PLANT_GL__ADD *   /l  Z   a   CARBON_PLANT_GL__ADD%HRU1 *   l  Z   a   CARBON_PLANT_GL__ADD%HRU2 %   ãl  {       CARBON_PLANT_GL_MULT *   ^m  Z   a   CARBON_PLANT_GL_MULT%HRU1 +   ¸m  8   a   CARBON_PLANT_GL_MULT%CONST $   ðm  {       CARBON_PLANT_GL_DIV )   kn  Z   a   CARBON_PLANT_GL_DIV%HRU1 *   Ån  8   a   CARBON_PLANT_GL_DIV%CONST 