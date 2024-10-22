      subroutine wet_initial (iihru)
      
      use reservoir_module
      use reservoir_data_module
      use reservoir_data_module
      use hydrograph_module
      use hru_module, only : hru, ihru
      use maximum_data_module
      use water_body_module
      use soil_module
      
      implicit none
      
      integer, intent (in) :: iihru     !none       |
      integer :: init_om
      integer :: iprop          !              | 
      integer :: ihyd           !none          |counter 
      integer :: init           !              |
      integer :: iweir          !              |
      real :: cnv               !none          |conversion factor (mm => m^3)
      real :: x1                !              |
      real :: wet_h             !              |
      real :: wet_h1            !              | 
      real :: wet_fr            !              | 
  
      iprop = hru(iihru)%dbs%surf_stor
      iweir = wet_ob(iihru)%iweir
      
      if (iprop > 0) then
        ihyd = wet_dat(iprop)%hyd
        if (wet_hyd(ihyd)%k > 0) then
          hru(iihru)%wet_hc = wet_hyd(ihyd)%k  !mm/hr
        else
          hru(iihru)%wet_hc = soil(iihru)%phys(1)%k 
        endif
        !! ha*mm*10. => m**3  - assume entire hru is wet and don't use fractional inputs (for simplicity)
        wet_ob(iihru)%evol = hru(iihru)%area_ha * wet_hyd(ihyd)%edep * 10.  ! * wet_hyd(ihyd)%esa
        wet_ob(iihru)%pvol = hru(iihru)%area_ha * wet_hyd(ihyd)%pdep * 10.  ! * wet_hyd(ihyd)%psa
        wet_ob(iihru)%psa = wet_hyd(ihyd)%psa * hru(iihru)%area_ha 
        wet_ob(iihru)%esa = wet_hyd(ihyd)%esa * hru(iihru)%area_ha
        !! set initial weir height to principal depth - m
        if (db_mx%res_weir > 0.and.wet_ob(iihru)%iweir>0) then !if available, read from weir.res Jaehak 2022
          wet_ob(iihru)%weir_hgt = res_weir(iweir)%h  !m weir height
          wet_ob(iihru)%weir_wid = res_weir(iweir)%w  !m, weir width
          
          !update pvol/evol according to weir height for paddy weir discharge. Jaehak 2023
          wet_ob(iihru)%pvol = hru(iihru)%area_ha * wet_ob(iihru)%weir_hgt * 10.**4  ! m3
          if (wet_ob(iihru)%evol < wet_ob(iihru)%pvol*1.2) then
            wet_ob(iihru)%evol = wet_ob(iihru)%pvol * 1.2   
          endif
          
        else
          wet_ob(iihru)%weir_hgt = wet_hyd(ihyd)%pdep / 1000.  !m 
          wet_ob(iihru)%weir_wid = 2.5  !m
        endif
         
        !Weir width reset to drain 90% in 3 days for 150mm weir. Special arrangement for RDA's Global Initiative
        wet_ob(iihru)%weir_wid =  hru(iihru)%area_ha * 10.
        
        
        !!set initial n and p concentrations --> (ppm) * (m^3) / 1000 = kg    !! ppm = t/m^3 * 10^6
        init = wet_dat(iprop)%init
        init_om = wet_init(init)%org_min
        cnv = om_init_water(init_om)%flo * wet_ob(iihru)%pvol / 1000.
        wet(iihru) = cnv * om_init_water(init_om)
        !wet(iihru)%flo = om_init_water(init_om)%flo * wet_ob(iihru)%pvol !Jaehak 2022
        wet_om_init(iihru) = wet(iihru)
  
        !! update surface area
        !! wetland on hru - solve quadratic to find new depth
        wet_wat_d(iihru)%area_ha = 0.
        if (wet(iihru)%flo > 0.) then
          x1 = wet_hyd(ihyd)%bcoef ** 2 + 4. * wet_hyd(ihyd)%ccoef * (1. - wet(iihru)%flo / wet_ob(iihru)%pvol)
          if (x1 < 1.e-6) then
            wet_h = 0.
          else
            wet_h1 = (-wet_hyd(ihyd)%bcoef - sqrt(x1)) / (2. * wet_hyd(ihyd)%ccoef)
            wet_h = wet_h1 + wet_hyd(ihyd)%bcoef
          end if
          wet_fr = (1. + wet_hyd(ihyd)%acoef * wet_h)
          wet_fr = min(wet_fr,1.)
          wet_wat_d(iihru)%area_ha = hru(iihru)%area_ha * wet_hyd(ihyd)%psa * wet_fr
        end if 
  
      end if
      

      return
      end subroutine wet_initial