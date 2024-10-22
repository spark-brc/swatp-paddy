      subroutine res_weir_release (jres, id, ihyd, pvol_m3, evol_m3, dep, weir_hgt)
      !! this subroutine calculates weir discharge from a wetland/rice paddy HRU called from manual operation schedule.  

      use reservoir_data_module
      use reservoir_module
      use conditional_module
      use climate_module
      use time_module
      use hydrograph_module
      use water_body_module
      use soil_module
      use hru_module
      use water_allocation_module
      use basin_module
      
      implicit none
      
      real,  intent (in) :: pvol_m3
      real,  intent (in) :: evol_m3
      real,  intent (in) :: dep       !m 
      real,  intent (in) :: weir_hgt  !m         |height of weir overflow crest from reservoir bottom
      integer,  intent (in) :: jres             !none      |hru number
      integer :: iweir             !none      |weir ID 
      integer :: nstep            !none      |counter
      integer :: tstep            !none      |hru number
      integer :: iac              !none      |counter 
      integer :: ic              !none      |counter
      integer,  intent (in) :: id               !none      |hru number
      integer :: ial              !none      |counter
      integer :: irel             !          |
      integer :: iob              !none      |hru or wro number
      integer,  intent (in) :: ihyd             !          |
      real :: vol,vol_above                 !          |
      real :: b_lo                !          |
      character(len=1) :: action  !          |
      real :: res_h               !m         |water depth
      real :: demand              !m3        |irrigation demand by hru or wro
      real :: wsa1                !m2        |water surface area 
      real :: qout                !m3        |weir discharge during short time step
      real :: hgt                 !m         |height of bottom of weir above bottom of impoundment
      real :: hgt_above           !m         |height of water above the above bottom of weir
      real :: sto_max             !m3        |maximum storage volume at the bank top
      
      !! store initial values
      vol = wbody%flo
      nstep = 1
      iweir = wet_ob(jres)%iweir
      vol_above = 0 !water storage above weir height
      
      wsa1 = wbody_wb%area_ha * 10000. !m2      
      hgt_above = max(0., dep - weir_hgt)  !m ponding depth above weir crest  

      do tstep = 1, nstep
          
        !! calculate weir discharge from scheduled management
        if (hgt_above > 0 .and. iweir > 0) then

          !emergency spillway discharge Jaehak 2023
          if (vol>evol_m3) then
            ht2%flo = ht2%flo + (wbody%flo - evol_m3) 
            ht2%flo = max(0.,ht2%flo)
            vol = evol_m3
            res_h = vol / wsa1 !m
            hgt_above = max(0.,res_h - weir_hgt)  !m 
          endif
          vol_above = hgt_above * wsa1 !m3

          if (nstep>1) then !revised by Jaehak 2023
            qout = res_weir(iweir)%c * res_weir(iweir)%w * hgt_above ** res_weir(iweir)%k !m3/s
            qout = max(0.,86400. / nstep * qout) !m3
            if (qout > vol_above) then
              ht2%flo = ht2%flo + vol_above !weir discharge volume for the day, m3
              vol = vol - vol_above
              vol_above = 0.
            else
              ht2%flo = ht2%flo + qout 
              vol = vol - qout
              vol_above = vol_above - qout
            end if
            res_h = vol / wsa1 !m
            hgt_above = max(0.,res_h - weir_hgt)  !m Jaehak 2022
            if (vol_above<=0.001.or.hgt_above<=0.0001) exit

          else
            do ic = 1, 24
              qout = res_weir(iweir)%c * res_weir(iweir)%w * hgt_above ** res_weir(iweir)%k !m3/s
              qout = 3600. * qout !m3
              if (qout > vol_above) then
                ht2%flo = ht2%flo + vol_above !weir discharge volume for the day, m3
                vol = vol - vol_above
                vol_above = 0.
              else
                ht2%flo = ht2%flo + qout 
                vol = vol - vol_above
                vol_above = vol_above - qout
              end if
                          
              if (wsa1 > 1.e-6) then
                res_h = vol / wsa1 !m
              else
                res_h = 0.
              end if
              hgt_above = max(0.,res_h - weir_hgt)  !m Jaehak 2022
              if (vol_above<=0.001.or.hgt_above<=0.0001) exit
            end do
          endif
        endif
        wbody%flo = vol !m3

      end do  

      return
      end subroutine res_weir_release