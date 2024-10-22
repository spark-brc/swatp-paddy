      subroutine res_hydro (jres, id, ihyd, pvol_m3, evol_m3, dep, weir_hgt)

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
      integer :: weir_flg=0        !none      |counter
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
      wsa1 = wbody_wb%area_ha * 10000. !m2
      
      do tstep = 1, nstep

        !calc release from decision table
        do iac = 1, d_tbl%acts
          action = "n"
          if (d_tbl%alts == 0) then
            action = "y"
          else
            do ial = 1, dtbl_res(id)%alts
              if (d_tbl%act_hit(ial) == "y" .and. d_tbl%act_outcomes(iac,ial) == "y") then
                action = "y"
                exit
              end if
            end do
          end if
          
          !condition is met - set the release rate
          if (action == "y") then
            select case (d_tbl%act(iac)%option)
            case ("rate")
              ht2%flo = d_tbl%act(iac)%const * 86400.
              
            case ("inflo_rate")
              ht2%flo = amax1 (ht1%flo, d_tbl%act(iac)%const * 86400.)
              
            case ("days")
              select case (dtbl_res(id)%act(iac)%file_pointer)
                case ("null")
                  b_lo = 0.
                case ("pvol")
                  b_lo = pvol_m3
                case ("evol")
                  b_lo = evol_m3
              end select
              ht2%flo = ht2%flo + (wbody%flo - b_lo) / d_tbl%act(iac)%const
              ht2%flo = max(0.,ht2%flo)
              wbody%flo = max(0.,wbody%flo - ht2%flo)
              vol = wbody%flo
            case ("dyrt")
              !for base volume for drawdown days, use condition associated with action
              select case (d_tbl%act(iac)%file_pointer)
                case ("con1")
                  ic = 1
                case ("con2")
                  ic = 2
                case ("con3")
                  ic = 3
                case ("con4")
                  ic = 4
              end select
              !perform operation on target variable to get target
              select case ((d_tbl%cond(ic)%lim_op))
              case ('=') 
                b_lo = pvol_m3 + (evol_m3 - pvol_m3) * d_tbl%cond(ic)%lim_const
              case ("*")
                b_lo = (evol_m3 - pvol_m3) * d_tbl%cond(ic)%lim_const
              case ("+")
                b_lo = (evol_m3 - pvol_m3) + d_tbl%cond(ic)%lim_const
              case ("-")
                b_lo = (evol_m3 - pvol_m3) - d_tbl%cond(ic)%lim_const
              case ("/")
                b_lo = (evol_m3 - pvol_m3) / d_tbl%cond(ic)%lim_const
              end select
              ht2%flo = (wbody%flo - b_lo) / d_tbl%act(iac)%const + d_tbl%act(iac)%const2 * pvol_m3 / 100.
              
            case ("irrig_dmd")
              iob = Int(d_tbl%act(iac)%const2)
              select case (d_tbl%act(iac)%file_pointer)
              case ("wro")    !demand from water resource object
                demand = wallo(iob)%tot%demand
              case ("hru")    !demand from single hru
                demand = irrig(iob)%demand
              end select
              !! const allows a fraction (usually > 1.0) of the demand (m3) released
              ht2%flo = demand * d_tbl%act(iac)%const
                 
            case ("weir")
              iweir = d_tbl%act_typ(iac)  !Jaehak 2023
              res_h = vol / wsa1 !m
              hgt_above = max(0.,res_h - weir_hgt)  !m Jaehak 2022

              if (nstep>1) then
                ht2%flo = res_weir(iweir)%c * res_weir(iweir)%w * hgt_above ** res_weir(iweir)%k !m3/s
                ht2%flo = max(0.,86400. / nstep * ht2%flo) !m3
                vol = vol - ht2%flo
              else
                do ic = 1, 24
                  vol_above = hgt_above * wsa1 !m3 water volume above weir height
                  qout = res_weir(iweir)%c * res_weir(iweir)%w * hgt_above ** res_weir(iweir)%k !m3/s
                  qout = 3600. * qout !m3
                  if (qout > vol_above) then
                    ht2%flo = ht2%flo + vol_above !weir discharge volume for the day, m3
                    vol = vol - vol_above
                  else
                    ht2%flo = ht2%flo + qout 
                    vol = vol - qout
                  end if
                              
                  res_h = vol / wsa1 !m
                  hgt_above = max(0.,res_h - weir_hgt)  !m Jaehak 2022
                  if (vol_above<=0.001.or.hgt_above<=0.0001) exit
                  
                end do

              endif
              res_h = vol / wsa1 !m
              wbody%flo = vol !m3

            case ("meas")
              irel = int(d_tbl%act_typ(iac))
              select case (recall(irel)%typ)
              case (1)    !daily
                ht2%flo = recall(irel)%hd(time%day,time%yrs)%flo
              case (2)    !monthly
                ht2%flo = recall(irel)%hd(time%mo,time%yrs)%flo
              case (3)    !annual
                ht2%flo = recall(irel)%hd(1,time%yrs)%flo
              end select
            end select
          end if    ! if action hit 
            
        end do      ! iac - actions loop
      end do    !tstep loop

      return
      end subroutine res_hydro