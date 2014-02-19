MODULE ICE_ALBEDO
  ! Compute the albedo of the ice using different parameterization that rely only 
  ! on the surface temperature [K] and the ice thickness [m].
  implicit none
  real, parameter :: fresh_water_freezing_point=273.15
  integer, parameter:: RP = selected_real_kind(12 ) 
  real(RP), parameter :: albedoo = 0.06   ! Ocean albedo

CONTAINS
  

  elemental  function albedo_Hibler_1980(temperature) result (albedo)
    ! Return the ice albedo given the surface temperature [K].
    !
    ! When the temperature is below the freezing point, the ice albedo is set to .75
    ! otherwise, it is at .616.
    ! 
    ! :Input:
    ! temperature : real
    !   Surface temperature [K].
    !
    ! :Output:
    ! albedo : real 
    ! Surface albedo [0,1]. 
    !
    ! Reference
    !
    ! Hibler, W.D, Modeling a variable thickness sea ice cover, Mon. Weather Review, 108, 1943-1973, 1980. 
    
    real(RP), intent(in) :: temperature
    real(RP) :: albedo

    if (temperature < fresh_water_freezing_point) then
       albedo = 0.75
    else
       albedo = .616
    end if
  end function albedo_hibler_1980
    

  elemental function albedo_Ingram(temperature) result (albedo)
    ! Return the ice albedo given the surface temperature [K]. 
    !
    ! This parameterization of albedo attempts to capture the effect of 
    ! snow ageing. 
    ! 
    ! :Input:
    ! temperature : real
    !   Surface temperature [K].
    !
    ! :Output:
    ! albedo : real 
    ! Surface albedo [0,1]. 
    !
    ! Notes
    !
    ! The freezing point used in this parameterization is the ocean water
    ! freezing temperature. 
    !
    ! Reference
    !
    ! Ingram, W.J., C.A. Wilson and J.F.B. Mitchell, Modeling climate change: An assessment of sea ice
    ! and surface albedo feedbacks, J. Geophys. REs., 94, 8609-8622, 1989. 
    
    real(RP), intent(in) :: temperature
    real(RP) :: albedo

    if (temperature <= 261.2) then
       albedo = 0.7
    elseif (temperature < 271.2) then
       albedo = .7 - 0.03*(temperature - 261.2)
    else
       albedo = .4
    end if
  end function albedo_Ingram
    
    


  elemental function albedo_Manabe(temperature, thickness) result (albedo)
    ! Return the ice albedo given the surface temperature [K] and sea ice thickness. 
    !
    ! This function is a modified version of the one given in the reference (although 
    ! I couldn't find the equations in the paper.)
    !
    !
    ! :Input:
    ! temperature : real
    !   Surface temperature [K].
    ! thickness : real
    !   Ice thickness [m].
    !
    !
    ! :Output:
    ! albedo : real 
    ! Surface albedo [0,1]. 
    !
    ! :Reference:
    ! Manabe, S., M.J. Spellman and R.J. Stouffer, Transient responses of a coupled ocean-atmosphere 
    ! model to gradual changes of the atmospheric CO_2, II, Seasonal response, J. Clim., 5, 105-126, 1992.
    
    real(RP), intent(in) :: temperature, thickness
    real(RP) :: albedo
    real(RP), parameter :: melting_temperature=273.15
    real(RP), parameter :: ice_albedo = 0.5
    real(RP), parameter :: snow_albedo = 0.8

    if (temperature >= melting_temperature) then
       albedo = ice_albedo
    elseif (temperature > melting_temperature - 10) then
       albedo = ice_albedo + (snow_albedo-ice_albedo)*(melting_temperature - temperature)/10.
    else
       albedo = snow_albedo
    end if


    if (thickness < 1.) then
       albedo = sqrt(thickness) * (albedo - albedoo) + albedoo
    end if
  end function albedo_Manabe
    


end MODULE ICE_ALBEDO



