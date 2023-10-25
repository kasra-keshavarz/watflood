    module area_debug

!   If debugOutput = false then in the release version the debug code marked with this will not be compiled
!   and this code will not slow down execution.
!   So set = true for debug compile,  false for release compile
    
        logical, parameter :: debug_output = .false.
!        logical, parameter :: debug_output = .true.
    
    end module area_debug
