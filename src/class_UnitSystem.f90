module class_UnitSystem
 use prec_def
 implicit none

 private :: MAX_NAME_SIZE

 integer, parameter :: MAX_NAME_SIZE = 10

 integer, parameter :: UNIT_SYSTEM_NUM_BASE_UNITS = 7

 character(len=11), dimension(UNIT_SYSTEM_NUM_BASE_UNITS), parameter :: &
      UNIT_SYSTEM_DIMENSION_NAMES = (/"length     ", &
                                      "mass       ", &
                                      "time       ", &
                                      "current    ", &
                                      "temperature", &
                                      "amount     ", &
                                      "intensity  " /)
         

 type UnitSystem

  !> name of system of units (SI, Nuclear, BEC, etc)
  character(MAX_NAME_SIZE) :: name

  !> scale factor compared to SI unit (e.g. femtometer would be 1e-15_Long)
  real (Long), dimension(UNIT_SYSTEM_NUM_BASE_UNITS) :: scaleSI

 end type UnitSystem

 type(UnitSystem), parameter :: &
      UNIT_SYSTEM_SI = UnitSystem("SI",(/1_Long,1_Long,1_Long &
                                        ,1_Long,1_Long,1_Long,1_Long/))

end module class_UnitSystem
