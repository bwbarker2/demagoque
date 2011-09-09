module class_PhysicalQuantity
 use class_UnitSystem
 use prec_def
 implicit none

 private :: MAX_SYMBOL_SIZE 

 integer, parameter :: MAX_SYMBOL_SIZE = 5

 !> This is the unit system that all quantities are defined in.
 !! By default, we use SI units.
 type(UnitSystem) :: theUnitSystem = UNIT_SYSTEM_SI

 !> \brief Stores a quantity with a value, uncertainty, and units
 !!
 type PhysicalQuantity

  character(MAX_SYMBOL_SIZE) :: symbol

  !> Includes units of quantity in terms of base SI units in order of
  !! UNIT_SYSTEM_DIMENSION_NAMES, currently given here:
  !!
  !! -# length
  !! -# mass
  !! -# time
  !! -# electric current
  !! -# thermodynamic temperature
  !! -# amount of substance
  !! -# luminous intensity
  !!
  !! For example, "centimeter" would have dimensions = (/1,0,0,0,0,0,0/)
  real (Long), dimension(UNIT_SYSTEM_NUM_BASE_UNITS) :: dimensions

  !> value of quantity
  real (Long) :: val

  !< uncertainty of quantity. 
  !! Assuming this is one standard deviation away from value in either
  !! direction
  real (Long) :: unc 
 end type physicalQuantity

 !base SI units

 type (PhysicalQuantity), parameter :: &
      METER    = PhysicalQuantity("m",  (/1,0,0,0,0,0,0/),1,0), &
      KILOGRAM = PhysicalQuantity("kg", (/0,1,0,0,0,0,0/),1,0), &
      SECOND   = PhysicalQuantity("s",  (/0,0,1,0,0,0,0/),1,0), &
      AMPERE   = PhysicalQuantity("A",  (/0,0,0,1,0,0,0/),1,0), &
      KELVIN   = PhysicalQuantity("K",  (/0,0,0,0,1,0,0/),1,0), &
      MOLE     = PhysicalQuantity("mol",(/0,0,0,0,0,1,0/),1,0), &
      CANDELA  = PhysicalQuantity("cd", (/0,0,0,0,0,0,1/),1,0)

contains

 !> Constructs new PhysicalQuantity.
 !! Assumes the value and uncertainty are given in the current
 !! unit system.
 function new_PhysicalQuantity(newSymbol, newDimensions, newVal, newUnc) result (this)
  use prec_def
  implicit none

  character(*), intent(in) :: newSymbol !< symbol of quantity
  real (Long), dimension(:), intent(in) :: newDimensions  !< new dimensions
  real (Long), optional,     intent(in) :: newVal    !< value of quantity
  real (Long), optional,     intent(in) :: newUnc    !< uncertainty of quantity
  type (physicalQuantity) :: this

  if(size(newDimensions) /= UNIT_SYSTEM_NUM_BASE_UNITS) then
   write(stderr,*)'class_PhysicalQuantity.f90:new_physicalQuantity: newDimensions array not correct size'
   return
  endif

  if(present(newVal)) then
   if(present(newUnc)) then
    this = PhysicalQuantity(newSymbol,newDimensions,newVal,newUnc)
   else
    this = PhysicalQuantity(newSymbol,newDimensions,newVal,0_Long)
   endif
  else
   this = PhysicalQuantity(newSymbol,newDimensions,1_Long,0_Long)
  endif

 end function new_physicalQuantity

end module class_PhysicalQuantity

