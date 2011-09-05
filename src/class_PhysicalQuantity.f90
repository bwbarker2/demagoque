module class_PhysicalQuantity
 use prec_def
 implicit none

 private :: MAX_SYMBOL_LENGTH, NUM_BASE_UNITS

 integer, parameter :: MAX_SYMBOL_LENGTH = 5
 integer, parameter :: NUM_BASE_UNITS = 7

 !> \brief Stores a quantity with a value, uncertainty, and units
 !!
 type PhysicalQuantity
  character(MAX_SYMBOL_LENGTH) :: symbol
  !> Includes units of quantity in terms of base SI units in following order:
  !!
  !! -# length [meter]
  !! -# mass [kilogram]
  !! -# time [second]
  !! -# electric current [ampere]
  !! -# thermodynamic temperature [kelvin]
  !! -# amount of substance [mole]
  !! -# luminous intensity [candela]
  real (Long), dimension(NUM_BASE_UNITS) :: units
  real (Long) :: val
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

 function new_PhysicalQuantity(newSymbol, newUnits, newVal, newUnc) result (this)
  use prec_def
  implicit none

  character(*), intent(in) :: newSymbol !< symbol of quantity
  real (Long), dimension(:), intent(in) :: newUnits  !< SI units of quantity
  real (Long), optional,     intent(in) :: newVal    !< value of quantity
  real (Long), optional,     intent(in) :: newUnc    !< uncertainty of quantity
  type (physicalQuantity) :: this

  if(size(newUnits) /= NUM_BASE_UNITS) then
   write(stderr,*)'class_PhysicalQuantity.f90:make_physicalQuantity: newUnits array not correct size'
   return
  endif

  if(present(newVal)) then
   if(present(newUnc)) then
    this = PhysicalQuantity(newSymbol,newUnits,newVal,newUnc)
   else
    this = PhysicalQuantity(newSymbol,newUnits,newVal,0_Long)
   endif
  else
   this = PhysicalQuantity(newSymbol,newUnits,1_Long,0_Long)
  endif

 end function new_physicalQuantity

end module class_PhysicalQuantity

