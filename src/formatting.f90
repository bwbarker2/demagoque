!> formatting settings for output
module formatting
 use iso_varying_string
 implicit none

 character(len=20), parameter :: fr5 = "(5E17.9)"

 !> prefix for output files that depend on mode of evolution
 type(varying_string) :: fout_mode_prefix

end module formatting


