!> formatting settings for output
module formatting
 use iso_varying_string
 implicit none

 character(len=20), parameter :: fr5 = "(5E17.9)"

 !> output directory, of form 'rel/path/' or '/abs/path/'
 type(varying_string) :: fout_dir

 !> prefix for output files that depend on mode of evolution
 type(varying_string) :: fout_mode_prefix

 !> complete prefix for evolution output files
 type(varying_string) :: fout_ev_pre

contains

 subroutine setModePrefix(pre)
  implicit none

  character(len=*), intent(in) :: pre

  fout_mode_prefix=pre

  fout_ev_pre=fout_dir//fout_mode_prefix

 end subroutine setModePrefix

end module formatting


