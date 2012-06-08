!> Tests class_WfSquareWell for expected energy and norm.
!!
!! Expectations are from a sample run at initial coding, so if this test fails,
!! it is different from original output.
program class_WfSquareWell_run
 use class_WfSquareWell
 use phys_cons
 use prec_def
 implicit none

 type(WfSquareWell) :: testState

! integer :: ix
! real(Long) :: xlo=-25
! real(Long) :: xhi=25
! integer :: nx=100
! real(Long) :: xx
 real(Long) :: eps=1e-5_Long !tolerance, in fractional error from expected
 real(Long) :: energyExp=7.85219860582229331e-036_Long !expected energy
 real(Long) :: normExp=95.301443332921053 !expected norm

 call phys_cons_initializeBEC

 testState=make_WfSquareWell(MASS_REL_RUBIDIUM_87,0,1e10_Long*hbar**2/(2e0_Long*MASS_REL_RUBIDIUM_87*meter**2),1e-4_Long*meter)

! write(*,*)testState

 ! test for correct energy. Divide by Joule to get dimensionless value.
 if(abs((testState%energy/joule-energyExp)/energyExp) > eps) then
  call throwException('class_WfSquareWell_run: Energy not correct',BEXCEPTION_FATAL)
 endif

! write(*,*)testState%norm*sqrt(meter)

 ! test for correct norm. Multiply by sqrt(meter) to get dimensionless value.
 if(abs((testState%norm*sqrt(meter)-normExp)/normExp) > eps) then
  call throwException('class_WfSquareWell_run: Norm not correct',BEXCEPTION_FATAL)
 endif



! xx=xlo

! do ix=1,nx
!  write(25,*)xx,testState%getWavefn(xx)
!  xx=xx+(xhi-xlo)/nx
! enddo
 
end program class_WfSquareWell_run


