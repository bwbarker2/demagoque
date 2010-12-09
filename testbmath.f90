program testbmath
 implicit none

 integer, parameter :: Num=3
 complex*16 :: zdet2d
 complex*16, dimension(Num,Num) :: cmat
 integer, dimension(Num) :: L   !index vector for gaussian elimination

! write(*,*)cmplx(2,-3,4)
! write(*,*)zabs(cmplx(2.d0,-3.d0,8))

! cmat=0.d0
 cmat(1,1)=cmplx(1.d0,0.d0,8)
 cmat(1,2)=cmplx(2.d0,0.d0,8)
 cmat(1,3)=cmplx(3.1d0,0.d0,8)
 cmat(2,1)=cmplx(2.d0,0.d0,8)
 cmat(2,2)=cmplx(5.d0,0.d0,8)
 cmat(2,3)=cmplx(6.d0,0.d0,8)
 cmat(3,1)=cmplx(3.1d0,0.d0,8)
 cmat(3,2)=cmplx(6.d0,0.d0,8)
 cmat(3,3)=cmplx(8.d0,0.d0,8)

 write(*,*)cmat

! call zGauss(Num,cmat,L)

! write(*,*)cmat

 write(*,*)'det=',zdet2d(cmat,Num)

end program testbmath
