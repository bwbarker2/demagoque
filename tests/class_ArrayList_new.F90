!> tests src/class_ArrayList.f90 to construct new ArrayList
program class_ArrayList_new
 use class_ArrayList

 type(iArrayList), dimension(10) :: myIArrayList
 type(dArrayList) :: myDArrayList

 myiArrayList(1) = new_iArrayList()

 write(*,*)'size:',myiArrayList(1)%size()
 write(*,*)'val(1):',myiArrayList(1)%add(1)


end program class_ArrayList_new
