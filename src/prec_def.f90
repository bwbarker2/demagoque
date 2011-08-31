MODULE prec_def
 ! defines precision of types
  IMPLICIT NONE
  
  INTEGER,PARAMETER :: long=8  !just use double precision for now
!  integer,parameter :: long=selected_real_kind(15,307)

  integer,parameter :: stderr=102  !file unit of standard error
  
END MODULE prec_def
