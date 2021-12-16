program main
  use floatPrecision
  use numerical_constants
  use coll_k2


  integer(kind=4)  :: icoll
  integer(kind=4)  :: iturn
  integer(kind=4)  :: ie
  real(kind=fPrec) :: c_length
  real(kind=fPrec) :: c_rotation
  real(kind=fPrec) :: c_aperture
  real(kind=fPrec) :: c_offset
  real(kind=fPrec) :: c_tilt(2)
  real(kind=fPrec) :: rcx(20000)
  real(kind=fPrec) :: rcxp(20000)
  real(kind=fPrec) :: rcy(20000)
  real(kind=fPrec) :: rcyp(20000)
  real(kind=fPrec) :: rcp(20000)
  real(kind=fPrec) :: rcs(20000)
  real(kind=fPrec) :: c_enom
  integer(kind=4)  :: part_hit_pos(20000)
  integer(kind=4)  :: part_hit_turn(20000)
  integer(kind=4)  :: part_abs_pos(20000)
  integer(kind=4)  :: part_abs_turn(20000)
  real(kind=fPrec) :: part_impact(20000)
  real(kind=fPrec) :: part_indiv(20000)
  real(kind=fPrec) :: part_linteract(20000)
  logical(kind=4)  :: onesided
  integer(kind=4)  :: nhit_stage(20000)
  integer(kind=4)  :: nabs_type(20000)
  logical(kind=4)  :: linside(20000)


  icoll = 31
  iturn = 1
  ie =1
  c_length = 0.59999999999999998
  c_rotation = 0
  c_aperture = 0.0025711021962573095
  c_offset = 0
  c_tilt = (0, 0)
  c_enom = 7000000
  onesided = .FALSE.
  linside(:) = .FALSE.

!   OPEN(UNIT=10, FILE="rcx.dump", ACTION="read", FORM="unformatted")
!     READ(10) rcx
!   CLOSE(UNIT=10)
!   OPEN(UNIT=11, FILE="rcxp.dump", ACTION="read", FORM="unformatted")
!     READ(11) rcxp
!   CLOSE(UNIT=11)
!   OPEN(UNIT=12, FILE="rcy.dump", ACTION="read", FORM="unformatted")
!     READ(12) rcy
!   CLOSE(UNIT=12)
!   OPEN(UNIT=13, FILE="rcyp.dump", ACTION="read", FORM="unformatted")
!     READ(13) rcyp
!   CLOSE(UNIT=13)
!   OPEN(UNIT=14, FILE="rcp.dump", ACTION="read", FORM="unformatted")
!     READ(14) rcp
!   CLOSE(UNIT=14)
!   OPEN(UNIT=15, FILE="rcs.dump", ACTION="read", FORM="unformatted")
!     READ(15) rcs
!   CLOSE(UNIT=15)
!   OPEN(UNIT=16, FILE="part_indiv.dump", ACTION="read", FORM="unformatted")
!     READ(16) part_indiv
!   CLOSE(UNIT=16)
!   OPEN(UNIT=17, FILE="part_hit_pos.dump", ACTION="read", FORM="unformatted")
!     READ(17) part_hit_pos
!   CLOSE(UNIT=17)
!   OPEN(UNIT=18, FILE="part_hit_turn.dump", ACTION="read", FORM="unformatted")
!     READ(18) part_hit_turn
!   CLOSE(UNIT=18)
!   OPEN(UNIT=19, FILE="part_abs_pos.dump", ACTION="read", FORM="unformatted")
!     READ(19) part_abs_pos
!   CLOSE(UNIT=19)
!   OPEN(UNIT=20, FILE="part_abs_turn.dump", ACTION="read", FORM="unformatted")
!     READ(20) part_abs_turn
!   CLOSE(UNIT=20)
!   OPEN(UNIT=21, FILE="part_impact.dump", ACTION="read", FORM="unformatted")
!     READ(21) part_impact
!   CLOSE(UNIT=21)
!   OPEN(UNIT=22, FILE="part_linteract.dump", ACTION="read", FORM="unformatted")
!     READ(22) part_linteract
!   CLOSE(UNIT=22)
!   OPEN(UNIT=23, FILE="nhit_stage.dump", ACTION="read", FORM="unformatted")
!     READ(23) nhit_stage
!   CLOSE(UNIT=23)
!   OPEN(UNIT=24, FILE="nabs_type.dump", ACTION="read", FORM="unformatted")
!     READ(24) nabs_type
!   CLOSE(UNIT=24)

! ! print *, c1m3
! print *, rcx
! print *, rcyp

!  k2coll_collimate(icoll, iturn, ie, c_length, c_rotation, c_aperture, c_offset, &
!      c_tilt, rcx, rcxp, rcy, rcyp, rcp, rcs, c_enom*c1m3, part_hit_pos, part_hit_turn, &
!      part_abs_pos, part_abs_turn, part_impact, part_indiv, part_linteract,             &
!      onesided, nhit_stage, 1, nabs_type, linside)



end program main

