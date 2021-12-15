
subroutine k2coll_collimate(icoll, iturn, ie, c_length, c_rotation, c_aperture, c_offset, c_tilt,  &
  x_in, xp_in, y_in, yp_in, p_in, s_in, enom, lhit_pos, lhit_turn, part_abs_pos_local,             &
  part_abs_turn_local, impact, indiv, lint, onesided, nhit_stage, j_slices, nabs_type, linside)

  use, intrinsic :: iso_fortran_env, only : int16
  use parpro
  use crcoall
  use coll_db
  use coll_common
  use coll_crystal, only : cry_doCrystal
  use coll_materials
  use mod_common, only : iexact, napx, unit208
  use mod_common_main, only : partID, naa
  use mathlib_bouncer
  use mod_ranlux


  integer,          intent(in)    :: icoll        ! Collimator ID
  integer,          intent(in)    :: iturn        ! Turn number
  integer,          intent(in)    :: ie           ! Structure element index

  real(kind=fPrec), intent(in)    :: c_length     ! Collimator length in m
  real(kind=fPrec), intent(in)    :: c_rotation   ! Collimator rotation angle vs vertical in radians
  real(kind=fPrec), intent(in)    :: c_aperture   ! Collimator aperture in m
  real(kind=fPrec), intent(in)    :: c_offset     ! Collimator offset in m
  real(kind=fPrec), intent(inout) :: c_tilt(2)    ! Collimator tilt in radians

  real(kind=fPrec), intent(inout) :: x_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: xp_in(npart) ! Particle coordinate
  real(kind=fPrec), intent(inout) :: y_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: yp_in(npart) ! Particle coordinate
  real(kind=fPrec), intent(inout) :: p_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: s_in(npart)  ! Particle coordinate

  real(kind=fPrec), intent(in)    :: enom         ! Reference momentum in GeV
  logical,          intent(in)    :: onesided

  integer,          intent(inout) :: lhit_pos(npart)
  integer,          intent(inout) :: lhit_turn(npart)
  integer,          intent(inout) :: part_abs_pos_local(npart)
  integer,          intent(inout) :: part_abs_turn_local(npart)
  integer,          intent(inout) :: nabs_type(npart)
  integer,          intent(inout) :: nhit_stage(npart)
  real(kind=fPrec), intent(inout) :: indiv(npart)
  real(kind=fPrec), intent(inout) :: lint(npart)
  real(kind=fPrec), intent(inout) :: impact(npart)
  logical,          intent(inout) :: linside(napx)

  logical isImp
  integer j,nabs,nhit,j_slices
  real(kind=fPrec) keeps,fracab,drift_length,mirror,tiltangle
  real(kind=fPrec) x00,z00,p,sp,s,s_impact
  real(kind=fPrec) x_flk,xp_flk,y_flk,yp_flk,s_flk,zpj
  real(kind=fPrec) x_Dump,xpDump,y_Dump,ypDump,s_Dump
  real(kind=fPrec) cRot,sRot,cRRot,sRRot
  real(kind=fPrec) xIn,xpIn,yIn,ypIn,xOut,xpOut,yOut,ypOut,sImp,sOut
  real(kind=fPrec) x_in0,xp_in0

  ! ien0,ien1: particle energy entering/leaving the collimator
  ! energy in MeV
  real(kind=fPrec)    :: ien0, ien1
  integer(kind=int16) :: nnuc0,nnuc1

  ! Initilaisation
  mat    = cdb_cMaterialID(icoll)
  length = c_length
  p0     = enom

  ! Initialise scattering processes
  call k2coll_scatin(p0)

  nhit   = 0
  nabs   = 0
  fracab = zero
  mirror = one

  ! Compute rotation factors for collimator rotation
  cRot   = cos_mb(c_rotation)
  sRot   = sin_mb(c_rotation)
  cRRot  = cos_mb(-c_rotation)
  sRRot  = sin_mb(-c_rotation)

  !Set energy and nucleon change variables as with the coupling
  nnuc0 = 0
  ien0  = zero
  nnuc1 = 0
  ien1  = zero

  do j=1,napx

    if(part_abs_pos_local(j) /= 0 .and. part_abs_turn_local(j) /= 0) then
      ! Don't do scattering process for particles already absorbed
      cycle
    end if

    impact(j) = -one
    lint(j)   = -one
    indiv(j)  = -one

    x      = x_in(j)
    xp     = xp_in(j)
    xp_in0 = xp_in(j)
    z      = y_in(j)
    zp     = yp_in(j)
    p      = p_in(j)
    sp     = zero
    dpop   = (p - p0)/p0
    x_flk  = zero
    y_flk  = zero
    xp_flk = zero
    yp_flk = zero

    ! Transform particle coordinates to get into collimator coordinate  system
    ! First do rotation into collimator frame
    x  =  x_in(j)*cRot + sRot*y_in(j)
    z  =  y_in(j)*cRot - sRot*x_in(j)
    xp = xp_in(j)*cRot + sRot*yp_in(j)
    zp = yp_in(j)*cRot - sRot*xp_in(j)

    ! For one-sided collimators consider only positive X. For negative X jump to the next particle
    if(onesided .and. x < zero .and. (icoll /= ipencil .or. iturn /= 1)) then
      cycle
    end if

! Log input energy + nucleons as per the FLUKA coupling
    nnuc0   = nnuc0 + naa(j)
    ien0    = ien0 + rcp(j) * c1e3


    ! Now mirror at the horizontal axis for negative X offset
    if(x < zero) then
      mirror    = -one
      tiltangle = -one*c_tilt(2)
    else
      mirror    = one
      tiltangle = c_tilt(1)
    end if
    x  = mirror*x
    xp = mirror*xp

    ! Shift with opening and offset
    x = (x - c_aperture/2D0) - mirror*c_offset

    ! Include collimator tilt
    if(tiltangle > zero) then
      xp = xp - tiltangle
    end if
    if(tiltangle < zero) then
      x  = x + sin_mb(tiltangle) * c_length
      xp = xp - tiltangle
    end if

    ! CRY Only: x_in0 has to be assigned after the change of reference frame
    x_in0 = x

    ! For selected collimator, first turn reset particle distribution to simple pencil beam
    if(icoll == ipencil .and. iturn == 1 .and. pencil_distr /= 3) then
      ! TW why did I set this to 0, seems to be needed for getting
      !    right amplitude => no "tilt" of jaw for the first turn !!!!
      c_tilt(1) = zero
      c_tilt(2) = zero

      ! Standard pencil beam as implemented by GRD
      if(pencil_rmsx == zero .and. pencil_rmsy == zero) then
        x  = pencil_dx(icoll)
        xp = zero
        z  = zero
        zp = zero
      end if

      ! Rectangular (pencil-beam) sheet-beam with
      ! pencil_offset is the rectangulars center
      ! pencil_rmsx defines spread of impact parameter
      ! pencil_rmsy defines spread parallel to jaw surface
      if(pencil_distr == 0 .and. (pencil_rmsx /= zero .or. pencil_rmsy /= zero)) then
        ! how to assure that all generated particles are on the jaw ?!
        x  = pencil_dx(icoll) + pencil_rmsx*(coll_rand() - half)
        xp = zero
        z  = pencil_rmsy*(coll_rand() - half)
        zp = zero
      end if

      ! Gaussian (pencil-beam) sheet-beam with ------- TW
      ! pencil_offset is the mean  gaussian distribution
      ! pencil_rmsx defines spread of impact parameter
      ! pencil_rmsy defines spread parallel to jaw surface
      if(pencil_distr == 1 .and. (pencil_rmsx /= zero .or. pencil_rmsy /= zero)) then
        x  = pencil_dx(icoll) + pencil_rmsx*ran_gauss(2D0)
        ! all generated particles are on the jaw now
        x  = sqrt(x**2)
        xp = zero
        z  = pencil_rmsy*ran_gauss(2D0)
        zp = zero
      end if

      ! Gaussian (pencil-beam) sheet-beam with ------- TW
      ! pencil_offset is the mean  gaussian distribution
      ! pencil_rmsx defines spread of impact parameter
      !             here pencil_rmsx is not gaussian!!!
      ! pencil_rmsy defines spread parallel to jaw surface
      if(pencil_distr == 2 .and. (pencil_rmsx /= zero .or. pencil_rmsy /= zero)) then
        x  = pencil_dx(icoll) + pencil_rmsx*(coll_rand() - half)
        ! all generated particles are on the jaw now
        x  = sqrt(x**2)
        xp = zero
        z  = pencil_rmsy*ran_gauss(2D0)
        zp = zero
      end if

      ! Selection of pos./neg. jaw  implemented by GRD ---- TW
      ! ensure that for onesided only particles on pos. jaw are created
      if(onesided) then
        mirror = one
      else
        if(coll_rand() < half) then
          mirror = -one
        else
          mirror = one
        end if
      end if

      ! TW if c_tilt is set to zero before entering pencil beam
      !    section the assigning of the tilt will result in assigning zeros
      if(mirror < zero) then
        tiltangle = c_tilt(2)
      else
        tiltangle = c_tilt(1)
      end if

      write(coll_pencilUnit,"(f10.8,4(2x,f10.8))") x, xp, z, zp, tiltangle
      flush(coll_pencilUnit)



    end if ! End pencil dist

    ! After finishing the coordinate transformation, or the coordinate manipulations in case of pencil beams,
    ! save the initial coordinates of the impacting particles
    xIn  = x
    xpIn = xp
    yIn  = z
    ypIn = zp

    ! particle passing above the jaw are discarded => take new event
    ! entering by the face, shorten the length (zlm) and keep track of
    ! entrance longitudinal coordinate (keeps) for histograms

    ! The definition is that the collimator jaw is at x>=0.

    ! 1) Check whether particle hits the collimator
    isImp = .false.
    s     = zero
    keeps = zero
    zlm   = -one*length

    if(cdb_isCrystal(icoll)) then ! This is a crystal collimator

      call cry_doCrystal(ie,iturn,j,mat,x,xp,z,zp,s,p,x_in0,xp_in0,zlm,sImp,isImp,nhit,nabs,lhit_pos,lhit_turn,&
        part_abs_pos_local,part_abs_turn_local,impact,indiv,c_length)

      if(nabs /= 0) then
        part_abs_pos_local(j)  = ie
        part_abs_turn_local(j) = iturn
        lint(j)                = zlm
      end if

      sImp  = (s - c_length) + sImp
      sOut  = s
      xOut  = x
      xpOut = xp
      yOut  = z
      ypOut = zp

    else ! "Normal" collimator

      if(x >= zero) then
        ! Particle hits collimator and we assume interaction length ZLM equal
        ! to collimator length (what if it would leave collimator after
        ! small length due to angle???)
        zlm       = length
        impact(j) = x
        indiv(j)  = xp
      else if(xp <= zero) then
        ! Particle does not hit collimator. Interaction length ZLM is zero.
        zlm = zero
      else
        ! Calculate s-coordinate of interaction point
        s = (-one*x)/xp
        if(s <= zero) then
          write(lerr,"(a)") "COLLK2> ERROR S <= zero. This should not happen!"
          call prror
        end if
        if(s < length) then
          zlm       = length - s
          impact(j) = zero
          indiv(j)  = xp
        else
          zlm = zero
        end if
      end if

      ! First do the drift part
      ! DRIFT PART
      drift_length = length - zlm
      if(drift_length > zero) then
        if(iexact) then
          zpj = sqrt(one-xp**2-zp**2)
          x   = x  + drift_length*(xp/zpj)
          z   = z  + drift_length*(zp/zpj)
          sp  = sp + drift_length
        else
          x  = x  + xp* drift_length
          z  = z  + zp * drift_length
          sp = sp + drift_length
        end if
      end if

      ! Now do the scattering part
      if(zlm > zero) then
        if(.not.linside(j)) then
          ! first time particle hits collimator: entering jaw
          linside(j) = .true.
          if(dowrite_impact) then
            if(tiltangle > zero) then
              x_Dump = (x + c_aperture/2D0 + tiltangle*sp)*mirror + c_offset
            else
              x_Dump = (x + c_aperture/2D0 + tiltangle*(sp - c_length))*mirror + c_offset
            end if
            xpDump = (xp + tiltangle)*mirror
            y_Dump = z
            ypDump = zp
            s_Dump = sp+real(j_slices-1,fPrec)*c_length
            write(coll_jawProfileUnit,"(3(1x,i7),5(1x,e17.9),1x,i1)") &
              icoll,iturn,partID(j),x_Dump,xpDump,y_Dump,ypDump,s_Dump,1
            flush(coll_jawProfileUnit)



          end if
        end if

        s_impact = sp
        nhit = nhit + 1
        call k2coll_jaw(s,nabs,icoll,iturn,partID(j))

        nabs_type(j) = nabs
        lhit_pos(j)  = ie
        lhit_turn(j) = iturn

        isImp = .true.
        sImp  = s_impact+(real(j_slices,fPrec)-one)*c_length
        sOut  = (s+sp)+(real(j_slices,fPrec)-one)*c_length
        xOut  = x
        xpOut = xp
        yOut  = z
        ypOut = zp

        ! Writeout should be done for both inelastic and single diffractive. doing all transformations
        ! in x_flk and making the set to 99.99 mm conditional for nabs=1
        if(dowrite_impact .or. nabs == 1 .or. nabs == 4) then
          ! Transform back to lab system for writeout.
          ! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

          x_flk  = xInt
          xp_flk = xpInt

          if(tiltangle > zero) then
            x_flk  = x_flk  + tiltangle*(sInt+sp)
            xp_flk = xp_flk + tiltangle
          else if(tiltangle < zero) then
            xp_flk = xp_flk + tiltangle
            x_flk  = x_flk  - sin_mb(tiltangle) * (length-(sInt+sp))
          end if

          x_flk  = (x_flk + c_aperture/2D0) + mirror*c_offset
          x_flk  = mirror*x_flk
          xp_flk = mirror*xp_flk
          y_flk  = (  yInt*cRRot -  x_flk*sRRot)*c1e3
          yp_flk = ( ypInt*cRRot - xp_flk*sRRot)*c1e3
          x_flk  = ( x_flk*cRRot +   yInt*sRRot)*c1e3
          xp_flk = (xp_flk*cRRot +  ypInt*sRRot)*c1e3
          s_flk  = (sInt+sp)+(real(j_slices,fPrec)-one)*c_length

          if(dowrite_impact) then
            ! Write out all impacts to all_impacts.dat
            write(coll_flukImpAllUnit,"(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))") &
              icoll,c_rotation,s_flk,x_flk,xp_flk,y_flk,yp_flk,nabs,partID(j),iturn
            flush(coll_flukImpAllUnit)



            if(nabs == 1 .or. nabs == 4) then
              ! Standard FLUKA_impacts writeout of inelastic and single diffractive
              write(coll_flukImpUnit,"(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))") &
                icoll,c_rotation,s_flk,x_flk,xp_flk,y_flk,yp_flk,nabs,partID(j),iturn
              flush(coll_flukImpUnit)



            end if
          end if

          ! Finally, the actual coordinate change to 99 mm
          if(nabs == 1) then
            fracab  = fracab + 1
            x       = 99.99e-3_fPrec
            z       = 99.99e-3_fPrec
            lint(j) = zlm
            part_abs_pos_local(j)  = ie
            part_abs_turn_local(j) = iturn
          end if
        end if
      end if ! Collimator jaw interaction

      if(nabs /= 1 .and. zlm > zero) then
        ! Do the rest drift, if particle left collimator early
        drift_length = (length-(s+sp))
        if(drift_length > c1m15) then
          linside(j) = .false.
          if(dowrite_impact) then
            if(tiltangle > zero) then
              x_Dump = (x + c_aperture/2D0 + tiltangle*(s+sp))*mirror + c_offset
            else
              x_Dump = (x + c_aperture/2D0 + tiltangle*(s+sp-c_length))*mirror + c_offset
            end if
            xpDump = (xp+tiltangle)*mirror
            y_Dump = z
            ypDump = zp
            s_Dump = s+sp+real(j_slices-1,fPrec)*c_length
            write(coll_jawProfileUnit,"(3(1x,i7),5(1x,e17.9),1x,i1)") &
              icoll,iturn,partID(j),x_Dump,xpDump,y_Dump,ypDump,s_Dump,2
            flush(coll_jawProfileUnit)



          end if
          if(iexact) then
            zpj = sqrt(one-xp**2-zp**2)
            x   = x  + drift_length*(xp/zpj)
            z   = z  + drift_length*(zp/zpj)
            sp  = sp + drift_length
          else
            x  = x  + xp * drift_length
            z  = z  + zp * drift_length
            sp = sp + drift_length
          end if
        end if
        lint(j) = zlm - drift_length
      end if

    end if ! Collimator isCrystal

    if(dowrite_impact .and. isImp .and. nhit_stage(j) == 0) then
      write(coll_fstImpactUnit,"(i8,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))") &
        partID(j),iTurn,iColl,nAbs,sImp,sOut,xIn,xpIn,yIn,ypIn,xOut,xpOut,yOut,ypOut
      flush(coll_fstImpactUnit)
    end if

    ! Transform back to particle coordinates with opening and offset
    if(x < 99.0e-3_fPrec) then
      ! Include collimator tilt
      if(tiltangle > zero) then
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
      else if(tiltangle < zero) then
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
        x  = x  - sin_mb(tiltangle) * c_length
      end if

      ! Transform back to particle coordinates with opening and offset
      z00 = z
      x00 = x + mirror*c_offset
      x   = (x + c_aperture/2D0) + mirror*c_offset

      ! Now mirror at the horizontal axis for negative X offset
      x  = mirror * x
      xp = mirror * xp

      ! Last do rotation into collimator frame
      x_in(j)  =  x*cRRot +  z*sRRot
      y_in(j)  =  z*cRRot -  x*sRRot
      xp_in(j) = xp*cRRot + zp*sRRot
      yp_in(j) = zp*cRRot - xp*sRRot

! Log output energy + nucleons as per the FLUKA coupling
! Do not log dead particles
      nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
      ien1        = ien1  + rcp(j) * c1e3                   ! outcoming energy

      if(icoll == ipencil .and. iturn == 1 .and. pencil_distr /= 3) then
        x00      = mirror * x00
        x_in(j)  = x00*cRRot + z00*sRRot
        y_in(j)  = z00*cRRot - x00*sRRot
        xp_in(j) = xp_in(j) + mirror*xp_pencil0
        yp_in(j) = yp_in(j) + mirror*yp_pencil0
        x_in(j)  = x_in(j)  + mirror*x_pencil(icoll)
        y_in(j)  = y_in(j)  + mirror*y_pencil(icoll)
      end if

      if(cdb_isCrystal(icoll)) then
        p_in(j) = p
        s_in(j) = s_in(j) + s
      else
        p_in(j) = (one + dpop) * p0
        s_in(j) = sp + (real(j_slices,fPrec)-one) * c_length
      end if
    else
      x_in(j) = x
      y_in(j) = z
    end if

  end do ! End of loop over all particles

! write out energy change over this collimator
  if((ien0-ien1) > one) then
    write(unit208,"(2(i6,1x),e24.16)") icoll, (nnuc0-nnuc1), c1m3*(ien0-ien1)
    flush(unit208)
  end if

end subroutine k2coll_collimate