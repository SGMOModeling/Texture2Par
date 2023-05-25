program Texture2Par
  use MakePar
  use fpath
  use IWFM
  use Krige
  use PilotPoints
  use errorhandle
  use str_utls
  use MFSupport
  implicit none

!-----------------------------------------------------------------------------!
!                             Texture2Par
!                         ********************
!                               V 1.0.0
!
! Converts well percent course with depth (1 or 0) data to aquifer parameters
! at specified nodes. Intended to work with input/output files of the
! Integrated Water Flow Model (IWFM) with added support for MODFLOW 2000/2005.
! Written by Leland Scantlebury and Marinko Karanovic of S.S. Papadopulos &
! Associates, based on original method and VBA program created by Timothy J.
! Durbin.
!
! Requires Input File: Texture2Par.in
!
! Copyright 2022 S.S. Papadopulos & Associates. All rights reserved.
!-----------------------------------------------------------------------------!

    integer                    :: ioerr, ilay, ifail, ipp, &
                                  i, ierr, iarr, twell, wpoint, zone, &
                                  modflowflag, zoneflag, geoflag
    real                       :: a_hmin
    character(150)             :: jnk,fpath_line
    character(256)             :: cerr
    character(1000)            :: line
    logical                    :: writenodefiles

! Intel-Fortran specific way of writing to same command line line
!    open(6, carriagecontrol='fortran')  ! Uses removed by VB

    ! Write Flag
    write(*,'(/,4x,a)') '             Texture2Par v1.0.0             '
    write(*,'(3x,a)') '----------------------------------------------'
    write(*,'(4x,a)')   'Copyright S.S. Papadopulos & Associates Inc.'

    ! Read input file
    ierr = 0
    writeNodeFiles = .false.
    write(*,'(/,2x,a)') 'Reading input file'
    open(10, file='Texture2Par.in', status='old', action='read', &
         iostat = ierr,  iomsg=cerr)
    call iomsghandler(ierr, cerr)

    ! Read in model type and well file
    call ReadToSymbol(10, '*')
    read(10,*) line
    read(10,*) wellfile
    read(10,*) geozonefile
    call ReadToSymbol(10, '*')

    ! Check if MODFLOW or IWFM
    if (index(to_upper(trim(line)),'IWFM')>0) then
      modflowflag = 0
    else if (index(to_upper(trim(line)), 'MODFLOW')>0) then
      modflowflag = 1
    else
      write(*,'(a)') 'Error - Invalid Model Type. Stopping.'
      write(*,'(a)') 'Valid Model Types: IWFM, MODFLOW'
      stop
    end if

    ! Different order/different files for each model
    if (modflowflag==0) then
      ! IWFM
      read(10,*) simFile
      read(10,*) fpath_line    !preSimFile
      call fpath_strip(fpath_line, prefolder, preSimFile)
      read(10,*) gwtemp
      read(10,*) zonefile
      ! MODFLOW Confining flag, set to 0 when IWFM
      ncbd = 0
    else if (modflowflag==1) then
      ! MODFLOW
      read(10,*) simFile
      read(10,*) layfiletemp
      read(10,*) zonefile
      read(10,*) xoff
      read(10,*) yoff
      read(10,*) mfrot
    end if

    ! Node output setting
    call ReadToSymbol(10, '*')
    read(10,*) line
    if (index(to_upper(trim(line)), 'TRUE')>0) writeNodeFiles = .true.

    ! Variogram Settings
    call ReadToSymbol(10, '*')
    read(10,*) itype
    read(10,*) sill
    read(10,*) range
    read(10,*) a_hmin
    read(10,*) ang1
    read(10,*) nugget
    read(10,*) pckrige_nwells

    ! Sediment Settings
    call ReadToSymbol(10, '*')
    read(10,*) KCk
    read(10,*) KFk
    read(10,*) KHp
    read(10,*) KVp
    read(10,*) Stp
    ! Handle p==0 - power-law function is discountinuous at 0
    if (KHp == 0.0) KHp = tiny(KHp)
    if (KVp == 0.0) KVp = tiny(KVp)

    ! Aquifer Pilot Point section
    ! TODO: Add descriptive error when not all pilot point parameters are present
    call ReadToSymbol(10, '*')
    call FindSectionLength(10, ppn, '*')
    call AquiferPPAllocate()
    do ipp = 1, ppn
      read(10,*) ppx(ipp), ppy(ipp), ppKCMin(ipp), ppdKC(ipp), ppKFMin(ipp), &
                 ppdKF(ipp), ppSsC(ipp), ppSsF(ipp), ppSyC(ipp),ppSyF(ipp), &
                 ppAnisoC(ipp), ppAnisoF(ipp), ppZone(ipp)
    end do

    ! Aquitard Pilot Point section
    !read(10,*)    !call ReadToSymbol(10, '*')
    call ReadToSymbol(10, '*')
    call FindSectionLength(10, atppn, '*')
    call AquitardPPAllocate()
    do ipp = 1, atppn
      read(10,*) ppatx(ipp), ppaty(ipp), ppatKCMin(ipp), ppatdKC(ipp), &
                 ppatKFMin(ipp), ppatdKF(ipp), ppatAnisoC(ipp), &
                 ppatAnisoF(ipp), ppatZone(ipp)
    end do
    close(10)

    ! MODEL SPECIFIC
    if (modflowflag == 0) then
      write(*,'(2x,a)') 'Reading IWFM Input Files'
      call readIWFM()
      call create_nodetree()
    else if (modflowflag == 1) then
      write(*,'(2x,a)') 'Reading MODFLOW Input Files'
      call read_modflow(simfile)
      call calcCellCenters_global()
      call create_nodetree()
    end if

    ! Pre-read Well file for file length, nwells, maxWPoints
    open(10, file = trim(wellfile), status='old', iostat = ierr,  iomsg=cerr)
    call iomsghandler(ierr, cerr)
    read(10,*)                                    ! Header
    ioerr = 0
    wfilelen = 0
    nwell = 0
    maxWPoints = 0
    do while (ioerr == 0)
      read(10,*, iostat = ioerr) jnk, twell, wpoint
      wfilelen = wfilelen + 1
      if (twell > nwell) nwell = twell
      if (wpoint > maxWPoints) maxWPoints = wpoint
    end do
    wfilelen = wfilelen - 1
    close(10)

    ! Read Node Zone file, unless zonefile was set to NONE
    allocate(pp_lay2arr(nlayers+natlayers))
    if (to_upper(trim(zonefile))=='NONE') then
      zoneflag = 0
    else
      ! External file
      zoneflag = 1
      if (modflowflag == 0) then
        call inputZones_IWFM()
      else
        call inputZones_modflow()
      end if
    end if

    ! Handle No Pilot Point Zones (same as one zone)
    if (zoneflag == 0) then
      call no_ppzone_setup(nnodes, nlayers, natlayers)
    end if

    ! Allocate now that many dimensions are known
    call WellAllocate()
    call NodeAllocate()

    ! Read geologic units file
    if (to_upper(trim(geozonefile))=='NONE') then
      ! No geologic zones
      ! perhaps a different subroutine should be used instead of dummy values?
      geoflag = 0
      ngeozones = 1
      ngeozones_lay(1:nlayers) = 1
      geozones(1:nnodes,1:nlayers) = '[NONE]'
      geonames_lay(1,1:nlayers) = '[NONE]'
    else
      geoflag = 1
      if (modflowflag == 0) then
        call readNodeGeoZones_IWFM()
      else if (modflowflag == 1) then
        call readNodeGeoZones_modflow()
      end if
    end if

    ! Write Read Summary
    write(*,'(2x,a)') 'Setup Summary:                          '
    write(*,'(4x,1a36,i8)') 'Number of aquifer pilot points: ', ppn
    write(*,'(4x,1a36,i8)') 'Number of aquitard pilot points: ', atppn
    write(*, '(4x,1a36,i8)') 'Number of Layers: ', nlayers
    if (modflowflag == 0) then
      write(*, '(4x,1a36,i8)') 'Number of Nodes: ', nnodes
      write(*, '(4x,1a36,i8)') 'Number of Elements: ', nelements
    else if (modflowflag == 1) then
      ! Perhaps "cells per layer"
      write(*, '(4x,1a36,i8)') 'Number of Cells: ', nnodes
      write(*, '(4x,1a36,i8)') 'Quasi-3D Confining Beds: ', ncbd
    end if
    !write(*,'(4x,1a36,i)') 'Well File Length: ', wfilelen
    write(*,'(4x,1a36,i8)') 'Number of Well Logs: ', nwell
    write(*,'(4x,1a36,i8)') 'Number of Pilot Point Zones: ', maxval(pp_arr_nzones)
    write(*,'(4x,1a36,i8,/)') 'Number of Geologic Zones: ', ngeozones

! Set anistropy and rotation matrix
    call setanis(a_hmin, a_hmin)
    call setrot()
! MODEL SPECIFIC
! Call input wells subroutine
    if (modflowflag==0) then
      call inputwells_IWFM(geoflag)
    else if (modflowflag==1) then
      call inputwells_modflow(geoflag)
    end if

! Call Layer Elevation subroutine
    write(*,'(2x,a)') 'Calculating Layer Elevations at Wells'
    if (modflowflag == 0) then
      call layerelevation(TopElev, BotElev, wellelevtop, wellelevbot)
      ! For aquitard
      call layerelevation(AqTardTopElev, AqTardBotElev, &
                          wellaqtardelevtop, wellaqtardelevbot)
    end if
    ! Modflow Well elevations are assigned during well read

! Move co-located wells (if kriging still crashes, well locations need to be revised)
    call colocate_corrector(xwell, ywell, nwell)

! Call Layer Texture subroutine
    write(*,'(2x,a)') 'Calculating Percent Coarse at Wells      '
    outfile = 'PcWellElem.out'
    call layertexture(wellelevtop, wellelevbot, PcWellElem)
    ! If IWFM or if there are MODFLOW confining beds
    if ((modflowflag == 0).or.(ncbd > 0 )) then
      outfile = 'PcWellAqTardElem.out'
      call layertexture(wellaqtardelevtop, wellaqtardelevbot, PcWellAqtardElem)
    end if

    ! Krige from pilot points to nodes
    ! Aquifer Parameters - KCMin, DeltaKC, KFMin, DeltaKF, SsC, SsF, SyC, Aniso
    write(*,'(2x,a,i3)') 'Kriging Aquifer Pilot Point Parameters'
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppKCMin,  nKCMin )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppdKC,    dKC    )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppKFMin,  nKFMin )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppdKF,    dKF    )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSsC,    nSsC   )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSsF,    nSsF   )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSyC,    nSyC   )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppSyF,    nSyF   )
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppAnisoC, nAnisoC)
    call ppaq_parkrige_arr(nnodes, nodex, nodey, ppAnisoF, nAnisoF)

    ! Krige Aquitard parameters from pilot points to nodes
    if ((modflowflag == 0).or.(ncbd > 0 )) then
      if (atppn > 0) then
        ! Aquitard Parameters - KCMin, DeltaKC, KFMin, DeltaKF, Aniso
        write(*,'(2x,a,i3)') 'Kriging Aquitard Pilot Point Parameters'
        call ppat_parkrige_arr(nnodes, nodex, nodey, ppatKCMin,  nKCMin )
        call ppat_parkrige_arr(nnodes, nodex, nodey, ppatdKC,    dKC    )
        call ppat_parkrige_arr(nnodes, nodex, nodey, ppatKFMin,  nKFMin )
        call ppat_parkrige_arr(nnodes, nodex, nodey, ppatdKF,    dKF    )
        call ppat_parkrige_arr(nnodes, nodex, nodey, ppatAnisoC, nAnisoC)
        call ppat_parkrige_arr(nnodes, nodex, nodey, ppatAnisoF, nAnisoF)
      else
        write(*,'(2x,a,/)') 'No Aquitard Pilot Points - Skipping       '
        nKCMin(:,ppaq_narr+1:pp_narr) = -999
        dKC   (:,ppaq_narr+1:pp_narr) = -999
        nKFMin(:,ppaq_narr+1:pp_narr) = -999
        dKF   (:,ppaq_narr+1:pp_narr) = -999
      end if
    end if

! Call interpolate subroutines (Well PC to Node Aquifer parameter)
    write(*,'(2x,a)') 'Calculating Node Aquifer Parameters           '
    call interpolate()
    if ((modflowflag == 0).or.(ncbd > 0 )) then
      call interpolate_aqtard()
    end if

! Write Output files, routines called depend on model type
    if (modflowflag == 0) then
      if (writeNodeFiles) then
        ! Only calculate if writing the file
        call calc_interbed_thickness()
        write(*,'(2x,a)') 'Writing Node Parameter Files'
        call writeWellFile_IWFM(PcWellElem,       't2p_WellPC.out')
        call writeWellFile_IWFM(PcWellAqtardElem, 't2p_WellAqTardPC.out')
        call writeNodeFile(PcNode, 't2p_NodePC.out')
        call writeNodeFile(KhB,    't2p_KhB.out'   )
        call writeNodeFile(KvB,    't2p_KvB.out'   )
        call writeNodeFile(SsB,    't2p_SsB.out'   )
        call writeNodeFile(SyB,    't2p_SyB.out'   )
        call writeNodeFile(interbed_thick,'t2p_InterbedThickness.out')
        call writeNodeFile(PcNode_aqtard, 't2p_NodeAqTardPC.out')
        call writeNodeFile(KvB_aqtard,    't2p_AqTardKvB.out'   )
      else
        write(*,'(2x,a)') 'WriteNodeFiles set to FALSE - no parameter files written'
      end if
        call writeIWFMgwfile()

    else if (modflowflag == 1) then
      if (writeNodeFiles) then
        ! Only calculate if writing the file
        call calc_interbed_thickness()
        write(*,'(2x,a)') 'Writing Cell Parameter Files'
        call writeWellFile_modflow(PcWellElem,       't2p_WellPC.out')
        call writeCellFile(PcNode, 't2p_NodePC.out')
        call writeCellFile(KhB,    't2p_KhB.out'   )
        call writeCellFile(KvB,    't2p_KvB.out'   )
        call writeCellFile(SsB,    't2p_SSB.out'   )
        call writeCellFile(SyB,    't2p_SyB.out'   )
        call writeCellFile(interbed_thick,'t2p_InterbedThickness.out')
        if (ncbd > 0 ) then
          call writeWellFile_modflow(PcWellAqtardElem, 't2p_WellAqTardPC.out')
          call writeCellFile(PcNode_aqtard, 't2p_AqTardNodePC.out')
          call writeCellFile(KvB_aqtard,    't2p_AqTardKvB.out'   )
        end if
      else
        write(*,'(2x,a)') 'WriteCellFiles set to FALSE - no parameter files written'
      end if
      call readwriteFlowPackage()
      !call clear_modflow_memory()  ! Unnecessary, crashes due to unallocated MF arrays
    end if

    write(*,'(/,2x,a)') 'Program finished!'

  end program Texture2Par
!-----------------------------------------------------------------------------!