MODULE IWFM
  use MakePar
  use errorhandle
  implicit none

!-----------------------------------------------------------------------------!
! IWFM Module for Texture2Par
! Combination of subroutines written by SSP&A and John Doherty
! used for the reading/writing of IWFM files
!
! Note on Nodes:
! - Nodes numbers ("NodeIDs") in IWFM input files can have gaps and might not
!   even be sequential (as of Oct 2020).
! - NodeIDs are what are read in from the input files (nid for short)
! - The sequential values used internally in T2P arrays are node indexes
!   (node for short, somewhat confusingly)
! - To get a nodeid from a node number, use the nodeid() array
! - To get a node index from a node id, use the nodebyid() array
!
! Note on Elements (12/27/2022)
! - Internally elements are numbered i:nelements within the array
!   elements(nelement,4[nodes])
! - Support non-contiguous (or otherwise) element IDs is provided through
!   an elemid(nelement) array in this module.
!
! Copyright 2022 S.S. Papadopulos & Associates. All rights reserved.
!-----------------------------------------------------------------------------!

  SAVE
  integer, parameter                  :: NUM_WORD_DIM=100
  integer, dimension(NUM_WORD_DIM)    :: left_word, right_word
  integer, allocatable                :: elemid(:)
  character(3000)                     :: cline
  character(200)                      :: nodefile, elementfile, prefolder, &
                                         preSimFile, simFile, stratfile, &
                                         gwfile, gwtemp

  CONTAINS

!-----------------------------------------------------------------------------!

subroutine readIWFM()
  use fpath
  implicit none

  integer                  :: i, k, ifail, ioerr, ierr, ipos, nreg, intd, &
                              idum, maxnodeid, nid, node
  integer                  :: ls(100), rs(100)
  real                     :: fact
  character(150)           :: jnk,currfile
  character(256)           :: cerr

  ! Read pre-processor file
!    write(*,'(2x,a)') 'Reading IWFM pre-processor file...'
    !currfile = adjustl(trim(prefolder)) // adjustl(trim(preSimFile))
    call fpath_join(prefolder, presimfile, currfile)
    open(12, file=trim(currfile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    do i = 1,4
      call readOneLine(preSimFile, 12, cline)
    end do
    read(12, *, iostat=ierr) elementfile
    read(12, *, iostat=ierr) nodefile
    read(12, *, iostat=ierr) stratfile
    close(12)

! Read node file
!    write(*,'(2x,a)') 'Reading IWFM node file...'
    !currfile = adjustl(trim(prefolder)) // adjustl(trim(nodefile))
    call fpath_join(prefolder, nodefile, currfile)
    open(12, file=trim(currfile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    call readOneLine(nodefile,12,cline)
    call multisplit(ifail,1,ls,rs,cline)
    call ifailhandler(ifail, nodefile, cline)
    call intread(ifail,cline(ls(1):rs(1)),nnodes)
    call ifailhandler(ifail, nodefile, cline)

    ! TODO : ALLOCATE NODE ROUTINE

    allocate(nodex(nnodes),nodey(nnodes), &
              nodeelem(nnodes, max_node_elem_members), &
              nnodeelem(nnodes),nodexy(nnodes,2), nodeid(nnodes))
    read(12,*) fact
    call readOneLine(elementfile,12,cline)
    BACKSPACE(12)

    maxnodeid = 0
    do i=1,nnodes
      read(12,*)intd,nodex(i),nodey(i)
      nodeid(i) = intd
      if (intd > maxnodeid) maxnodeid = intd
      nodex(i) = nodex(i) * fact
      nodey(i) = nodey(i) * fact
      ! Initialize node-element membership array to 0
      nnodeelem(i) = 0
    end do
    close(12)

    call setup_nodeindex_byid(maxnodeid)

! Read elements file
!    write(*,'(2x,a)') 'Reading IWFM elements file...'
    !currfile = adjustl(trim(prefolder)) // adjustl(trim(elementfile))
    call fpath_join(prefolder, elementfile, currfile)
    open(12, file=trim(currfile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    do
      call readOneLine(elementfile,12,cline)
      !cline= adjustl(trim(cline))
      IPOS=SCAN(cline,'NE')
      if (IPOS.GT.0) THEN
        EXIT
      end if
    end do

    call multisplit(ierr,1,ls,rs,cline)
    call intread(ierr,cline(ls(1):rs(1)),nelements)

    ! TODO : ALLOCATE ELEMENTS ROUTINE

    allocate(elements(nelements,4), elemid(nelements))

    call readOneLine(elementfile,12,cline)
    call multisplit(ifail,1,ls,rs,cline)
    call ifailhandler(ifail, elementfile, cline)
    call intread(ifail,cline(ls(1):rs(1)),nreg)
    call ifailhandler(ifail, elementfile, cline)
    call readOneLine(elementfile,12,cline)

    ! Subregions
    do i=1,nreg-1
      call readOneLine(elementfile,12,cline)
    end do

    ! Element array - element_id, node1, node2, node3, [node4], subregion
    call readOneLine(elementfile,12,cline)
    BACKSPACE(12)
    do i=1,nelements
      READ(12,*) elemid(i), (elements(i,k),k=1,4)
      ! Translate nodeid to node index for each connection LS 10/28/2020
      do k=1,4
        if (elements(i,k) > 0) then
          if (nodebyid(elements(i,k)) < 1) then
            write(*,*) 'ERROR! Invalid node ID in Element File. Element Number = ', elemid(i)  !elements(intd,k)
            stop
          end if
          elements(i,k) = nodebyid(elements(i,k))
        end if
      end do
      ! Add to node-element relationship array
      ! So we can identify by node what element we're in
      call add_elem2nodeelem(i)
    end do
    close(12)

! Read Stratigraphy file
!    write(*,'(2x,a)') 'Reading IWFM stratigraphy file'
    !currfile = adjustl(trim(prefolder)) // adjustl(trim(stratfile))
    call fpath_join(prefolder, stratfile, currfile)
    open(12, file=trim(currfile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    call readOneLine(nodefile,12,cline)
    call multisplit(ifail,1,ls,rs,cline)
    call ifailhandler(ifail, stratfile, cline)
    call intread(ifail,cline(ls(1):rs(1)),nlayers)
    ! In IWFM, every layer has an aquitard
    natlayers = nlayers
    call ifailhandler(ifail, stratfile, cline)
    allocate(Elevation    (nnodes),         &
             TopElev      (nnodes,nlayers), &
             BotElev      (nnodes,nlayers), &
             Thick        (nnodes,nlayers), &
             AqTardThick  (nnodes,nlayers), &
             AqTardTopElev(nnodes,nlayers), &
             AqTardBotElev(nnodes,nlayers)  )
    topelev = 0.0
    botelev = 0.0
    AqTardTopElev = 0.0
    AqTardBotElev = 0.0
    !call readOneLine(nodefile,12,cline)
    read(12, *) fact
    call readOneLine(nodefile,12,cline)
    BACKSPACE(12)
    do i=1,nnodes
      call readOneLine(nodefile,12,cline)
      read(cline,*) nid, Elevation(nodebyid(nid)),((AqTardThick(nodebyid(nid),k),Thick(nodebyid(nid),k)),k=1,nlayers)
      ! Adjust for Factor
      node = nodebyid(nid)
      if (node < 1) then
        write(*,*) 'ERROR! Invalid node ID in Stratigraphy File = ', elements(intd,k)
        stop
      end if
      Elevation(node)             = Elevation(node) * fact
      AqTardThick(node,1:nlayers) = AqTardThick(node,1:nlayers) * fact
      Thick(node,1:nlayers)       = Thick(node,1:nlayers)       * fact
      do k = 1,nlayers
        ! Aquitard layer is above aquifer layer
        if (k.eq.1) then
          AqTardTopElev(node,k) = Elevation(node)
          AqTardBotElev(node,k) = AqTardTopElev(node,k)-AqTardThick(node,k)
          TopElev(node,k) = Elevation(i)-AqTardThick(node,k)
          BotElev(node,k) = TopElev(node,k)-Thick(node,k)
        else
          AqTardTopElev(node,k) = BotElev(node,k-1)
          AqTardBotElev(node,k) = AqTardTopElev(node,k) - AqTardThick(node,k)
          TopElev(node,k) = BotElev(node,k-1)-AqTardThick(node,k)
          BotElev(node,k) = TopElev(node,k)-Thick(node,k)
        endif
      end do
    end do
    close(12)

! Read simulation file
!    write(*,'(2x,a)') 'Reading simulation file'
    currfile = adjustl(trim(simfile))
    open(12, file=trim(currfile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    do i= 1,4
      call readOneLine(simfile,12,cline)
    end do
    read(12,*,iostat=ierr) GWfile
    close(12)

end subroutine readIWFM

!-----------------------------------------------------------------------------!
SUBROUTINE inputwells_IWFM(readzones)
    use MakePar
    use kd_tree
    use str_utls
    use errorhandle
    implicit none

!-----------------------------------------------------------------------------!
! Inputs well information.
! Uses NodeTree to locate closest nodes and pass this information to
! subroutine that locates cell/element of well.
!
! Well GeoZones added 7/1/2019 -LS
!   Argument "readzones" is a flag for whether GeoZones exist (1 = yes, 0 = no)
!-----------------------------------------------------------------------------!

    integer, parameter       :: nnear=5
    integer                  :: i, iwell, izone, igeo, ilay, readzones, &
                                currwellelem, node, nwellgeozones, ierr
    integer                  :: closenodes(nnear)
    real                     :: depth, tx, ty, tz, td, tpc
    real                     :: nodedists(nnear), queryxy(2)
    character(30)            :: wellname, geobylayer(nlayers), &
                                wellgeonames(max_zones)
    character(256)           :: cerr
    logical                  :: newzone

    ! Initialize (to prevent issues where compiler does not initilize to 0)
    nwellgeozones = 0
    wellelem = NODATA

    ! Read well file
    write(*,'(2x,2a)') 'Reading Well File ', trim(wellfile)
    open(10, file =trim(wellfile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    read(10,*)                        ! Header

    if (readzones == 0) then
      ! No GeoZones
      nwellgeozones = 1
      wellgeozones(1:nwell,1:nlayers) = '[NONE]'
      do i=1, wfilelen
        read(10,*) wellname, iwell, izone, tpc, tx, ty, tz, td
        ! Well "zone" is discrete depth interval counter and should not be
        ! confused with other uses of zone (e.g. pilot points, geologic)
        PcWellZone(iwell, izone) = tpc
        numZone(iwell) = izone
        if (izone == 1) then
          Xwell(iwell) = tx
          Ywell(iwell) = ty
          Zland(iwell) = tz
          ! Use node tree to find closest node
          queryxy(1) = tx
          queryxy(2) = ty
          call n_nearest_to(nodetree, queryxy, nnear, closenodes, nodedists)
          wellelem(iwell) = get_well_element(tx, ty, iwell, closenodes, nnear)
        end if
        Zwell(iwell, izone) = tz - td
      end do
    else
      ! GeoZones
      do i=1, wfilelen
        read(10,*) wellname, iwell, izone, tpc, tx, ty, tz, td, geobylayer
        ! Well "zone" is discrete depth interval counter and should not be
        ! confused with other uses of zone (e.g. pilot points, geologic)
        PcWellZone(iwell, izone) = tpc
        numZone(iwell) = izone
        if (izone == 1) then
          Xwell(iwell) = tx
          Ywell(iwell) = ty
          Zland(iwell) = tz

          ! Well geologic zones
          wellgeozones(iwell,:) = geobylayer(:)
          ! Counts geozones for reporting and comparing
          do ilay=1, nlayers
            call strtracker(nwellgeozones, geobylayer(ilay), &
                  wellgeonames, newzone)
          end do

          ! Use node tree to find closest node
          queryxy(1) = tx
          queryxy(2) = ty
          call n_nearest_to(nodetree, queryxy, nnear, closenodes, nodedists)
          wellelem(iwell) = get_well_element(tx, ty, iwell, closenodes, nnear)

        end if
        Zwell(iwell, izone) = tz - td
      end do
    end if
    close(10)

    ! Report if mismatch in model vs well geo zone count
    if (nwellgeozones > ngeozones) then
      ! Just a warning - the extra wells will not be used
      write(*,'(4x,a)') 'Warning - more well geologic zones than node geologic zones'
      write(*,'(4x,i3,a,i3)') nwellgeozones, ' vs ', ngeozones
    else if (nwellgeozones < ngeozones) then
      ! All node geologic zones MUST have wells - otherwise can't calc PC
      write(*,'(4x,a)') 'Error - more node geologic zones than well geologic zones'
      write(*,'(4x,i3,a,i3)') ngeozones, ' vs ', nwellgeozones
      stop
    end if

    !write(*,'(4x,a)') 'Well file read successfully.'

end subroutine inputwells_IWFM
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine readNodeGeoZones_IWFM()
  use MakePar
  use PilotPoints
  use str_utls
  implicit none
!-----------------------------------------------------------------------------!
! Added 7/1/2019 to facilitate limiting percent coarse interpolation to wells
! and nodes of the same assigned geologic unit. File is assumed to have
! a heading and the node number followed by the geologic zones of each
! corresponding layer.
!
! Requires: strtracker (see Tools.f90) to increment zone arrays & counters
!
! Author: Leland Scantlebury of S.S. Papadopulos & Associates
!
! Updated 10/28/2020 by LS to use nodeids
!-----------------------------------------------------------------------------!

  integer                  :: nid, node, i, ilay, izone, ierr
  character(30)            :: geobylayer(nlayers)
  character(256)           :: cerr
  logical                  :: newzone, newzone_lay

  ! Initialize
  ngeozones = 0
  ngeozones_lay = 0

  write(*,'(2x,2a)') 'Reading Model Geologic Zone File ', trim(geozonefile)
  open(10, file =trim(geozonefile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
  read(10, *) ! Header

  do i=1, nnodes
    read(10, *) nid, geobylayer
    node = nodebyid(nid)
    if (node < 1) then
      write(*,*) 'ERROR! Invalid node ID in Geologic Zone File = ', nid
      stop
    end if
    geozones(node,:) = geobylayer
    ! Track geozones for looping, reporting and comparing
    do ilay=1, nlayers
      ! Is this zone new to the layer?
      call strtracker(ngeozones_lay(ilay), geobylayer(ilay), &
                      geonames_lay(:,ilay), newzone_lay)
      ! Handle if new zone for layer
      if (newzone_lay) then
        ! Is it a zone we haven't seen in any layer?
        call strtracker(ngeozones, geobylayer(ilay), &
                        geonames, newzone)
      end if
    end do
  end do

  close(10)

end subroutine readNodeGeoZones_IWFM
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine inputZones_IWFM()
    use MakePar
    use PilotPoints
    use errorhandle
    implicit none
!-----------------------------------------------------------------------------!
! Reads zone file consisting of node number and corresponding pilot point zone
!
! Updated 10/28/2020 by LS to use nodeids
! Modified 7/22/2022 by LS for different zone arrays per layer - but NOT fully
!   implemented! MODFLOW capability implemented first.
!-----------------------------------------------------------------------------!
    integer                :: i, nid, node, zone, ierr, iarr, &
                              zones_temp(max_zones,2), &
                              nzones_temp(2), &             ! eventually, allocatable
                              arr_temp(nnodes,1)            !arr_temp(nnodes,nlayers)
    logical                :: newzone
    character(256)         :: cerr

    iarr = 1 ! Placeholder until zone arrays are fully implemented
    pp_lay2arr(:) = 1  ! All share an array
    pp_lay2arr(nlayers+1:nlayers+natlayers) = 2  ! But aquitard arrays are separate

    open(10, file=trim(zonefile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    read(10,*)   ! Header
    do while (ierr == 0)
      read(10,*,iostat=ierr) nid, zone
      node = nodebyid(nid)
      if (node < 1) then
        write(*,*) 'ERROR! Invalid node ID in Pilot Point Zone File = ',nid
        stop
      end if
      arr_temp(node,iarr) = zone
      ! Is this a new zone?
      newzone = .true.
      do i=1, nzones_temp(iarr)
        if (zone == zones_temp(i,iarr)) newzone = .false.
      end do
      if (newzone) then
        nzones_temp(iarr) = nzones_temp(iarr) + 1
        ! Handle parameter exceedence if necessary
        if (nzones_temp(iarr) > max_zones) then
          call parexceed(max_zones, 'max_zones', &
                          'Failed during list creation of node zones', &
                          finalint=zone)
        end if
        zones_temp(nzones_temp(iarr),iarr) = zone
      end if
    end do

    ! Placeholder initialization until zones-by-layer is fully implemented
    ppaq_narr = 1
    ppat_narr = 1
    pp_narr = ppaq_narr + ppat_narr

    call PPZoneAllocate(nzones_temp(1),nnodes) ! nzones, nnodes

    ! Copy values over to permanent arrays
    pp_node_zone = arr_temp
    ! Note: for now aquifer zones == aquitard zones
    do i=1, ppaq_narr
      pp_arr_nzones(i) = nzones_temp(i)
      pp_arr_zones(:,i) = zones_temp(1:nzones_temp(i),i)
    end do
    nzones_temp(2) = nzones_temp(1)  ! Aquitard
    zones_temp(1:nzones_temp(2),2) = zones_temp(1:nzones_temp(1),1)
    do i=ppaq_narr+1, pp_narr
      pp_arr_nzones(i) = nzones_temp(i)
      pp_arr_zones(:,i) = zones_temp(1:nzones_temp(i),i)
    end do

    close(10)

end subroutine inputZones_IWFM
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine writeWellFile_IWFM(v, filename)
  use MakePar
  implicit none

  real,intent(in)          :: v(nwell,nlayers)
  character(*),intent(in)  :: filename
  integer                  :: iwell, ilay
  character(60)            :: fmt(2)

  write(*,'(4x,a,a35)') 'Writing:', filename

  write(fmt(1),'("(2a10,",i5,"i14)")') nlayers
  write(fmt(2),'("(2i10,",i5,"es14.5)")') nlayers

  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  ! Header
  write(20,fmt(1)) 'Well', 'Element', (ilay, ilay=1, nlayers)
  ! Values
  do iwell = 1, nwell
    if (wellelem(iwell) > 0) then  ! Make sure well exists
      write(20,fmt(2)) iwell, elemid(wellelem(iwell)), (v(iwell,ilay), ilay=1, nlayers)
    end if
  end do

  close(20)

end subroutine writeWellFile_IWFM
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine writeNodeFile(v, filename)
  use MakePar
  use PilotPoints
  implicit none

  real,intent(in)          :: v(nnodes,nlayers)
  character(*),intent(in)  :: filename
  integer                  :: inode, ilay
  character(60)            :: fmt(2)

  write(*,'(4x,a,a35)') 'Writing:', filename

  write(fmt(1),'("(2a8,2a18,",i5,"i14)")') nlayers
  write(fmt(2),'("(2i8,2f18.5,",i5,"e14.5)")') nlayers

  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  ! Header
  write(20,fmt(1)) 'Node', 'Zone', 'X', 'Y', (ilay, ilay=1, nlayers)
  ! Values
  do inode = 1, nnodes
    write(20,fmt(2)) nodeid(inode), pp_node_zone(inode,1), nodex(inode), nodey(inode), &
                     (v(inode, ilay), ilay=1,nlayers)
  end do

  close(20)

end subroutine writeNodeFile
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine writeIWFMgwfile()
  implicit none

  !integer, intent(in)      :: nskip
  integer                  :: i, ilay, ierr, nouth, noutf, ngroup
  character(256)           :: cerr
  character(1000)          :: line

  write(*,'(2x,a)') 'Reading gw temp file, writing new gw file'

  open(13, file=gwtemp, status='old', iostat = ierr, iomsg=cerr)   ! Temp File
  call iomsghandler(ierr, cerr)
  open(12, file=GWFile, status='replace', iostat = ierr, iomsg=cerr)   ! New File
  call iomsghandler(ierr, cerr)

  ! Copy Section Number
  call CopySettingsSection(13, 12, 1)
  ! Copy over the file header
  call CopyCommentSection(13, 12)
  ! Copy over main groundwater file settings (19 lines in V4.0)
  call CopySettingsSection(13, 12)
  ! Copy Debugging
  call CopyCommentSection(13, 12)
  call CopySettingsSection(13, 12)
  ! Copy Groundwater Hydrograph Output Settings
  call CopyCommentSection(13, 12)
  read(13, '(a1000)', iostat = ierr) line
  read(line, *) nouth
  write(12,'(a)') trim(line)
  call CopySettingsSection(13, 12)
  ! Copy Hydrogrph data
  call CopyCommentSection(13, 12)
  call CopySettingsSection(13, 12, nouth)
  ! Copy Element Face Flow Output Settings
  call CopyCommentSection(13, 12)
  read(13, '(a1000)', iostat = ierr) line
  read(line, *) noutf
  write(12,'(a)') trim(line)
  call CopySettingsSection(13, 12)
  ! Copy Element Face Data
  call CopyCommentSection(13, 12)
  call CopySettingsSection(13, 12, noutf)
  ! Copy Aquifer Parameters, Handle NGROUP setting
  call CopyCommentSection(13, 12)
  read(13, '(a1000)', iostat = ierr) line
  write(12,'(a)') trim(line)
  read(line, *) ngroup
  if (ngroup > 0) then
    write(*,'(a)') 'ERROR - NGROUP > 0 in Groundwater file (parametric grid)'
    write(*,'(a)') 'Only NGROUP = 0 is supported (node parameters)'
    stop
  end if
  ! Copy all the way down to aquifer parameter section
  call CopyCommentSection(13, 12)
  call CopySettingsSection(13, 12) ! Conversion Factors
  call CopyCommentSection(13, 12)
  call CopySettingsSection(13, 12) ! Units
  call CopyCommentSection(13, 12)

  ! Finally ready to write the aquifer parameters!
  do i=1, nnodes
    do ilay=1, nlayers
      if (ilay.eq.1) then
        write(12,'(i13,f14.8,e12.3,f14.8,f14.8,f14.8)') nodeid(i), KhB(i, ilay), SsB(i, ilay), SyB(i, ilay), KvB_aqtard(i, ilay), KvB(i, ilay)
      else
        write(12,'(13x,f14.8,e12.3,f14.8,f14.8,f14.8)') KhB(i, ilay), SsB(i, ilay), SyB(i, ilay), KvB_aqtard(i, ilay), KvB(i, ilay)
      end if
      read(13,*)
    end do
  end do

  ! Copy the rest of the file
  do
    read(13, '(a1000)', iostat = ierr) line
    if (ierr /= 0) exit
    write(12,'(a)') adjustl(trim(line))
  end do

  close(12)
  close(13)

end subroutine writeIWFMgwfile

!-----------------------------------------------------------------------------!
! IWFM File copying subroutines
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine CopyCommentSection(copyunit, newunit)
  implicit none

  integer,intent(in)       :: copyunit, newunit
  character(1)             :: compare
  character(1000)          :: line

  ! Copy lines until there are no more comment lines
  do
    read(copyunit, '(a1000)') line
    compare = ADJUSTL(line)
    if ((compare == 'C').or.(compare == 'c').or.&
        (compare == '*').or.(compare == '#')) then
      write(newunit,'(a)') trim(line)
    else
      ! Non comment found, back one line and exit
      backspace(copyunit)
      exit
    end if
  end do

end subroutine
!-----------------------------------------------------------------------------!
subroutine CopySettingsSection(copyunit, newunit, nlines)
  implicit none

  integer,intent(in)       :: copyunit, newunit
  integer                  :: i
  character(1)             :: compare
  character(1000)          :: line
  integer,intent(in),optional  :: nlines

  ! If nlines is passed, just copy over that many lines
  if (present(nlines)) then
    do i=1, nlines
      read(copyunit, '(a1000)') line
      write(newunit,'(a)') trim(line)
    end do
  else
    ! Copy lines until a comment line is found
    do
      read(copyunit, '(a1000)') line
      compare = adjustl(line)
      if ((compare == 'C').or.(compare == 'c').or.&
          (compare == '*').or.(compare == '#')) then
        ! Non comment found, back one line and exit
        backspace(copyunit)
        exit
      else
        write(newunit,'(a)') trim(line)
      end if
    end do
  end if

end subroutine

!-----------------------------------------------------------------------------!

subroutine add_elem2nodeelem(intd)
  use MakePar
  use errorhandle
  implicit none

  integer                  :: i, j, k, node
  integer, intent(in)      :: intd
  logical                  :: stored

  do j=1, 4
    stored = .false.
    node = elements(intd, j)
    if (node > 0) then
      ! Increment node-element membership array
      nnodeelem(node) = nnodeelem(node) + 1
      ! Ensure we're not going over our max
      if (nnodeelem(node) > max_node_elem_members) then
        write(*,*) 'ERROR - max_node_elem_members exceeded!'
        write(*,'(a,i7,a,i3)') 'node ', intd, 'is in >', max_node_elem_members
        stop
      end if
      ! Otherwise, add to the array
      nodeelem(node, nnodeelem(node)) = intd
    end if
  end do

end subroutine add_elem2nodeelem

!-----------------------------------------------------------------------------!

subroutine setup_nodeindex_byid(maxnodeid)
!-----------------------------------------------------------------------------!
! Sets up the nodebyid array so that it is easy to find the node index
! of any node by the nodeid. This comes in handy when reading subsequent files
! that refer to the node by the nodeid. Prevents having to loop through the
! whole array. If Fortran had 'maps' (C++) or 'dictionaries' (python) that
! would be a more ideal (efficient) data structure
! This structure could potentially have large gaps wasting memory.
! - LS 10/28/2020
!-----------------------------------------------------------------------------!
  use MakePar
  use errorhandle
  implicit none

  integer, intent(in)        :: maxnodeid
  integer                    :: i

  ! We're allocating an array based on the maximum number
  ! So a quick check to make sure it's not too big of a number!
  if (maxnodeid > max_nodeid) then
    call parexceed(maxnodeid, 'max_nodeid', &
                          'Failed during creation of nodebyid')
  end if

  allocate(nodebyid(maxnodeid))

  ! Initialize all values to a negative number so we can easily check whether
  ! a nodeid is valid
  nodebyid = -1

  ! Fill array
  do i=1, nnodes
    nodebyid(nodeid(i)) = i
  end do

end subroutine setup_nodeindex_byid

!-----------------------------------------------------------------------------!

subroutine intread(IFAIL,CLINE,iTEMP)
      integer, intent(out)            ::IFAIL
      character (len=*), intent(in)   ::cline
      integer, intent(out)            ::iTEMP

! -- Subroutine intREAD reads an integer number from a string.


! -- Subroutine arguments are as follows:-
!       ifail:    returned as non-zero in case of failure
!       cline:    character string
!       itemp:    return integer

       CHARACTER*6 AFMT

       IFAIL=0
       AFMT='(i   )'
       WRITE(AFMT(3:5),'(I3)') LEN(CLINE)
       READ(CLINE,AFMT,ERR=100) iTEMP

       RETURN

100    IFAIL=1
       RETURN
end subroutine intread

!-----------------------------------------------------------------------------!

subroutine readOneLine(infile,iinfile,cline)
!   read and return one line from a file, skipping lines that are commented out

    implicit none

    character (len=*), intent(in)  :: infile
    integer,           intent(in)  :: iinfile
    character (len=*), intent(out) :: cline

100  read(iinfile,'(a)',end=900) cline
        if((cline(1:1).eq.'c').or.(cline(1:1).eq.'C').or.&
            (cline(1:1).eq.'*').or.(cline(1:1).eq.'#')) go to 100    ! skip comment line

    go to 999 ! skip past error statements to end of subroutine

    !================================================================================
    ! error statements
900 write(*,901) trim(infile)
901 format(/,' *** Unexpected end to file ',a,' ***',/)
      stop
999 end subroutine readOneLine

!-----------------------------------------------------------------------------!

    SUBROUTINE multisplit(IFAIL,NUM,LW,RW,CLINE)

! -- Subroutine multisplit splits a string into blank-delimited fragments.

! -- Subroutine arguments are as follows:-
!       ifail:    returned as non-zero in case of failure
!       num:      number of ...
!       lw:       number of ...
!       rw:       number of ...
!       cline:    character string

! -- Author:-
!       John Doherty

       INTEGER IFAIL,NW,NBLC,J,I
       INTEGER NUM,NBLNK
       INTEGER LW(NUM),RW(NUM)
       CHARACTER*(*) CLINE
       IFAIL=0
       NW=0
       NBLC=LEN_TRIM(CLINE)
       IF((NBLC.NE.0).AND.(INDEX(CLINE,CHAR(9)).NE.0)) THEN
         CALL TABREM(CLINE)
         NBLC=LEN_TRIM(CLINE)
       ENDIF
       IF(NBLC.EQ.0) THEN
         IFAIL=-1
         RETURN
       END IF
       J=0
5      IF(NW.EQ.NUM) RETURN
       DO 10 I=J+1,NBLC
         IF((CLINE(I:I).NE.' ').AND.(CLINE(I:I).NE.',').AND.&
         (ICHAR(CLINE(I:I)).NE.9)) GO TO 20
10     CONTINUE
       IFAIL=1
       RETURN
20     NW=NW+1
       LW(NW)=I
       DO 30 I=LW(NW)+1,NBLC
         IF((CLINE(I:I).EQ.' ').OR.(CLINE(I:I).EQ.',').OR.&
         (ICHAR(CLINE(I:I)).EQ.9)) GO TO 40
30     CONTINUE
       RW(NW)=NBLC
       IF(NW.LT.NUM) IFAIL=1
       RETURN
40     RW(NW)=I-1
       J=RW(NW)
       GO TO 5

    END subroutine multisplit

!-----------------------------------------------------------------------------!

  subroutine TABREM(CLINE)

! -- Subroutine TABREM removes tabs from a string.

! -- Subroutine arguments are as follows:-
!       cline:    character string


       INTEGER I
       CHARACTER*(*) CLINE

       DO 10 I=1,LEN(CLINE)
10     IF(ICHAR(CLINE(I:I)).EQ.9) CLINE(I:I)=' '

       RETURN
  end subroutine tabrem

!-----------------------------------------------------------------------------!

end module IWFM