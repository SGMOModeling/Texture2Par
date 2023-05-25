module MFSupport
  implicit none

  SAVE
  real                      :: xoff, yoff, mfrot
  real,allocatable          :: mftop(:,:,:), mfbot(:,:,:), locx(:), locy(:)
  integer                   :: mfmaxunit, layfile, ncbd
  character(3)              :: gwftyp
  character(200)            :: layfiletemp
  logical                   :: lop

  CONTAINS

!-----------------------------------------------------------------------------!
  subroutine read_modflow(namfile)

! MODFLOW globals to import
    use GLOBAL, ONLY : IUNSTR,NLAY,NROW,NCOL,DELR,DELC,TOP,BOT,LAYCBD,IUNIT
    use PARAMMODULE, ONLY: NMLTAR, NZONAR
    use GWFBCFMODULE, ONLY: VKA, SC1, SC2

! Texture2Par modules
    use MakePar
    use errorhandle

    implicit none

! Local variables
    character(200)         :: namfile
    integer                :: ierr, i, j, k

! Open modflow files
    call init_modflow_usgs(namfile,ierr,mfmaxunit)
    call iostathandler(ierr)
    if (IUNSTR==0) then
!      write(*,*) 'Structured Model'
    else if (IUNSTR==1) then
      call USGunimplemented()
    end if

! Set MakePar module values
    nnodes = nrow * ncol
    nlayers = nlay
    natlayers = 0

! Handle MODFLOW settings
    ! UPW or LPF (or none)
    layfile = 0
    if (IUNIT(23) > 0) then
      layfile = 23
      gwftyp = 'LPF'
    end if
    if (IUNIT(45) > 0) then
      layfile = 45
      gwftyp = 'UPW'
    end if
    if (layfile==0) call NotInNam('LPF or UPW')

! Allocate
    call allocate_modflow()

! Get Top & Bot arrays into node format
    ncbd = 0
    Elevation(1:nnodes) = top(1:nnodes)
    do k=1, nlay
      topelev(1:nnodes,k) = top((k-1)*nnodes+1 : k*nnodes)
      botelev(1:nnodes,k) = bot((k-1)*nnodes+1 : k*nnodes)
      ! Quasi-3d confining beds == Aquitards
      if (laycbd(k) > 0) then
        ncbd = ncbd + 1
        natlayers = natlayers + 1
        aqtardtopelev(1:nnodes,k) = bot((k-1)*nnodes+1 : k*nnodes)
        aqtardbotelev(1:nnodes,k) = top(k*nnodes+1 : (k+1)*nnodes)
      else
        ! Fill with junk data
        aqtardtopelev(1:nnodes,k) = NODATA
        aqtardbotelev(1:nnodes,k) = NODATA
      end if
    end do

  end subroutine read_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine calcCellCenters_global()
    use GLOBAL, ONLY : NROW, NCOL, DELR, DELC
    use MakePar
    ! Calculates cell centers in global coordinates

    integer                :: i, j, n
    real                   :: xbase, ybase

    ybase = sum(delc)
    do i=1, nrow
      xbase = 0.0
      do j=1, ncol
        n = (i-1)*ncol+j
        locx(n) = xbase + 0.5 * delr(j)
        locy(n) = ybase - 0.5 * delc(i)
        call loc2glo(locx(n), locy(n), xoff, yoff, mfrot, nodex(n), nodey(n))
        xbase = xbase + delr(j)
      end do
      ybase = ybase - delc(i)
    end do

  end subroutine calcCellCenters_global
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
SUBROUTINE inputwells_modflow(readzones)
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
! Well GeoZones added 7/3/2019 -LS
!   Argument "readzones" is a flag for whether GeoZones exist (1 = yes, 0 = no)
!-----------------------------------------------------------------------------!

    integer, parameter       :: nnear=4
    integer                  :: i, iwell, izone, igeo, ilay, readzones, &
                                currwellelem, node, ierr, nwellgeozones
    integer                  :: closenodes(nnear)
    real                     :: depth, tx, ty, tz, td, tpc
    real                     :: nodedists(nnear), queryxy(2)
    character(30)            :: wellname, geobylayer(nlayers), &
                                wellgeonames(max_zones)
    character(265)           :: cerr
    logical                  :: newzone

    ! Initialize (to prevent issues where compiler does not initilize to 0)
    nwellgeozones = 0

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
        ! Well "zone" is discrete depth interval and should not be confused with
        ! other uses of zone (e.g. pilot points)
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
          node = get_well_cell(tx, ty, closenodes, nnear)
          wellelem(iwell) = node
          ! While looping over wells, might as well assign modflow
          ! elevations to wells, as no interpolation is required
          wellelevtop(iwell,1:nlayers) = topelev(node,1:nlayers)
          wellelevbot(iwell,1:nlayers) = botelev(node,1:nlayers)
          if (ncbd > 0) then
            wellaqtardelevtop(iwell,1:nlayers) = aqtardtopelev(node,1:nlayers)
            wellaqtardelevbot(iwell,1:nlayers) = aqtardbotelev(node,1:nlayers)
          end if

        end if
        Zwell(iwell, izone) = tz - td
      end do

    else
      ! GeoZones
      do i=1, wfilelen
        read(10,*) wellname, iwell, izone, tpc, tx, ty, tz, td, geobylayer(:)
        ! Well "zone" is discrete depth interval and should not be confused with
        ! other uses of zone (e.g. pilot points, geologic)
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
          node = get_well_cell(tx, ty, closenodes, nnear)
          wellelem(iwell) = node
          ! While looping over wells, might as well assign modflow
          ! elevations to wells, as no interpolation is required
          wellelevtop(iwell,1:nlayers) = topelev(node,1:nlayers)
          wellelevbot(iwell,1:nlayers) = botelev(node,1:nlayers)
          if (ncbd > 0) then
            wellaqtardelevtop(iwell,1:nlayers) = aqtardtopelev(node,1:nlayers)
            wellaqtardelevbot(iwell,1:nlayers) = aqtardbotelev(node,1:nlayers)
          end if

        end if
        Zwell(iwell, izone) = tz - td
      end do
    end if
    close(10)

    !write(*,'(4x,a)') 'Well file read successfully.'

end subroutine inputwells_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine readNodeGeoZones_modflow()
  use GLOBAL, ONLY : IUNSTR, NLAY, NROW, NCOL
  use MakePar
  use str_utls
  use errorhandle
  implicit none
!-----------------------------------------------------------------------------!
! Added 7/1/2019 to facilitate limiting percent coarse interpolation to wells
! and cells of the same assigned geologic unit. File is assumed to have
! a heading and row, col followed by the geologic zones of each
! corresponding layer.
!
!
! Requires: strtracker (see Tools.f90) to increment zone arrays & counters
!
! Author: Leland Scantlebury of S.S. Papadopulos & Associates
!-----------------------------------------------------------------------------!

  integer                  :: node, i, ilay, izone, ierr, row, col
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
    read(10, *) row, col, geobylayer
    geozones((row-1)*ncol + col,:) = geobylayer
    ! Track geozones for looping, reporting and comparing
    do ilay=1, nlayers
      ! Is this zone new to the layer?
      call strtracker(ngeozones_lay(ilay), geobylayer(ilay), &
                      geonames_lay(:, ilay), newzone_lay)
      ! Handle if new zone for layer
      if (newzone_lay) then
        ! Is it a zone we haven't seen in any layer?
        call strtracker(ngeozones, geobylayer(ilay), &
                        geonames, newzone)
      end if
    end do
  end do

  close(10)

end subroutine readNodeGeoZones_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine inputZones_modflow()
    use GLOBAL, ONLY : NLAY
    use PilotPoints
    use errorhandle
    use MakePar, ONLY : nlayers,natlayers,nnodes,max_zones
!-----------------------------------------------------------------------------!
! Reads input file with pilot point zones for nodes
!
! Modified 7/21/2022 by LS:
!  - Custom zone file instead of using MODFLOW zones
!  - PP Zones can be specified by layer
! Author: Leland Scantlebury of S.S. Papadopulos & Associates
!-----------------------------------------------------------------------------!
    integer                :: ierr, i, k, temp, nzones, n
    integer                :: arr_temp(nnodes,nlayers)
    integer,allocatable    :: temp_ints(:)
    character(200)         :: line
    character(256)         :: cerr

    ! Zones are input in MODFLOW-style arrays
    nzones = 0
    pp_narr = 0
    ppaq_narr = 0
    ppat_narr = 0
    pp_lay2arr = 1
    
    open(10, file =trim(zonefile), status='old', action='read', iostat=ierr, iomsg=cerr)
    call iomsghandler(ierr, cerr)
    
    ! Arrays are by Layer
    do k=1, nlayers
      ! Skip any comment lines
      do
        read(10, *) line
        if (line(1:1) == '#') cycle
        backspace(10)
        exit
      end do
      ! Check if first number is negative -> reuse previous layer zones
      read(10, *) temp
      if (temp < 0) then
        if (k==1) then
          write(*,*) 'ERROR - First layer pilot point zone set to reuse non-existant previous array'
          STOP
        else
          ! Handle reuse array
          pp_lay2arr(k) = pp_lay2arr(k-1)
        end if
      else
        pp_narr = pp_narr + 1
        ppaq_narr = ppaq_narr + 1
        pp_lay2arr(k) = pp_narr
        backspace(10)
        read(10, *) arr_temp(:,pp_narr)
      end if
    end do

    ! Aquitard Arrays are also by Layer (well, lack of layer)
    do k=1, natlayers
      ! Skip any comment lines
      do
        read(10, *) line
        if (line(1:1) == '#') cycle
        backspace(10)
        exit
      end do
      ! Check if first number is negative -> reuse previous layer zones
      read(10, *) temp
      if (temp < 0) then
        if (k==1) then
          write(*,*) 'ERROR - First layer pilot point zone set to reuse non-existant previous array'
          STOP
        else
          ! Handle reuse array
          pp_lay2arr(k) = pp_lay2arr(k-1)
        end if
      else
        pp_narr = pp_narr + 1
        ppat_narr = pp_narr + 1
        pp_lay2arr(k) = pp_lay2arr(pp_narr)
        backspace(10)
        read(10, *) arr_temp(:,pp_narr)
      end if
    end do
    
    ! Now we can allocate the zone arrays
    nzones = maxval(arr_temp(:,1:pp_narr))  ! Perhaps too much work for a temporary value
    if (nzones > max_zones) then
      write(*,*) 'Too many unique zones specified in ',trim(zonefile)
      STOP
    end if    
    call PPZoneAllocate(nzones, nnodes)

    ! Copy over the zones:
    pp_node_zone = arr_temp

    ! Now we do some counting
    do i=1, ppaq_narr
      ! Get unique values
      call get_unique_ints(nnodes, arr_temp(:,i), nzones, pp_arr_zones(:,i), pp_arr_nzones(i))
    end do
    ! Aquitards
    do i=ppaq_narr+1, pp_narr
      ! Get unique values
      call get_unique_ints(nnodes, arr_temp(:,i), nzones, pp_arr_zones(:,i), pp_arr_nzones(i))
    end do

  end subroutine inputZones_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
integer function get_well_cell(x, y, closenodes, nnear) result(cell)
  implicit none

  real, intent(in)         :: x, y
  integer, intent(in)      :: nnear, closenodes(nnear)
  integer                  :: i
  real                     :: xl, yl, xmin, xmax, ymin, ymax

  cell = 0
  ! Check each cell center to see if well belongs
  do i=1, nnear
    call getCellBounds_local(closenodes(i), xmin, xmax, ymin, ymax)
    call glo2loc(x, y, xoff, yoff, mfrot, xl, yl)
    if ((xl >= xmin).and.(xl < xmax).and.(yl >= ymin).and.(yl < ymax)) then
      cell = closenodes(i)
      exit
    end if
  end do

  if (cell == 0) then
    ! Maybe this should just be a warning?
    write(*,*) 'ERROR well not found in a MODFLOW cell'
    write(*,*) x, y
    stop
  end if

end function get_well_cell
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine getCellBounds_local(node, xmin, xmax, ymin, ymax)
    use MakePar
    use GLOBAL, ONLY : DELR, DELC

    integer,intent(in)     :: node
    real,intent(out)       :: xmin,xmax,ymin,ymax
    integer                :: row, col

    call node2rowcol(node, row, col)

    xmin = locx(node) - 0.5 * delr(col)
    xmax = locx(node) + 0.5 * delr(col)
    ymin = locy(node) - 0.5 * delc(row)
    ymax = locy(node) + 0.5 * delc(row)

  end subroutine getCellBounds_local
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine getCellBounds_global(node, xmin, xmax, ymin, ymax)
    use MakePar
    use GLOBAL, ONLY : DELR, DELC

    ! Not compatible with rotated grids - not currently used in Texture2Par

    integer,intent(in)     :: node
    real,intent(out)       :: xmin,xmax,ymin,ymax
    integer                :: row, col

    call node2rowcol(node, row, col)

    xmin = nodex(node) - 0.5 * delr(col)
    xmax = nodex(node) + 0.5 * delr(col)
    ymin = nodey(node) - 0.5 * delc(row)
    ymax = nodey(node) + 0.5 * delc(row)

  end subroutine getCellBounds_global
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine node2rowcol(node, row, col)
    use GLOBAL, ONLY : NROW, NCOL

    integer,intent(in)     :: node
    integer,intent(out)    :: row, col

    row = CEILING(real(node)/real(ncol))
    col = node - (row-1)*ncol

  end subroutine node2rowcol
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine allocate_modflow()
    use GLOBAL, ONLY : IUNSTR, NLAY, NROW, NCOL, DELR, DELC, TOP, BOT, NODES
    use MakePar
    implicit none

    allocate(locx(nnodes), &
             locy(nnodes)  )

    ! TODO move to MakePar routine for allocation
    ! Since these are all use for both models
    allocate(nodex(nnodes),              &
             nodey(nnodes),              &
             nodexy(nnodes,2),           &
             Elevation(nnodes),          &
             topelev(nnodes,nlay),       &
             botelev(nnodes,nlay),       &
             aqtardtopelev(nnodes,nlay), &
             aqtardbotelev(nnodes,nlay)  )

  end subroutine allocate_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine writeWellFile_modflow(v, filename)
  use MakePar
  implicit none

  real,intent(in)          :: v(nwell,nlayers)
  character(*),intent(in)  :: filename
  integer                  :: iwell, ilay, row, col
  character(60)            :: fmt(2)

  write(*,'(4x,a,a35)') 'Writing:', filename

  write(fmt(1),'("(3a10,",i5,"i14)")') nlayers
  write(fmt(2),'("(3i10,",i5,"es14.5)")') nlayers

  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  ! Header
  write(20,fmt(1)) 'Well', 'Row', 'Col', (ilay, ilay=1, nlayers)
  ! Values
  do iwell = 1, nwell
    call node2rowcol(wellelem(iwell), row, col)
    write(20,fmt(2)) iwell, row, col, (v(iwell,ilay), ilay=1, nlayers)
  end do

  close(20)

end subroutine writeWellFile_modflow
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
subroutine writeCellFile(v, filename)
  use PilotPoints
  use MakePar
  implicit none

  real,intent(in)          :: v(nnodes, nlayers)
  character(*),intent(in)  :: filename
  integer                  :: inode, ilay, row, col
  character(60)            :: fmt(2)
  character(60)            :: zonetxt

  write(*,'(4x,a,a35)') 'Writing:', filename

  write(zonetxt,"(i5,a)") nlayers, "(' Zone',i3)"
  write(fmt(1),'("(2a7,",a,",2a18,",i5,"i14)")') trim(zonetxt), nlayers
  write(fmt(2),'("(2i7,",i5,"i8,2f18.5,",i5,"e14.5)")') nlayers, nlayers

  open(20, file=trim(filename), status='replace', action='write', buffered='YES')
  ! Header
  write(20,fmt(1)) 'Row', 'Col', (ilay, ilay=1, nlayers), 'X', 'Y', (ilay, ilay=1, nlayers)
  ! Values
  do inode = 1, nnodes
    call node2rowcol(inode, row, col)
    write(20,fmt(2)) row, col, (pp_node_zone(inode,pp_lay2arr(ilay)), ilay=1, nlayers), &
                     nodex(inode), nodey(inode), v(inode,:)
  end do

  close(20)

end subroutine writeCellFile
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
  subroutine readwriteFlowPackage()
    use GLOBAL, ONLY:IUNIT
    use errorhandle
    implicit none

    integer                :: ierr, upwflag
    character(300)         :: fname, cerr
    logical                :: lop

    write(*,'(2x,2(a,3a))') 'Reading temp ', gwftyp, ', writing new ', gwftyp

    ! If LPF/UPW is open, close it. We need to do a bait-and-switch on MODFLOW
    inquire(unit=IUNIT(layfile), opened=lop, name=fname)
    if (lop) close(IUNIT(layfile))

    ! Initialize MODFLOW LPF values with template LPF file (user provided)
    open(IUNIT(layfile),file=trim(layfiletemp),iostat=ierr,status='old',&
         iomsg=cerr)
    call iomsghandler(ierr, cerr)
    upwflag = 0
    if (gwftyp=='UPW') upwflag = 1
    CALL GWF2BCFU1AR(IUNIT(1),IUNIT(22),IUNIT(layfile),IUNIT(64), upwflag)
    close(IUNIT(layfile))

    ! Re-open actual LPF file and OVERWRITE
    open(IUNIT(layfile),file=trim(fname),iostat=ierr,status='replace',&
         iomsg=cerr)
    call iomsghandler(ierr, cerr)
    call writeGWFP(IUNIT(layfile))
    close(IUNIT(layfile))

  end subroutine readwriteFlowPackage
!-----------------------------------------------------------------------------!

end module MFSupport