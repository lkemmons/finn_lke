
   program finn_fire_tst
!-----------------------------------------------------------------------------
! GLOBAL FIRE EMISSIONS ESTIMATES FOR MOZART MODEL
! This program will calculate global daily fire emisisons
! This is the prototype for the MIRAGE and INTEX campaigns
! Developed by Christine Wiedinmyer, December 30, 2005
!
!Inputs to this model include 1 file with overlaid data:
! 1) Fire detections/locations/times from Rapid Response web site
! 2) MODIS Land Cover Type (MOD12C, using the IGBP land cover classification
!    available online at http://edcdaac.usgs.gov/modis/mod12c1v4.asp)
! 3) VCF % tree and % herbaceous data
! 4) Model domains for the MOZART and WRF-Chem simulations
!
! DECEMBER 20, 2011
! - Edited the comments in this code for Stacy
!-----------------------------------------------------------------------------

   implicit none

!-----------------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------------
   integer, parameter :: i_unit = 10           ! raw fire file
   integer, parameter :: o_unit = 11           ! intermediate fire output file
   integer, parameter :: l_unit = 12           ! intermediate fire log file
   integer, parameter :: c_unit = 13           ! speciated fire log file
   integer, parameter :: s_unit = 14           ! speciated fire output file
   integer :: yearnum
   integer :: astat, ios
   integer :: nlines, tokcnt, ntokens, slen, ntoday, novrlap
   integer :: nfires, ngoodfires, nfire1, ntropics, nconfgt20
   integer :: wrk_int, reg
   integer :: yr, mnth, day, ntotdays
   integer :: i, il, iu, j, jl, ju, m, n, n1, ndx
   integer :: genveg
   integer :: lct0 = 0, spixct = 0, antarc = 0, allbare = 0
   integer :: genveg0 = 0, bmass0 = 0, vcfcount = 0
   integer :: vcflt50 = 0, overlapct = 0, urbnum = 0
   integer :: spec_ndx = 0
   integer :: speciation_start_day = -999
   integer :: speciation_end_day   = 999
   integer, pointer :: lct(:)
   integer, pointer :: globreg(:)
   integer, pointer :: jd(:)
   integer, pointer :: tod(:)
   integer, allocatable :: toklen(:)
   integer, allocatable :: iwrk(:)
   integer, allocatable :: flag(:)
   integer, allocatable :: today_ndx(:),ovrlap_ndx(:)
   integer :: max_loc(1)
   integer(1), pointer  :: Modis_LCT(:,:), glb_region(:,:)

   real, parameter :: rearthkm = 6371.    ! km
   real, parameter :: scale    = 2.
   real, parameter :: CO2_mw  = 44.01
   real, parameter :: CO_mw   = 28.01
   real, parameter :: CH4_mw  = 16.04
   real, parameter :: NH3_mw  = 17.03
   real, parameter :: NO_mw   = 30.01
   real, parameter :: NO2_mw  = 46.01
   real, parameter :: SO2_mw  = 64.06
   real, parameter :: H2_mw   = 2.02

   real    :: FAC = 0.5                   ! Factor to apply to duplicate tropical fires 
   real    :: wrk_real, wrk_reali, wrk_lat
   real    :: pi, dtor, area, bmass, bmass1, bmassburn
   real    :: xtrack, atrack
   real    :: dxdump, dydump, dxdumm, dydumm
   real    :: CF1, CF3, CF4, CF5, CF6 
   real    :: pctherb, pcttree, herbbm, coarsebm
!-----------------------------------------------------------------------------
! For the total biomass burned in each genveg for output file 
! Added 08/24/2010
!-----------------------------------------------------------------------------
   real :: TOTTROP = 0.0
   real :: TOTTEMP = 0.0
   real :: TOTBOR  = 0.0
   real :: TOTSHRUB = 0.0
   real :: TOTCROP  = 0.0
   real :: TOTGRAS  = 0.0
 
!-----------------------------------------------------------------------------
! For the total area in each genveg for output log file
! added 06/21/2011
!-----------------------------------------------------------------------------
   real :: TOTTROParea = 0.0
   real :: TOTTEMParea = 0.0
   real :: TOTBORarea  = 0.0
   real :: TOTSHRUBarea = 0.0
   real :: TOTCROParea  = 0.0
   real :: TOTGRASarea  = 0.0

   real :: CO2
   real :: CO
   real :: CH4
   real :: NMHC
   real, pointer :: VOC
   real, target  :: NMOC
   real :: H2
   real :: NOX
   real :: NO
   real :: NO2
   real :: SO2
   real :: PM25
   real :: TPM
   real :: TPC
   real :: OC
   real :: BC
   real :: NH3
   real :: PM10

   real :: CO2total  = 0.0
   real :: COtotal   = 0.0
   real :: CH4total  = 0.0
   real :: NMHCtotal = 0.0
   real :: NMOCtotal = 0.0
   real :: H2total   = 0.0
   real :: NOXtotal  = 0.0
   real :: NOtotal   = 0.0
   real :: NO2total  = 0.0
   real :: SO2total  = 0.0
   real :: PM25total = 0.0
   real :: TPMtotal  = 0.0
   real :: TPCtotal  = 0.0
   real :: OCtotal   = 0.0
   real :: BCtotal   = 0.0
   real :: NH3total  = 0.0
   real :: PM10total = 0.0
   real :: AREAtotal = 0.0 ! added 06/21/2011
   real :: BMASStotal= 0.0 ! Addded 06/21/2011

   real, allocatable :: rwrk(:)
   real, pointer :: lon(:)
   real, pointer :: lat(:)
   real, pointer :: spix(:)
   real, pointer :: tpix(:)
   real, pointer :: CONF(:)
   real, pointer :: tree(:)
   real, pointer :: herb(:)
   real, pointer :: bare(:)
   real, pointer :: factortrop(:)
   real, pointer :: totcov(:)
   real, allocatable :: xearth(:)
   real, allocatable :: yearth(:)
   real, allocatable :: speciated_emissions(:)
   real, allocatable :: speciation(:,:)

   real(8), pointer :: herb_lons(:), herb_lats(:)
   real(8), pointer :: tree_lons(:), tree_lats(:)

   character(len=1), parameter  :: comma = ','
   character(len=128) :: buffer
   character(len=128) :: fuel_load_filespec
   character(len=128) :: emis_factor_filespec
   character(len=128) :: intermediate_Outfile_filespec = ' '
   character(len=128) :: intermediate_logfile_filespec = ' '
   character(len=128) :: speciated_logfile_filespec = ' '
   character(len=128) :: speciated_Outfile_filespec = ' '
   character(len=128) :: raw_infile_filespec
   character(len=128) :: speciation_infile_filespec
   character(len=128) :: bare_frac_filespec = ' '
   character(len=128) :: herb_frac_filespec = ' '
   character(len=128) :: tree_frac_filespec = ' '
   character(len=128) :: lct_filespec = ' '
   character(len=128) :: globreg_filespec = ' '
   character(len=128) :: i_frmt, s_frmt
   character(len=32)  :: var_name
   character(len=16)  :: speciation_case = ' '
   character(len=10)  :: speciation_start_date = ' '
   character(len=10)  :: speciation_end_date   = ' '
   character(len=32), allocatable :: tokens(:)
   character(len=10), allocatable :: date(:)
   character(len=10), allocatable :: cwrk(:)

   logical :: in_tropics
   logical :: has_intermediate_output
   logical :: has_intermediate_log
   logical :: has_speciated_output
   logical :: has_speciated_log
   logical :: has_speciation = .false.
   logical, allocatable :: mask(:)

!-----------------------------------------------------------------------------
!  ... defined types
!-----------------------------------------------------------------------------
   type fuel_load_type
     integer :: globreg2
     real :: tffuel                  ! tropical forest fuels
     real :: tefuel                  ! temperate forest fuels
     real :: bffuel                  ! boreal forest fuels
     real :: wsfuel                  ! woody savanna fuels
     real :: grfuel                  ! grassland and savanna fuels
   end type

   type emis_factor_type
     integer :: lctemis              ! LCT Type (Added 10/20/2009)
     integer :: vegemis              ! generic vegetation type --> this is ignored in model
     real :: CO2EF                   ! CO2 emission factor
     real :: COEF                    ! CO emission factor
     real :: CH4EF                   ! CH4 emission factor
     real :: NMHCEF                  ! NMHC emission factor
     real :: NMOCEF                  ! NMOC emission factor (added 10/20/2009)
     real :: H2EF                    ! H2 emission factor
     real :: NOXEF                   ! NOx emission factor
     real :: NOEF                    ! NO emission factors (added 10/20/2009)
     real :: NO2EF                   ! NO2 emission factors (added 10/20/2009)
     real :: SO2EF                   ! SO2 emission factor
     real :: PM25EF                  ! PM2.5 emission factor
     real :: TPMEF                   ! TPM emission factor
     real :: TCEF                    ! TC emission factor
     real :: OCEF                    ! OC emission factor
     real :: BCEF                    ! BC emission factor
     real :: NH3EF                   ! NH3 emission factor
     real :: PM10EF                  ! PM10 emission factor (added 08/18/2010)
   end type emis_factor_type
   
   type(fuel_load_type), allocatable   :: fuel_load(:)
   type(emis_factor_type), allocatable :: emis_factor(:)

!-----------------------------------------------------------------------------
!  ... namelist
!-----------------------------------------------------------------------------
   namelist /control/ yearnum, speciation_case, fuel_load_filespec, emis_factor_filespec, FAC, &
                      intermediate_Outfile_filespec, intermediate_logfile_filespec, raw_infile_filespec, &
                      speciation_infile_filespec, speciated_Outfile_filespec, speciated_logfile_filespec, &
                      speciation_start_date, speciation_end_date, bare_frac_filespec, herb_frac_filespec, &
                      tree_frac_filespec, lct_filespec, globreg_filespec

!-------------------------------------------------------------------------
!  read control variables
!-------------------------------------------------------------------------
   read(*,nml=control,iostat=ios)
   if( ios /= 0 ) then
     write(*,*) 'finn_fire: failed to read namelist; error = ',ios
     stop
   end if

   has_intermediate_output = trim( intermediate_Outfile_filespec ) /= ' '
   has_intermediate_log    = trim( intermediate_logfile_filespec ) /= ' '
   has_speciated_output    = trim( speciated_Outfile_filespec ) /= ' '
   has_speciated_log       = trim( speciated_logfile_filespec ) /= ' '

   VOC => NMOC

   select case( trim(speciation_case) )
     case ( 'MOZ4' )
       spec_ndx = 1
       has_speciation = .true.
     case ( 'SAPRC99' )
       spec_ndx = 2
       has_speciation = .true.
     case ( 'GEOSCHEM' )
       spec_ndx = 3
       has_speciation = .true.
   end select
   has_speciation = has_speciation .and. (has_speciated_output .or. has_speciated_log)

!-----------------------------------------------------------------------------
! ASSIGN FUEL LOADS, EMISSION FACTORS FOR GENERIC LAND COVERS AND REGIONS
!   - created tables to be read in, instead of hardwiring the values here
! 08/25/08: Edited pathways here
! 10/19/2009: created new input files
!-----------------------------------------------------------------------------
!   read in the fuel loading file
!   NOTE: Input fuel loads have units of g/m2 DM
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!  ... read in the fuel load file
!-----------------------------------------------------------------------------
   call read_fuel_load_file
!-----------------------------------------------------------------------------
!  ... read in the emission factors file
!-----------------------------------------------------------------------------
   call read_emis_factors_file

   write(*,*) ' '
   write(*,*) 'fuel_load(1)'
   write(*,*) fuel_load(1)
   write(*,*) ' '
   write(*,*) 'fuel_load(nlines)'
   write(*,*) fuel_load(nlines)
   write(*,*) ' '
   write(*,*) 'emis_factor(1)'
   write(*,*) emis_factor(1)
   write(*,*) ' '
   write(*,*) 'emis_factor(nlines)'
   write(*,*) emis_factor(nlines)

   write(*,*) ' '
   write(*,*) 'Finished reading in fuel and emission factor files'

!-----------------------------------------------------------------------------
!  ... read in the speciation
!-----------------------------------------------------------------------------
   if( has_speciation ) then
     call read_speciation_file
   endif

!-----------------------------------------------------------------------------
!  ... setup the output and log files
!-----------------------------------------------------------------------------
   call setup_output_and_log_files

!-----------------------------------------------------------------------------
!  ... the raw fire input file
!-----------------------------------------------------------------------------
   call read_raw_fire_file
   write(*,*) ' '
   write(*,*) 'Finished reading raw input file ',trim(raw_infile_filespec)

!-----------------------------------------------------------------------------
!  ... the global region
!-----------------------------------------------------------------------------
   var_name = 'glob_reg'
   call get_lct_fractions( globreg_filespec, int_fraction=globreg, do_average=.false. )
   write(*,*) 'min,max globreg        = ',minval(globreg(:)),maxval(globreg(:))
!-----------------------------------------------------------------------------
!  ... the Modis LCT
!-----------------------------------------------------------------------------
   var_name = 'ModisLCT_rcls'
   call get_lct_fractions( lct_filespec, int_fraction=lct, do_average=.false. )
   write(*,*) 'min,max lct            = ',minval(lct(:)),maxval(lct(:))
!-----------------------------------------------------------------------------
!  ... the bare land fractional coverage
!-----------------------------------------------------------------------------
   var_name = 'bare_globe'
   call get_lct_fractions( bare_frac_filespec, real_fraction=bare, do_average=.true. )
!-----------------------------------------------------------------------------
!  ... the herb land fractional coverage
!-----------------------------------------------------------------------------
   var_name = 'herb_globe2'
   call get_lct_fractions( herb_frac_filespec, real_fraction=herb, do_average=.true. )
!-----------------------------------------------------------------------------
!  ... the tree land fractional coverage
!-----------------------------------------------------------------------------
   var_name = 'tree_globe'
   call get_lct_fractions( tree_frac_filespec, real_fraction=tree, do_average=.true. )

   write(*,*) '# fires with bare+herb+tree > 100 = ',count((bare(:)+herb(:)+tree(:)) > 100.)
   max_loc(:) = maxloc( bare(:)+herb(:)+tree(:) )
   write(*,*) 'min,max bare+herb+tree = ',minval((bare(:)+herb(:)+tree(:))), &
                                          maxval((bare(:)+herb(:)+tree(:)))
   write(*,*) 'min,max bare+herb+tree fire index = ',max_loc(:)
   write(*,*) 'min,max lct            = ',minval(lct(:)),maxval(lct(:))

   i_frmt = '(D20.10,",",D20.10,",",(8(I10,",")),2(D20.10,","),17(D25.10,","),F6.3)'
!-----------------------------------------------------------------------------
!  ... check confidence levels, if ALL are zero then reset to 100
!-----------------------------------------------------------------------------
   if( all( CONF(:) == 0. ) ) then
     CONF(:) = 100.
   end if

   allocate( mask(nfires),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate mask; error = ',astat
     stop 'Alloc err'
   end if

!-----------------------------------------------------------------------------
!  ... flag fires with confidence level <= 20
!-----------------------------------------------------------------------------
   where( CONF(:) > 20. )
     mask(:) = .true.
   elsewhere
     mask(:) = .false.
   endwhere

   ngoodfires = count( mask(:) )
   nconfgt20  = ngoodfires

   write(*,*) 'finn_fire: total,good fire cnts = ',nfires,ngoodfires
   write(*,*) 'finn_fire: ',nfires-ngoodfires,' low confidence fires'
   write(*,*) 'finn_fire: min,max lct     = ',minval(lct(:)),maxval(lct(:))
   write(*,*) 'finn_fire: min,max globreg = ',minval(globreg(:)),maxval(globreg(:))
   if( any(lct(:) == 0) ) then
     write(*,*) 'finn_fire: ',count(lct(:) == 0),' fires with lct == 0'
     write(*,*) 'finn_fire: index min lct = ',minloc(lct(:))
   endif
   if( any(globreg(:) == 0 ) ) then
     write(*,*) 'finn_fire: ',count(globreg(:) == 0),' fires with globreg == 0'
     write(*,*) 'finn_fire: index min globreg = ',minloc(globreg(:))
   endif

   if( has_intermediate_log ) then
     write(unit=l_unit,fmt=*,iostat=ios) 'The original number of fires is: ',nfires
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to write log file; error = ',ios
       stop 'Write err'
     endif
   endif

!-----------------------------------------------------------------------------
!  ... remove low confidence fires
!-----------------------------------------------------------------------------
remove_fires : &
   if( ngoodfires < nfires ) then
     allocate( iwrk(ngoodfires),rwrk(ngoodfires),cwrk(ngoodfires),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate iwrk ... cwrk; error = ',astat
       stop 'Alloc err'
     end if
     rwrk(:) = pack( lat,mask )
     lat(:ngoodfires) = rwrk(:)
     rwrk(:) = pack( lon,mask )
     lon(:ngoodfires) = rwrk(:)
     rwrk(:) = pack( spix,mask )
     spix(:ngoodfires) = rwrk(:)
     rwrk(:) = pack( tpix,mask )
     tpix(:ngoodfires) = rwrk(:)
     rwrk(:) = pack( CONF,mask )
     CONF(:ngoodfires) = rwrk(:)
     rwrk(:) = pack( tree,mask )
     tree(:ngoodfires) = rwrk(:)
     rwrk(:) = pack( herb,mask )
     herb(:ngoodfires) = rwrk(:)
     rwrk(:) = pack( bare,mask )
     bare(:ngoodfires) = rwrk(:)
     iwrk(:) = pack( lct,mask )
     lct(:ngoodfires) = iwrk(:)
     iwrk(:) = pack( globreg,mask )
     globreg(:ngoodfires) = iwrk(:)
     iwrk(:) = pack( tod,mask )
     tod(:ngoodfires) = iwrk(:)
     cwrk(:) = pack( date,mask )
     date(:ngoodfires) = cwrk(:)

     deallocate( mask, iwrk, rwrk, cwrk )
   endif remove_fires

   write(*,*) 'finn_fire: Finished reading Input file'

   allocate( jd(ngoodfires),totcov(ngoodfires),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate totcov; error = ',astat
     stop 'Alloc err'
   end if

!-----------------------------------------------------------------------------
!  ... just in case limit vcf to zero
!-----------------------------------------------------------------------------
   write(*,*) 'finn_fire: nummisstree = ',count( tree(:ngoodfires) < 0. )
   write(*,*) 'finn_fire: nummissherb = ',count( herb(:ngoodfires) < 0. )
   write(*,*) 'finn_fire: nummissbare = ',count( bare(:ngoodfires) < 0. )

   where( tree(:ngoodfires) < 0. )
     tree(:ngoodfires) = 0.
   endwhere
   where( herb(:ngoodfires) < 0. )
     herb(:ngoodfires) = 0.
   endwhere
   where( bare(:ngoodfires) < 0. )
     bare(:ngoodfires) = 0.
   endwhere

   totcov(:) = tree(:ngoodfires) + herb(:ngoodfires) + bare(:ngoodfires)

   write(*,*) 'finn_fire: nummissvcf = ',count( totcov(:) < 98. )
   write(*,*) 'finn_fire: min,max totcov = ',minval(totcov(:)),maxval(totcov(:))

   nfire1 = ngoodfires

!-----------------------------------------------------------------------------
!  ... set day of year; NOTE, the ntotdays setting depends on the
!                             fire dates all being for the same year
!-----------------------------------------------------------------------------
   do n = 1,nfire1
     jd(n) = doy( date(n) )
     if( n == nfire1 ) then
       if( is_leap_year( date(n) ) ) then
         ntotdays = 366
       else
         ntotdays = 365
       endif
     endif
   end do

   if( speciation_start_date /= ' ' ) then
     speciation_start_day = doy( speciation_start_date )
   endif
   if( speciation_end_date /= ' ' ) then
     speciation_end_day = doy( speciation_end_date )
   endif
   speciation_end_day = min( speciation_end_day,ntotdays )

   write(*,*) 'finn_fire: first date,doy = ',date(1),jd(1)
   write(*,*) 'finn_fire: last  date,doy = ',date(nfire1),jd(nfire1)
   write(*,*) 'finn_fire: min,max doy    = ',minval(jd(:)),maxval(jd(:))
   write(*,*) 'finn_fire: Finished calculating Julian day of year'

!-----------------------------------------------------------------------------
!  ... handle tropical fires
!-----------------------------------------------------------------------------
   ntropics = count( lat(:nfire1) > -30. .and. lat(:nfire1) < 30. )

   allocate( factortrop(nfire1),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate factortrop; error = ',astat
     stop 'Alloc err'
   end if
   factortrop(:) = 1.

   if( has_intermediate_log ) then
     write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires added (because in tropics) is ',ntropics
     if( ios == 0 ) then
       write(unit=l_unit,fmt=*,iostat=ios) 'The new number of fires = ',ngoodfires+ntropics
     endif
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to write log file; error = ',ios
       stop 'Write err'
     endif
   endif

has_tropics : &
   if( ntropics > 0 ) then
     write(*,'('' finn_fire: '',i7,'' tropical fires'')') ntropics
     allocate( mask(nfire1),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate mask; error = ',astat
       stop 'Alloc err'
     end if
     mask(:) = lat(:nfire1) > -30. .and. lat(:nfire1) < 30.
     allocate( rwrk(ntropics),iwrk(ntropics),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate rwrk; error = ',astat
       stop 'Alloc err'
     end if
     rwrk(:) = pack( lat,mask )
     call concat_reals( nfire1, ntropics, lat, rwrk )
     rwrk(:) = pack( lon,mask )
     call concat_reals( nfire1, ntropics, lon, rwrk )
     rwrk(:) = pack( spix,mask )
     call concat_reals( nfire1, ntropics, spix, rwrk )
     rwrk(:) = pack( tpix,mask )
     call concat_reals( nfire1, ntropics, tpix, rwrk )
     rwrk(:) = pack( tree,mask )
     call concat_reals( nfire1, ntropics, tree, rwrk )
     rwrk(:) = pack( herb,mask )
     call concat_reals( nfire1, ntropics, herb, rwrk )
     rwrk(:) = pack( bare,mask )
     call concat_reals( nfire1, ntropics, bare, rwrk )
     rwrk(:) = pack( totcov,mask )
     call concat_reals( nfire1, ntropics, totcov, rwrk )
     rwrk(:) = FAC
     call concat_reals( nfire1, ntropics, factortrop, rwrk )
     iwrk(:) = pack( lct,mask )
     call concat_ints( nfire1, ntropics, lct, iwrk )
     iwrk(:) = pack( jd,mask )
     iwrk(:) = iwrk(:) + 1
     call concat_ints( nfire1, ntropics, jd, iwrk )
     iwrk(:) = pack( tod,mask )
     call concat_ints( nfire1, ntropics, tod, iwrk )
     iwrk(:) = pack( globreg,mask )
     call concat_ints( nfire1, ntropics, globreg, iwrk )
     ngoodfires = nfire1 + ntropics
     deallocate( mask, rwrk, iwrk )
   endif has_tropics

   write(*,'('' finn_fire: '',i7,'' total fires'')') ngoodfires

   allocate( flag(ngoodfires),rwrk(ngoodfires),iwrk(ngoodfires),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate flag ... iwrk; error = ',astat
     stop 'Alloc err'
   end if
   allocate( xearth(ngoodfires),yearth(ngoodfires),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate xearth, yearth; error = ',astat
     stop 'Alloc err'
   end if
   flag(:) = 1

   pi   = 4.*atan(1.)
   dtor = pi/180.

!-----------------------------------------------------------------
!  ... reorder fires by decreasing area
!-----------------------------------------------------------------
   write(*,*) 'finn_fire: min,max spix = ',minval(spix(:)),maxval(spix(:))
   write(*,*) 'finn_fire: min,max tpix = ',minval(tpix(:)),maxval(tpix(:))
   rwrk(:) = tpix(:) * spix(:)
   iwrk(:) = (/ (n,n=1,ngoodfires) /)
   call heapsort( ngoodfires, rwrk, iwrk )
   call reverse_array( ngoodfires, iwrk )
   lat(:)  = lat(iwrk(:))
   lon(:)  = lon(iwrk(:))
   spix(:) = spix(iwrk(:))
   tpix(:) = tpix(iwrk(:))
   tree(:) = tree(iwrk(:))
   herb(:) = herb(iwrk(:))
   bare(:) = bare(iwrk(:))
   lct(:)  = lct(iwrk(:))
   jd(:)   = jd(iwrk(:))
   tod(:)  = tod(iwrk(:))
   totcov(:)     = totcov(iwrk(:))
   globreg(:)    = globreg(iwrk(:))
   factortrop(:) = factortrop(iwrk(:))

   xearth(:) = rearthkm*lon(:)*dtor*cos( lat(:)*dtor )
   yearth(:) = rearthkm*lat(:)*dtor

   write(*,*) 'finn_fire: min,max doy = ',minval(jd(:)),maxval(jd(:))
   write(*,*) 'finn_fire: min,max lon 1st day = ',minval(lon(:),mask=jd(:)==1),maxval(lon(:),mask=jd(:)==1)
   write(*,*) 'finn_fire: min,max lat 1st day = ',minval(lat(:),mask=jd(:)==1),maxval(lat(:),mask=jd(:)==1)

!-----------------------------------------------------------------
!  ... flag spatially overlapping fires
!-----------------------------------------------------------------
   write(*,*) 'finn_fire: removing overlapping fires'
   write(*,*) 'finn_fire: this can take several minutes'
!$omp parallel do private( m ), schedule(dynamic,1)
   do m = 1,ntotdays
     call overlap( m, ngoodfires, ntotdays, jd, flag, &
                   tpix, spix, xearth, yearth, scale )
   end do
!$omp end parallel do

   deallocate( xearth, yearth )

   write(*,*) 'Finished calculating overlapping fires'

!-----------------------------------------------------------------
!  ... sort fires by day of year
!-----------------------------------------------------------------
   rwrk(:) = real( jd(:) )
   iwrk(:) = (/ (n,n=1,ngoodfires) /)
   call heapsort( ngoodfires, rwrk, iwrk )
   lat(:)  = lat(iwrk(:))
   lon(:)  = lon(iwrk(:))
   spix(:) = spix(iwrk(:))
   tpix(:) = tpix(iwrk(:))
   tree(:) = tree(iwrk(:))
   herb(:) = herb(iwrk(:))
   bare(:) = bare(iwrk(:))
   lct(:)  = lct(iwrk(:))
   jd(:)   = jd(iwrk(:))
   tod(:)  = tod(iwrk(:))
   flag(:) = flag(iwrk(:))
   totcov(:)     = totcov(iwrk(:))
   globreg(:)    = globreg(iwrk(:))
   factortrop(:) = factortrop(iwrk(:))

   deallocate( rwrk, iwrk )
   write(*,*) 'finn_fire: min,max doy    = ',minval(jd(:)),maxval(jd(:))
   write(*,*) 'finn_fire: first day popcnt = ',count(jd(:) == 1)

   write(*,*) 'Calculating Emissions'

   overlapct = count( flag(:) == -999 )
   antarc    = count( globreg(:) == 0 )

   write(*,*) 'first day popcnt = ',count(jd(:) == 1)
   write(*,*) 'overlap popcnt   = ',overlapct
   write(*,*) 'antarc popcnt    = ',antarc

!-----------------------------------------------------------------
!  ... calculate and write out emissions
!-----------------------------------------------------------------

emiss_loop : &
   do j = 1,ngoodfires
     if( flag(j) == -999 .or. globreg(j) == 0 ) then
       cycle emiss_loop
     elseif( lct(j) < 1 .or. lct(j) > 16 .or. lct(j) == 15 ) then
       lct0 = lct0 + 1
       cycle emiss_loop
     endif
     wrk_real = totcov(j)
     if( wrk_real /= 100. ) then
       if( wrk_real >= 1. .and. wrk_real < 240.) then
         wrk_reali = 100./wrk_real
         tree(j) = tree(j)*wrk_reali
         herb(j) = herb(j)*wrk_reali
         bare(j) = bare(j)*wrk_reali
         totcov(j) = tree(j) + herb(j) + bare(j)
         if( wrk_real < 50. ) then
           vcflt50 = vcflt50 + 1
         else
           vcfcount = vcfcount + 1
         endif
       endif
     endif
     if( totcov(j) < 1. .or. totcov(j) >= 240. .or. bare(j) == 100. ) then
       allbare = allbare + 1
       if( lct(j) >= 15 ) then
         cycle emiss_loop
       elseif( lct(j) <= 5 ) then
         tree(j) = 60. ; herb(j) = 40. ; bare(j) = 0.
       elseif( (lct(j) >= 6 .and. lct(j) <= 8) .or. lct(j) == 11 .or. lct(j) == 14 ) then
         tree(j) = 50. ; herb(j) = 50. ; bare(j) = 0.
       else
         tree(j) = 20. ; herb(j) = 80. ; bare(j) = 0.
       endif
     endif
     wrk_int = lct(j)
     wrk_lat = lat(j)
     in_tropics = wrk_lat > -30. .and. wrk_lat < 30.
     if( wrk_int > 5 ) then
       select case( wrk_int )
         case( 6:8 )
           genveg = 2
         case( 9:11 )
           genveg = 1
         case( 12 )
           genveg = 9
         case( 13 )
           urbnum = urbnum + 1
           if( tree(j) < 40. ) then
             genveg = 1
             lct(j) = 10
           elseif( tree(j) >= 40. .and. tree(j) < 60. ) then
             genveg = 2
             lct(j) = 8
           elseif( tree(j) >= 60. ) then
             if( wrk_lat > 50. ) then
               genveg = 5
               lct(j) = 1
             else
               if( in_tropics ) then
                 genveg = 3
               else
                 genveg = 4
               endif
               lct(j) = 5
             endif
           endif
         case( 14,16 )
           genveg = 1
       end select
     else
       select case( wrk_int )
         case( 1,3 )
           if( wrk_lat > 50. ) then
             genveg = 5
           else
             genveg = 4
           endif
         case( 2 )
           genveg = 3
         case( 4 )
           genveg = 4
         case( 5 )
           if( in_tropics ) then
             genveg = 3
           elseif( wrk_lat > 50. ) then
             genveg = 5
           else
             genveg = 4
           endif
       end select
     endif

     if( genveg == 1 .or. genveg == 9 ) then
       area = .75
     else
       area = 1.
     endif
     area = area*(100. - bare(j))*.01
     reg = globreg(j)
     if( reg < 1 ) then
       if( has_intermediate_log ) then
         write(unit=l_unit,fmt=*,iostat=ios) 'Fire number: ',j, &
          ' removed. Something is WRONG with global regions and fuel loads. Globreg =', globreg(j)
         if( ios /= 0 ) then
           write(*,*) 'finn_fire: Failed to write log file; error = ',ios
           stop 'Write err'
         endif
       endif
       cycle emiss_loop
     endif
     if( genveg == 0 ) then
       genveg0 = genveg0 + 1
       if( has_intermediate_log ) then
         write(unit=l_unit,fmt=*,iostat=ios) 'Fire number: ',j, &
          ' removed. Something is WRONG with generic vegetation. genveg = 0'
         if( ios /= 0 ) then
           write(*,*) 'finn_fire: Failed to write log file; error = ',ios
           stop 'Write err'
         endif
       endif
       cycle emiss_loop
     endif
     select case( genveg )
       case ( 1 )
         bmass1 = fuel_load(reg)%grfuel
       case ( 2 )
         bmass1 = fuel_load(reg)%wsfuel
       case ( 3 )
         bmass1 = fuel_load(reg)%tffuel
       case ( 4 )
         bmass1 = fuel_load(reg)%tefuel
       case ( 5 )
         if( globreg(j) /= 11 ) then
           bmass1 = fuel_load(reg)%bffuel
         else
           bmass1 = fuel_load(reg)%tefuel
         endif
       case ( 9 )
         if( (lon(j) <= -47.323 .and. lon(j) >= -49.156) .and. &
             (lat(j) <= -20.356 .and. lat(j) >= -22.708) ) then
           bmass1 = 1100.
         else
           bmass1 = 500.
         endif
       case default
         bmass0 = bmass0 + 1
         if( has_intermediate_log ) then
           write(unit=l_unit,fmt=*,iostat=ios) 'Fire number:',j,' removed. bmass not assigned'
           if( ios == 0 ) then
             write(unit=l_unit,fmt=*,iostat=ios) 'genveg =', genveg, ' and globreg = ', globreg(j), ' and reg = ', reg
           endif
           if( ios /= 0 ) then
             write(*,*) 'finn_fire: Failed to write log file; error = ',ios
             stop 'Write err'
           endif
         endif
         cycle emiss_loop
     end select

     pctherb  = .01*herb(j)
     pcttree  = .01*tree(j)
     herbbm   = fuel_load(reg)%grfuel
     coarsebm = bmass1
     if( tree(j) > 60. ) then
       CF1 = .3 ; CF3 = .9
       bmass = (pctherb + pcttree)*herbbm*CF3 + pcttree*coarsebm*CF1
     elseif( tree(j) > 40. ) then
       CF1 = .3
       CF3 = exp( -.013*pcttree )
       bmass = (pctherb + pcttree)*herbbm*CF3 + pcttree*coarsebm*CF1
     else
       CF3 = .98
       bmass = (pctherb + pcttree)*herbbm*CF3
     endif
     ndx = lct(j)
     if( lct(j) == 14 ) then
       ndx = 13
     elseif( lct(j) == 16 ) then
       ndx = 14
     endif
     wrk_real = area*bmass*factortrop(j)
     area  = area*1.e6
     bmass = bmass*.001
     CO2 = emis_factor(ndx)%CO2EF*wrk_real
     CO  = emis_factor(ndx)%COEF*wrk_real
     CH4 = emis_factor(ndx)%CH4EF*wrk_real
     NMHC = emis_factor(ndx)%NMHCEF*wrk_real
     NMOC = emis_factor(ndx)%NMOCEF*wrk_real
     H2  = emis_factor(ndx)%H2EF*wrk_real
     NOX = emis_factor(ndx)%NOXEF*wrk_real
     NO  = emis_factor(ndx)%NOEF*wrk_real
     NO2 = emis_factor(ndx)%NO2EF*wrk_real
     SO2 = emis_factor(ndx)%SO2EF*wrk_real
     TPM = emis_factor(ndx)%TPMEF*wrk_real
     TPC = emis_factor(ndx)%TCEF*wrk_real
     PM10 = emis_factor(ndx)%PM10EF*wrk_real
     PM25 = emis_factor(ndx)%PM25EF*wrk_real
     OC  = emis_factor(ndx)%OCEF*wrk_real
     BC  = emis_factor(ndx)%BCEF*wrk_real
     NH3 = emis_factor(ndx)%NH3EF*wrk_real

     if( has_intermediate_output ) then
       write(unit=o_unit,fmt=trim(i_frmt),iostat=ios) lon(j),lat(j),jd(j),tod(j),lct(j),genveg,globreg(j),int(tree(j)), &
         int(herb(j)),int(bare(j)),area,bmass,CO2,CO,CH4,H2,NOX,NO,NO2,NH3,SO2,NMHC,NMOC,PM25,TPM,OC,BC,TPC,PM10,Factortrop(j)
       if( ios /= 0 ) then
         write(*,*) 'finn_fire: failed to write intermediate output file; error = ',ios
         stop 'Write err'
       endif
     endif

!--------------------------
! LKE 7/31/2013
!  apparently have to manually match the order of these arrays with the speciation input file
!  and the output file header
!  spec_ndx=1: MOZ4, 2:SAPRC99, 3:GEOSchem
!  'MOZ4': MOZ4_chem_speciation_table.csv
!  'DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,H2,NO,NO2,SO2,NH3,CH4,' // &
!                                                 'NMOC,BIGALD,BIGALK,BIGENE,C10H16,C2H4,C2H5OH,C2H6,C3H6,C3H8,CH2O,' // &
!                                                 'CH3CHO,CH3CN,CH3COCH3,CH3COCHO,CH3COOH,CH3OH,CRESOL,GLYALD,HCN,' // &
!                                                 'HYAC,ISOP,MACR,MEK,MVK,TOLUENE,HCOOH,C2H2,OC,BC,PM10,PM25'
! 'SAPRC99': 'DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,NO,NO2,SO2,NH3,CH4,VOC,ACET,' // &
!                                                 'ALK1,ALK2,ALK3,ALK4,ALK5,ARO1,ARO2,BALD,CCHO,CCO_OH,ETHENE,HCHO,HCN,' // &
!                                                 'HCOOH,HONO,ISOPRENE,MEK,MEOH,METHACRO,MGLY,MVK,OLE1,OLE2,PHEN,PROD2,' // &
!                                                 'RCHO,TRP1,OC,BC,PM10,PM25'
! 'GEOSCHEM': Geos_chem_speciation_table.csv (7/31/2013)
!  Output order:  'DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,NO,NO2,SO2,NH3,CH4,ACET,ALD2,' // &
!                                                 'ALK4,BENZ,C2H2,C2H4,C2H6,C3H8,CH2O,GLYC,GLYX,HAC,MEK,MGLY,PRPE,' // &
!                                                 'TOLU,XYLE,OC,BC,PM25'
!--------------------------
has_speciated_emissions : &
     if( has_speciated_output ) then
       speciated_emissions(1) = CO2*1.e3/CO2_mw
       speciated_emissions(2) = CO*1.e3/CO_mw
       ndx = genveg
       if( genveg == 9 ) then
         ndx = 6
       endif
       if( spec_ndx > 1 ) then
         speciated_emissions(3) = NO*1.e3/NO_mw
         speciated_emissions(4) = NO2*1.e3/NO2_mw
         speciated_emissions(5) = SO2*1.e3/SO2_mw
         speciated_emissions(6) = NH3*1.e3/NH3_mw
         speciated_emissions(7) = CH4*1.e3/CH4_mw
         if( spec_ndx == 3 ) then
           speciated_emissions(3)     = speciated_emissions(3) + VOC*speciation(18,ndx)    ! NO+HONO (called NO in spec.file)
           speciated_emissions(8:24)  = VOC*speciation(:17,ndx)                            ! ACET ... XYLE
           speciated_emissions(25:27) = (/ OC, BC, PM25 /)
         else
           speciated_emissions(8)     = VOC
           speciated_emissions(9:36)  = VOC*speciation(:28,ndx)      ! ACET ...  TRP1
           speciated_emissions(37:40) = (/ OC, BC, PM10, PM25 /)
         endif
       elseif( spec_ndx == 1 ) then
         speciated_emissions(3) = H2*1.e3/H2_mw
         speciated_emissions(4) = NO*1.e3/NO_mw + VOC*speciation(25,ndx)  !NO+HONO
         speciated_emissions(5) = NO2*1.e3/NO2_mw
         speciated_emissions(6) = SO2*1.e3/SO2_mw
         speciated_emissions(7) = NH3*1.e3/NH3_mw
         speciated_emissions(8) = CH4*1.e3/CH4_mw
         speciated_emissions(9) = VOC
         speciated_emissions(10:33) = VOC*speciation(:24,ndx)        ! BIGALD ... MVK
         speciated_emissions(34:36) = VOC*speciation(26:28,ndx)      ! TOLUENE,HCOOH,C2H2
         speciated_emissions(37:40) = (/ OC, BC, PM10, PM25 /)
       endif
       if( speciation_start_day <= jd(j) .and. jd(j) <= speciation_end_day ) then
         write(s_unit,fmt=trim(s_frmt),iostat=ios) jd(j),tod(j),genveg,lat(j),lon(j),area,speciated_emissions(:)
         if( ios /= 0 ) then
           write(*,*) 'finn_fire: failed to write speciated output file; error = ',ios
           stop 'Write err'
         endif
       endif
     endif has_speciated_emissions

     if( has_intermediate_log ) then
       bmassburn = bmass*area
       bmasstotal = bmasstotal + bmassburn
       select case( genveg )
       case( 1 )
         totgras = totgras + bmassburn
         totgrasarea = totgrasarea + area
       case( 2 )
         totshrub = totshrub + bmassburn
         totshrubarea = totshrubarea + area
       case( 3 )
         tottrop = tottrop + bmassburn
         tottroparea = tottroparea + area
       case( 4 )
         tottemp = tottemp + bmassburn
         tottemparea = tottemparea + area
       case( 5 )
         totbor = totbor + bmassburn
         totborarea = totborarea + area
       case( 9 )
         totcrop = totcrop + bmassburn
         totcroparea = totcroparea + area
       end select
!-----------------------------------------------------------------
!  ... global sums
!-----------------------------------------------------------------
       co2total  = co2total + co2
       cototal   = cototal + co
       ch4total  = ch4total + ch4
       nmhctotal = nmhctotal + nmhc
       nmoctotal = nmoctotal + nmoc
       h2total   = h2total + h2
       noxtotal  = noxtotal + nox
       nototal   = nototal + no
       no2total  = no2total + no2
       so2total  = so2total + so2
       pm10total = pm10total + pm10
       pm25total = pm25total + pm25
       tpmtotal  = tpmtotal + tpm
       tpctotal  = tpctotal + tpc
       octotal   = octotal + oc
       bctotal   = bctotal + bc
       nh3total  = nh3total + nh3
       areatotal = areatotal + area
     endif
   end do emiss_loop

   write(*,*) 'Finished calculating Emissions'

!-----------------------------------------------------------------
!  ... log file output
!-----------------------------------------------------------------
   if( has_intermediate_log ) then
     call write_intermediate_log
   endif

   deallocate( lat, lon, spix, tpix, tree, herb, bare, lct, &
               globreg, date, tod, CONF, totcov, jd, flag, factortrop, &
               fuel_load, emis_factor )
   if( has_speciated_output ) then
     deallocate( speciation, speciated_emissions )
   endif

!-----------------------------------------------------------------
!  ... close units
!-----------------------------------------------------------------
   if( has_intermediate_output ) then
     close( o_unit )
   endif
   if( has_intermediate_log ) then
     close( l_unit )
   endif
   if( has_speciated_output ) then
     close( s_unit )
   endif
   if( has_speciated_log ) then
     close( c_unit )
   endif

   CONTAINS

   integer function get_file_size( unitno )
!-----------------------------------------------------------------
!  get file size
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  dummy arguments
!-----------------------------------------------------------------
   integer, intent(in) :: unitno
!-----------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------
   integer :: nl, istat
   character :: c

   nl = 0
   do
     read(unitno,fmt=*,iostat=istat) c
     if( istat == 0 ) then
       nl = nl + 1
     else
       exit
     endif
   end do
   
   rewind( unitno )
   get_file_size = nl

   end function get_file_size

   subroutine gettokens( string, ls, delim, maxlen, tokens, &
                         toklen, maxtok, tokcnt )
!-----------------------------------------------------------------------     
!     Input arguments:
!        string - character string to crack into tokens
!        ls     - length of string
!        delim  - token delimiter character
!        maxlen - maximum length of any single token
!        maxtok - maximum number of tokens
!     Output arguments:
!        tokcnt - number of actual tokens
!                 < 0 => hit maxtok before end of string
!                 = 0 => error in input string
!        toklen - array containing length of each token
!        tokens - character array of tokens
!-----------------------------------------------------------------------     

   integer, intent(in)  ::  ls, maxlen, maxtok
   integer, intent(out) ::  tokcnt
   integer, intent(out) ::  toklen(*)
      
   character(len=*), intent(in)  :: string
   character(len=*), intent(out) :: tokens(*)
   character(len=1), intent(in)  :: delim
      
!-----------------------------------------------------------------------     
!  ... local variables
!-----------------------------------------------------------------------     
   integer  ::   marker, i, length, endpos

   tokcnt = 0
   marker = 1
character_loop : &
   do i = 1,ls
have_delimiter : &
     if( string(i:i) == delim .or. i == ls ) then
       if( i == ls ) then
         if( string(i:i) == delim ) then
           tokcnt = 0
           exit character_loop
         end if
         length = i - marker + 1
         endpos = i
       else
         length = i - marker
         endpos = i - 1
       end if
       if( length < 1 .or. length > maxlen ) then
         tokcnt = 0
         exit character_loop
       end if
       tokcnt = tokcnt + 1
       if( tokcnt > maxtok ) then
         tokcnt = -tokcnt
         exit character_loop
       end if
       tokens(tokcnt) = ' '
       tokens(tokcnt)(:length) = string(marker:endpos)
       toklen(tokcnt) = length
       marker = i + 1
     endif have_delimiter
   end do character_loop
      
   end subroutine gettokens

   subroutine heapsort( n, ra, ndx )
!-----------------------------------------------------------------------     
!  ... heap sort
!-----------------------------------------------------------------------     
 
   integer, intent(in) :: n
   integer, intent(inout) :: ndx(n)
   real, intent(inout)    :: ra(n)

   integer :: i, ir, j, l
   integer :: nh
   real    :: rra

   if( n >= 2 ) then
     l  = n/2 + 1
     ir = n
outer_loop : &
     do
       if( l > 1 ) then
         l = l - 1
         rra = ra(l)
         nh = ndx(l)
       else
         rra = ra(ir)
         nh = ndx(ir)
         ra(ir) = ra(1)
         ndx(ir) = ndx(1)
         ir = ir - 1
         if( ir == 1 ) then
           ra(1) = rra
           ndx(1) = nh
           exit outer_loop
         endif
       endif
       i = l
       j = l + l
inner_loop : &
       do
         if( j <= ir ) then
           if( j < ir .and. ra(j) < ra(j+1) ) then
             j = j + 1
           endif
           if( rra < ra(j) ) then
             ra(i) = ra(j)
             ndx(i) = ndx(j)
             i = j
             j = 2*j
           else
             j = ir + 1
           endif
         else
           exit inner_loop
         endif
         ra(i) = rra
         ndx(i) = nh
       end do inner_loop
     end do outer_loop
   endif

   end subroutine heapsort

   logical function is_leap_year( date )
!-----------------------------------------------------------------------------
!   ... determine if date is for a leap year
!-----------------------------------------------------------------------------

   character(len=*), intent(in) :: date

   integer :: my1, my2, my3, yr

   read(date(7:10),fmt=*,iostat=ios) yr
   if( ios /= 0 ) then
     write(*,*) 'is_leap_year: failed to read year from date; error = ',ios
     stop 'Rd_err'
   endif

   my1 = mod(yr,4)
   my2 = mod(yr,100)
   my3 = mod(yr,400)

   is_leap_year = (MY1 == 0 .AND. MY2 /= 0) .OR. MY3 == 0

   end function is_leap_year

   subroutine concat_ints( n1, n2, vec1, vec2 )
!-----------------------------------------------------------------------
!  ... concatenate vec2 to vec1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ... dummy arguments
!-----------------------------------------------------------------------
   integer, intent(in) :: n1, n2
   integer, intent(in) :: vec2(:)
   integer, pointer    :: vec1(:)
!-----------------------------------------------------------------------
!  ... local variables
!-----------------------------------------------------------------------
   integer :: astat
   integer, allocatable :: wrk(:)

   allocate( wrk(n1),stat=astat)
   if( astat /= 0 ) then
     write(*,*) 'concat_ints: failed to allocate wrk space; error = ',astat
     stop 'Alloc err'
   endif

   wrk(:) = vec1(:)
   deallocate( vec1,stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'concat_ints: failed to deallocate vec1; error = ',astat
     stop 'Dealloc err'
   endif

   allocate( vec1(n1+n2),stat=astat)
   if( astat /= 0 ) then
     write(*,*) 'concat_ints: failed to reallocate vec1; error = ',astat
     stop 'Alloc err'
   endif

   vec1(:n1) = wrk(:)
   vec1(n1+1:n1+n2) = vec2(:)

   deallocate( wrk )

   end subroutine concat_ints

   subroutine concat_reals( n1, n2, vec1, vec2 )
!-----------------------------------------------------------------------
!  ... concatenate vec2 to vec1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ... dummy arguments
!-----------------------------------------------------------------------
   integer, intent(in) :: n1, n2
   real, intent(in) :: vec2(:)
   real, pointer    :: vec1(:)
!-----------------------------------------------------------------------
!  ... local variables
!-----------------------------------------------------------------------
   integer :: astat
   real, allocatable :: wrk(:)

   allocate( wrk(n1),stat=astat)
   if( astat /= 0 ) then
     write(*,*) 'concat_reals: failed to allocate wrk space; error = ',astat
     stop 'Alloc err'
   endif

   wrk(:) = vec1(:)
   deallocate( vec1,stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'concat_reals: failed to deallocate vec1; error = ',astat
     stop 'Dealloc err'
   endif

   allocate( vec1(n1+n2),stat=astat)
   if( astat /= 0 ) then
     write(*,*) 'concat_reals: failed to reallocate vec1; error = ',astat
     stop 'Alloc err'
   endif

   vec1(:n1) = wrk(:)
   vec1(n1+1:n1+n2) = vec2(:)

   deallocate( wrk )

   end subroutine concat_reals

   integer function doy( date )
!-----------------------------------------------------------------------
!  ... Compute day of year
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ... dummy args
!-----------------------------------------------------------------------
   character(len=*), intent(in) :: date

   integer :: yr, mon, cday

   integer, save :: jdbase(12) = &
         (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)
   integer :: jdcon(12)

   read(date(7:10),fmt=*,iostat=ios) yr
   if( ios == 0 ) then
     read(date(1:2),fmt=*,iostat=ios) mnth
     if( ios == 0 ) then
       read(date(4:5),fmt=*,iostat=ios) day
     endif
   endif
   if( ios /= 0 ) then
     write(*,*) 'doy: character date to yr,mnth,day conversion failed; error = ',ios
     stop 'Conversion err'
   endif

   jdcon(:) = jdbase(:)
   if( is_leap_year( date ) ) then
     jdcon(3:) = jdcon(3:) + 1
   endif

   doy = jdcon(mnth) + day

   end function doy

   subroutine overlap( m, ngf, ntotdays, jd, flag, &
                       tpix, spix, xearth, yearth, scale )
!-----------------------------------------------------------------------
!  ... identify and flag overlapping fires
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ... dummy args
!-----------------------------------------------------------------------
   integer, intent(in) :: m
   integer, intent(in) :: ngf
   integer, intent(in) :: ntotdays
   integer, pointer    :: jd(:)
   integer, intent(inout)       :: flag(:)
   real, intent(in)    :: scale
   real, intent(in)    :: xearth(:)
   real, intent(in)    :: yearth(:)
   real, pointer       :: tpix(:)
   real, pointer       :: spix(:)

!-----------------------------------------------------------------------
!  ... local variables
!-----------------------------------------------------------------------
   integer :: j, n, n1, ntoday, novrlap
   integer :: astat, cntr
   integer, allocatable :: today_ndx(:), ovrlap_ndx(:)
   real    :: xtrack, atrack
   real    :: dxdumm, dxdump, dydumm, dydump

   ntoday = count( jd(:) == m )
is_today : &
   if( ntoday > 0 ) then
       allocate( today_ndx(ntoday),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'finn_fire: failed to allocate today_ndx; error = ',astat
         stop 'Alloc err'
       end if
       today_ndx(:) = pack( (/ (j,j=1,ngf) /),mask=(jd(:) == m) )
       cntr = max( 1,count(flag(today_ndx(:)) == -999) )
match_loop : &
       do n = 1,ntoday-1
         n1 = n + 1
         j = today_ndx(n)
         if( flag(j) /= -999 ) then
           xtrack = tpix(j)*scale
           atrack = spix(j)*scale
           dxdump = xearth(j) + .5*xtrack
           dydump = yearth(j) + .5*atrack
           dxdumm = xearth(j) - .5*xtrack
           dydumm = yearth(j) - .5*atrack
           novrlap = count( flag(today_ndx(n1:ntoday)) == 1 .and. &
                            yearth(today_ndx(n1:ntoday)) >= dydumm .and. &
                            yearth(today_ndx(n1:ntoday)) <= dydump .and. &
                            xearth(today_ndx(n1:ntoday)) >= dxdumm .and. &
                            xearth(today_ndx(n1:ntoday)) <= dxdump )
has_ovrlap : &
           if( novrlap > 0 ) then
             allocate( ovrlap_ndx(novrlap),stat=astat )
             if( astat /= 0 ) then
               write(*,*) 'finn_fire: failed to allocate ovrlap_ndx; error = ',astat
               stop 'Alloc err'
             end if
             ovrlap_ndx(:) = pack( today_ndx(n1:ntoday),mask= &
                                   flag(today_ndx(n1:ntoday)) == 1 .and. &
                                   yearth(today_ndx(n1:ntoday)) >= dydumm .and. &
                                   yearth(today_ndx(n1:ntoday)) <= dydump .and. &
                                   xearth(today_ndx(n1:ntoday)) >= dxdumm .and. &
                                   xearth(today_ndx(n1:ntoday)) <= dxdump )
               flag(ovrlap_ndx(:)) = -999
             deallocate( ovrlap_ndx )
           endif has_ovrlap
         endif
       end do match_loop
!      write(*,'('' finn_fire: on day '',i3,'' examined '',i5,'' fires; '',i5,'' ovrlapping fires'')') &
!           m,ntoday,count(flag(today_ndx(:)) == -999) - cntr + 1
       deallocate( today_ndx )
   endif is_today

   end subroutine overlap

   subroutine read_speciation_file
!-----------------------------------------------------------------------
!  ... read the speciation file
!-----------------------------------------------------------------------

     open( unit=i_unit,file=trim(speciation_infile_filespec),iostat=ios )
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to open ',trim(speciation_infile_filespec),'; error = ',ios
       stop 'Open err'
     end if
     nlines = get_file_size( unitno=i_unit ) - 1
     if( nlines > 0 ) then
       select case( trim(speciation_case) )
         case ( 'MOZ4' )
           n = 40
         case ( 'SAPRC99' )
           n = 40
         case ( 'GEOSCHEM' )
           n = 27
       end select
       allocate( speciation(nlines,6),speciated_emissions(n),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'finn_fire: failed to allocate speciation type; error = ',astat
         stop 'Alloc err'
       end if
     else
       write(*,*) 'finn_fire: No data in ',trim(speciation_infile_filespec)
       stop 'Data err'
     endif
   
     read(unit=i_unit,fmt=*,iostat=ios) buffer
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to read ',trim(speciation_infile_filespec),'; error = ',ios
       stop 'Read err'
     end if
     ntokens = 7
     allocate( tokens(ntokens),toklen(ntokens),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate speciation token arrays; error = ',astat
       stop 'Alloc err'
     end if
     do n = 1,nlines
       read(unit=i_unit,fmt='(a)',iostat=ios) buffer
       if( ios /= 0 ) then
         write(*,*) 'finn_fire: failed to read ',trim(speciation_infile_filespec),' header; error = ',ios
         stop 'Read err'
       end if
       slen = len_trim(buffer)
       if( iachar(buffer(slen:slen)) == 13 ) then
         slen = slen - 1
       endif
       if( buffer(slen:slen) == ',' ) then
         slen = slen - 1
       endif
       call gettokens( buffer, slen, comma, 32, tokens, &
                       toklen, ntokens, tokcnt )
       if( tokcnt /= ntokens ) then
         write(*,*) 'finn_fire: speciation file should have ',ntokens,' fields per line'
         write(*,*) 'finn_fire: but has ',tokcnt,' fields from input file'
         stop 'Data err'
       endif
       do m = 2,tokcnt
         read(tokens(m),fmt=*,iostat=ios) speciation(n,m-1)
         if( ios /= 0 ) then
           write(*,*) 'finn_fire: failed to read ',trim(speciation_infile_filespec),'; error = ',ios
           stop 'Read err'
         end if
       end do
     end do

     close( i_unit )
     deallocate( tokens, toklen )

   end subroutine read_speciation_file

   subroutine read_fuel_load_file
!-----------------------------------------------------------------------
!  ... read the fuel load file
!-----------------------------------------------------------------------

   open( unit=i_unit,file=trim(fuel_load_filespec),iostat=ios )
   if( ios /= 0 ) then
     write(*,*) 'finn_fire: failed to open ',trim(fuel_load_filespec),'; error = ',ios
     stop 'Open err'
   end if
   nlines = get_file_size( unitno=i_unit ) - 1
   if( nlines > 0 ) then
     allocate( fuel_load(nlines),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate fuel_load type; error = ',astat
       stop 'Alloc err'
     end if
   else
     write(*,*) 'finn_fire: No data in ',trim(fuel_load_filespec)
     stop 'Data err'
   endif
   
   read(unit=i_unit,fmt=*,iostat=ios) buffer
   if( ios /= 0 ) then
     write(*,*) 'finn_fire: failed to read ',trim(fuel_load_filespec),'; error = ',ios
     stop 'Read err'
   end if
   ntokens = 6
   allocate( tokens(ntokens),toklen(ntokens),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate fuel_load token arrays; error = ',astat
     stop 'Alloc err'
   end if
   do n = 1,nlines
     read(unit=i_unit,fmt='(a)',iostat=ios) buffer
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to read ',trim(fuel_load_filespec),' header; error = ',ios
       stop 'Read err'
     end if
     slen = len_trim(buffer)
     if( iachar(buffer(slen:slen)) == 13 ) then
       slen = slen - 1
     endif
     if( buffer(slen:slen) == ',' ) then
       slen = slen - 1
     endif
     call gettokens( buffer, slen, comma, 32, tokens, &
                     toklen, ntokens, tokcnt )
     if( tokcnt /= ntokens ) then
       write(*,*) 'finn_fire: fuel load file should have ',ntokens,' fields per line'
       write(*,*) 'finn_fire: but has ',tokcnt,' fields from input file'
       stop 'Data err'
     endif
     do m = 1,tokcnt
       if( m == 1 ) then
         read(tokens(m),fmt=*,iostat=ios) wrk_int
       else
         read(tokens(m),fmt=*,iostat=ios) wrk_real
       end if
       if( ios /= 0 ) then
         write(*,*) 'finn_fire: failed to read ',trim(fuel_load_filespec),'; error = ',ios
         stop 'Read err'
       end if
       select case( m )
         case( 1 )
           fuel_load(n)%globreg2 = wrk_int
         case( 2 )
           fuel_load(n)%tffuel = wrk_real
         case( 3 )
           fuel_load(n)%tefuel = wrk_real
         case( 4 )
           fuel_load(n)%bffuel = wrk_real
         case( 5 )
           fuel_load(n)%wsfuel = wrk_real
         case( 6 )
           fuel_load(n)%grfuel = wrk_real
       end select
     end do
   end do

   close( i_unit )

   deallocate( tokens, toklen )

   end subroutine read_fuel_load_file

   subroutine read_emis_factors_file
!-----------------------------------------------------------------------
!  ... read the emission factors file
!-----------------------------------------------------------------------

   open( unit=i_unit,file=trim(emis_factor_filespec),iostat=ios )
   if( ios /= 0 ) then
     write(*,*) 'finn_fire: failed to open ',trim(emis_factor_filespec),'; error = ',ios
     stop 'Open err'
   end if
   nlines = get_file_size( unitno=i_unit ) - 1
   if( nlines > 0 ) then
     allocate( emis_factor(nlines),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate emis_factor type; error = ',astat
       stop 'Alloc err'
     end if
   else
     write(*,*) 'finn_fire: No data in ',trim(emis_factor_filespec)
     stop 'Data err'
   endif
   
   read(unit=i_unit,fmt=*,iostat=ios) buffer
   if( ios /= 0 ) then
     write(*,*) 'finn_fire: failed to read ',trim(emis_factor_filespec),' header; error = ',ios
     stop 'Read err'
   end if
   ntokens = 20
   allocate( tokens(ntokens),toklen(ntokens),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate emis_factor token arrays; error = ',astat
     stop 'Alloc err'
   end if
   do n = 1,nlines
     read(unit=i_unit,fmt='(a)',iostat=ios) buffer
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to read ',trim(emis_factor_filespec),'; error = ',ios
       stop 'Read err'
     end if
     slen = len_trim(buffer)
     if( iachar(buffer(slen:slen)) == 13 ) then
       slen = slen - 1
     endif
     if( buffer(slen:slen) == ',' ) then
       slen = slen - 1
     endif
     call gettokens( buffer, slen, comma, 32, tokens, &
                     toklen, ntokens, tokcnt )
     if( tokcnt /= ntokens ) then
       write(*,*) 'finn_fire: emis factor file should have ',ntokens,' fields per line'
       write(*,*) 'finn_fire: but has ',tokcnt,' fields from input file'
       stop 'Data err'
     endif
     do m = 1,tokcnt
       if( m < 3 ) then
         read(tokens(m),fmt=*,iostat=ios) wrk_int
       elseif( m > 3 ) then
         read(tokens(m),fmt=*,iostat=ios) wrk_real
       else
         cycle
       end if
       if( ios /= 0 ) then
         write(*,*) 'finn_fire: failed to read ',trim(emis_factor_filespec),'; error = ',ios
         stop 'Read err'
       end if

       select case( m )
         case( 1 )
           emis_factor(n)%lctemis = wrk_int
         case( 2 )
           emis_factor(n)%vegemis = wrk_int
         case( 4 )
           emis_factor(n)%co2ef = wrk_real
         case( 5 )
           emis_factor(n)%coef = wrk_real
         case( 6 )
           emis_factor(n)%ch4ef = wrk_real
         case( 7 )
           emis_factor(n)%nmocef = wrk_real
         case( 8 )
           emis_factor(n)%h2ef = wrk_real
         case( 9 )
           emis_factor(n)%noxef = wrk_real
         case( 10 )
           emis_factor(n)%so2ef = wrk_real
         case( 11 )
           emis_factor(n)%pm25ef = wrk_real
         case( 12 )
           emis_factor(n)%tpmef = wrk_real
         case( 13 )
           emis_factor(n)%tcef = wrk_real
         case( 14 )
           emis_factor(n)%ocef = wrk_real
         case( 15 )
           emis_factor(n)%bcef = wrk_real
         case( 16 )
           emis_factor(n)%nh3ef = wrk_real
         case( 17 )
           emis_factor(n)%noef = wrk_real
         case( 18 )
           emis_factor(n)%no2ef = wrk_real
         case( 19 )
           emis_factor(n)%nmhcef = wrk_real
         case( 20 )
           emis_factor(n)%pm10ef = wrk_real
       end select
     end do
   end do

   close( i_unit )

   deallocate( tokens, toklen )

   end subroutine read_emis_factors_file

   subroutine read_raw_fire_file
!-----------------------------------------------------------------------
!  ... read the raw fire input file
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ... local variables
!-----------------------------------------------------------------------
   integer :: hour,min
   character(len=2) :: digits

   open( unit=i_unit,file=trim(raw_infile_filespec),iostat=ios )
   if( ios /= 0 ) then
     write(*,*) 'finn_fire: failed to open ',trim(raw_infile_filespec),'; error = ',ios
     stop 'Open err'
   else
     write(*,*) 'finn_fire: opened in file: ',trim(raw_infile_filespec)
   end if

!-----------------------------------------------------------------------
!  ... new format has a header, reduce count by one
!-----------------------------------------------------------------------
   nfires = get_file_size( unitno=i_unit ) - 1
   if( nfires > 0 ) then
     allocate( lat(nfires),lon(nfires),spix(nfires),tpix(nfires), &
               tree(nfires),herb(nfires),bare(nfires),lct(nfires), &
               globreg(nfires),date(nfires),tod(nfires),CONF(nfires),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate lat ... CONF; error = ',astat
       stop 'Alloc err'
     end if
   else
     write(*,*) 'finn_fire: No data in ',trim(raw_infile_filespec)
     stop 'Data err'
   endif
   
   ntokens = 12
   allocate( tokens(ntokens),toklen(ntokens),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate infile token arrays; error = ',astat
     stop 'Alloc err'
   end if

!-----------------------------------------------------------------------
!  ... read header
!-----------------------------------------------------------------------
   read(unit=i_unit,fmt='(a)',iostat=ios) buffer
   if( ios /= 0 ) then
     write(*,*) 'finn_fire: line # = ',n
     write(*,*) 'finn_fire: failed to read ',trim(raw_infile_filespec),' header; error = ',ios
     stop 'Read err'
   end if

!-----------------------------------------------------------------------
!  ... read raw data
!-----------------------------------------------------------------------
fire_read_loop : &
   do n = 1,nfires
     read(unit=i_unit,fmt='(a)',iostat=ios) buffer
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: line # = ',n
       write(*,*) 'finn_fire: failed to read ',trim(raw_infile_filespec),' header; error = ',ios
       stop 'Read err'
     end if
     slen = len_trim(buffer)
     if( iachar(buffer(slen:slen)) == 13 ) then
       slen = slen - 1
     endif
     if( buffer(slen:slen) == ',' ) then
       slen = slen - 1
     endif
     call gettokens( buffer, slen, comma, 32, tokens, &
                     toklen, ntokens, tokcnt )
     if( tokcnt /= ntokens ) then
       write(*,*) 'finn_fire: infile file should have ',ntokens,' fields per line'
       write(*,*) 'finn_fire: but has ',tokcnt,' fields from input file'
       write(*,*) 'finn_fire: line # = ',n
       stop 'Data err'
     endif
     do m = 1,9
       if( m == 3 .or. m == 8 ) then
         cycle
       endif
       if( m /= 6 .and. m /= 7 ) then
         read(tokens(m),fmt=*,iostat=ios) wrk_real
         if( ios /= 0 ) then
           write(*,*) 'finn_fire: line,token # = ',n,m
           write(*,*) 'finn_fire: failed to read ',trim(raw_infile_filespec),'; error = ',ios
           stop 'Read err'
         end if
       endif
       select case( m )
         case( 1 )
           lat(n) = wrk_real
         case( 2 )
           lon(n) = wrk_real
         case( 4 )
           spix(n) = wrk_real
         case( 5 )
           tpix(n) = wrk_real
         case( 6 )
           date(n) = tokens(m)(6:7) // '/' // tokens(m)(9:10) // '/' // tokens(m)(1:4)
         case( 7 )
           read(tokens(m)(2:3),fmt=*,iostat=ios) hour
           if( ios /= 0 ) then
             write(*,*) 'read_raw_fire_file: time of day ',tokens(m)(2:3),' is invalid hour'
             stop 'Read err'
           endif
           read(tokens(m)(5:6),fmt=*,iostat=ios) min
           if( ios /= 0 ) then
             write(*,*) 'read_raw_fire_file: time of day ',tokens(m)(5:6),' is invalid minutes'
             stop 'Read err'
           endif
           tod(n) = 100*hour + min
         case( 9 )
           CONF(n) = wrk_real
       end select
     end do
   end do fire_read_loop

   close( i_unit )
   deallocate( tokens, toklen )

   end subroutine read_raw_fire_file

   subroutine write_intermediate_log

   write(unit=l_unit,fmt=*,iostat=ios) ' '
   write(unit=l_unit,fmt=*,iostat=ios) 'The Input file was: ',trim(raw_infile_filespec)
   write(unit=l_unit,fmt=*,iostat=ios) 'The emissions file was: ',trim(emis_factor_filespec)
   write(unit=l_unit,fmt=*,iostat=ios) 'The fuel load file was: ',trim(fuel_load_filespec)
   write(unit=l_unit,fmt=*,iostat=ios) ' '
   write(unit=l_unit,fmt=*,iostat=ios) 'The total number of fires input was:', nfire1
   write(unit=l_unit,fmt=*,iostat=ios) 'The total number of fires removed with confidence < 20:',nfires-nconfgt20
   write(unit=l_unit,fmt=*,iostat=ios) ' '
   write(unit=l_unit,fmt=*,iostat=ios) 'The total number of fires in the tropics was: ', ntropics
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires processed (ngoodfires):', ngoodfires
   write(unit=l_unit,fmt=*,iostat=ios) ' '
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of urban fires: ',urbnum
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires removed for overlap:', overlapct
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires skipped due to lct<= 0 or lct > 17:', lct0
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires skipped due to Global Region = Antarctica:', antarc
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires skipped due to 100% bare cover:', allbare
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires skipped due to problems with genveg:', genveg0
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires skipped due to bmass assignments:', bmass0
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires scaled to 100:', vcfcount
   write(unit=l_unit,fmt=*,iostat=ios) 'The number of fires with vcf < 50:', vcflt50
   write(unit=l_unit,fmt=*,iostat=ios) 'Total number of fires skipped:', lct0+antarc+allbare+genveg0+bmass0
   write(unit=l_unit,fmt=*,iostat=ios) ' '
   write(unit=l_unit,fmt=*,iostat=ios) 'Global Totals (Tg) of biomass burned per vegetation type'
   write(unit=l_unit,fmt=*,iostat=ios) 'GLOBAL TOTAL (Tg) biomass burned (Tg),', BMASStotal/1.e9
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Temperate Forests (Tg),', TOTTEMP*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Tropical Forests (Tg),', TOTTROP*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Boreal Forests (Tg),', TOTBOR*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Shrublands/Woody Savannah(Tg),', TOTSHRUB*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Grasslands/Savannas (Tg),', TOTGRAS*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Croplands (Tg),', TOTCROP*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) ' '
   write(unit=l_unit,fmt=*,iostat=ios) 'Global Totals (km2) of area per vegetation type'
   write(unit=l_unit,fmt=*,iostat=ios) 'TOTAL AREA BURNED (km2),', AREATOTAL
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Temperate Forests (km2),', TOTTEMParea
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Tropical Forests (km2),', TOTTROParea
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Boreal Forests (km2),', TOTBORarea
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Shrublands/Woody Savannah(km2),', TOTSHRUBarea
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Grasslands/Savannas (km2),', TOTGRASarea
   write(unit=l_unit,fmt=*,iostat=ios) 'Total Croplands (km2),', TOTCROParea
   write(unit=l_unit,fmt=*,iostat=ios) ' '
   write(unit=l_unit,fmt=*,iostat=ios) 'GLOBAL TOTALS (Tg)'
   write(unit=l_unit,fmt=*,iostat=ios) 'CO2 = ', CO2total*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'CO = ', COtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'CH4 = ', CH4total*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'NMHC = ', NMHCtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'NMOC = ', NMOCtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'H2 = ', H2total*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'NOx = ', NOXtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'NO = ', NOtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'NO2 = ', NO2total*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'SO2 = ', SO2total*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'PM2.5 = ', PM25total*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'TPM = ', TPMtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'TPC = ', TPCtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'OC = ', OCtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'BC = ', BCtotal*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'NH3 = ', NH3total*1.e-9
   write(unit=l_unit,fmt=*,iostat=ios) 'PM10 = ', PM10total*1.e-9

   end subroutine write_intermediate_log

   subroutine setup_output_and_log_files

!-----------------------------------------------------------------------------
!  ... the intermediate output file
!-----------------------------------------------------------------------------
   if( has_intermediate_output ) then
     open( unit=o_unit,file=trim(intermediate_Outfile_filespec),iostat=ios )
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to open ',trim(intermediate_Outfile_filespec),'; error = ',ios
       stop 'Open err'
     else
       write(*,*) 'finn_fire: opened output file: ',trim(intermediate_Outfile_filespec)
       write(unit=o_unit,fmt='(a)',iostat=ios) 'longi,lat,day,TIME,lct,genLC,globreg,pct_tree,pct_herb,pct_bare,' // &
              'area,bmass,CO2,CO,CH4,H2,NOx,NO,NO2,NH3,SO2,NMHC,NMOC,PM25,TPM,OC,BC,TPC,PM10,FACTOR'
       if( ios /= 0 ) then
         write(*,*) 'finn_fire: failed to write header to ',trim(intermediate_Outfile_filespec),'; error = ',ios
         stop 'Write err'
       end if
     end if
   end if
!-----------------------------------------------------------------------------
!  ... the intermediate log file
!-----------------------------------------------------------------------------
   if( has_intermediate_log ) then
     open( unit=l_unit,file=trim(intermediate_logfile_filespec),iostat=ios )
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to open ',trim(intermediate_logfile_filespec),'; error = ',ios
       stop 'Open err'
     else
       write(*,*) 'finn_fire: opened log file: ',trim(intermediate_logfile_filespec)
     end if
   end if
!-----------------------------------------------------------------------------
!  ... the speciated output file
!-----------------------------------------------------------------------------
   if( has_speciated_output ) then
     open( unit=s_unit,file=trim(speciated_Outfile_filespec),iostat=ios )
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to open ',trim(speciated_Outfile_filespec),'; error = ',ios
       stop 'Open err'
     else
       write(*,*) 'finn_fire: opened output file: ',trim(speciated_Outfile_filespec)
       select case( trim(speciation_case) )
       case( 'MOZ4' )
         s_frmt = '(I6,",",I6,",",I6,",",43(D20.10,","))'
         write(unit=s_unit,fmt='(a)',iostat=ios) 'DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,H2,NO,NO2,SO2,NH3,CH4,' // &
                                                 'NMOC,BIGALD,BIGALK,BIGENE,C10H16,C2H4,C2H5OH,C2H6,C3H6,C3H8,CH2O,' // &
                                                 'CH3CHO,CH3CN,CH3COCH3,CH3COCHO,CH3COOH,CH3OH,CRESOL,GLYALD,HCN,' // &
                                                 'HYAC,ISOP,MACR,MEK,MVK,TOLUENE,HCOOH,C2H2,OC,BC,PM10,PM25'
       case( 'SAPRC99' )
         s_frmt = '(I6,",",I6,",",I6,44(",",D20.10))'
         write(unit=s_unit,fmt='(a)',iostat=ios) 'DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,NO,NO2,SO2,NH3,CH4,VOC,ACET,' // &
                                                 'ALK1,ALK2,ALK3,ALK4,ALK5,ARO1,ARO2,BALD,CCHO,CCO_OH,ETHENE,HCHO,HCN,' // &
                                                 'HCOOH,HONO,ISOPRENE,MEK,MEOH,METHACRO,MGLY,MVK,OLE1,OLE2,PHEN,PROD2,' // &
                                                 'RCHO,TRP1,OC,BC,PM10,PM25'
       case( 'GEOSCHEM' )
         s_frmt = '(I6,",",I6,",",I6,31(",",D20.10))'
         write(unit=s_unit,fmt='(a)',iostat=ios) 'DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,NO,NO2,SO2,NH3,CH4,ACET,ALD2,' // &
                                                 'ALK4,BENZ,C2H2,C2H4,C2H6,C3H8,CH2O,GLYC,GLYX,HAC,MEK,MGLY,PRPE,' // &
                                                 'TOLU,XYLE,OC,BC,PM25'
       end select
       if( ios /= 0 ) then
         write(*,*) 'finn_fire: failed to write header to ',trim(speciated_Outfile_filespec),'; error = ',ios
         stop 'Write err'
       end if
     end if
   end if
!-----------------------------------------------------------------------------
!  ... the speciated log file
!-----------------------------------------------------------------------------
   if( has_speciated_log ) then
     open( unit=c_unit,file=trim(speciated_logfile_filespec),iostat=ios )
     if( ios /= 0 ) then
       write(*,*) 'finn_fire: failed to open ',trim(speciated_logfile_filespec),'; error = ',ios
       stop 'Open err'
     else
       write(*,*) 'finn_fire: opened log file: ',trim(speciated_logfile_filespec)
     end if
     write(*,*) 'finn_fire: Opened intermediate log file'
   end if

   end subroutine setup_output_and_log_files

   subroutine read_netcdf_file( filespec, nlons, nlats, lons, lats, data )
!-----------------------------------------------------------------------------
!  ... read netcdf input file
!-----------------------------------------------------------------------------
 
!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(out)         :: nlons, nlats
   integer(1), pointer          :: data(:,:)
   real(8), pointer             :: lons(:), lats(:)
   character(len=*), intent(in) :: filespec

!-----------------------------------------------------------------------------
!  ... read netcdf input file
!-----------------------------------------------------------------------------
   integer :: ncid
   integer :: dimid, varid
   character(len=17)  :: hdr = 'read_netcdf_file:'
   character(len=132) :: message

!---------------------------------------------------------------------
!  ... include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   write(*,*) ' '
   write(*,*) 'Reading netcdf file ' // trim(filespec)
!---------------------------------------------------------------------
!   open dataset file
!---------------------------------------------------------------------
   message = hdr // ' Failed to open ' // trim(filespec)
   call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get lon,lat dimesions
!---------------------------------------------------------------------
   message = hdr // ' Failed to get lon dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid ), message )
   message = hdr // ' Failed to get lon dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlons ), message )
   message = hdr // ' Failed to get lat dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid ), message )
   message = hdr // ' Failed to get lat dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlats ), message )
   write(*,*) hdr // '  nlons, nlats = ',nlons,nlats

!---------------------------------------------------------------------
!   allocate lons
!---------------------------------------------------------------------
   if( associated( lons ) ) then
     deallocate( lons )
   endif
   allocate( lons(nlons),stat=astat )
   if( astat /= 0 ) then
     write(*,*) hdr // ' Failed to allocate lons; error = ',astat
     stop 'Alloc err'
   endif
!---------------------------------------------------------------------
!   read longitudes
!---------------------------------------------------------------------
   message = hdr // ' Failed to get lon variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), message )
   message = hdr // ' Failed to read lon variable'
   call handle_ncerr( nf_get_var_double( ncid, varid, lons ), message )

!---------------------------------------------------------------------
!   allocate lats
!---------------------------------------------------------------------
   if( associated( lats ) ) then
     deallocate( lats )
   endif
   allocate( lats(nlats),stat=astat )
   if( astat /= 0 ) then
     write(*,*) hdr // ' Failed to allocate lats; error = ',astat
     stop 'Alloc err'
   endif
!---------------------------------------------------------------------
!   read latitudes
!---------------------------------------------------------------------
   message = hdr // ' Failed to get lat variable id'
   call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), message )
   message = hdr // ' Failed to read lat variable'
   call handle_ncerr( nf_get_var_double( ncid, varid, lats ), message )
!---------------------------------------------------------------------
!   allocate data
!---------------------------------------------------------------------
   if( associated( data ) ) then
     deallocate( data )
   endif
   allocate( data(nlons,nlats),stat=astat )
   if( astat /= 0 ) then
     write(*,*) hdr // ' Failed to allocate data; error = ',astat
     stop 'Alloc err'
   endif
!---------------------------------------------------------------------
!   read data
!---------------------------------------------------------------------
   message = hdr // ' Failed to get ' // trim(var_name) // ' data variable id'
   call handle_ncerr( nf_inq_varid( ncid, trim(var_name), varid ), message )
   message = hdr // ' Failed to read ' // trim(var_name) // ' data variable'
   call handle_ncerr( nf_get_var_int1( ncid, varid, data ), message )
!---------------------------------------------------------------------
!   close file
!---------------------------------------------------------------------
   message = hdr // ' Failed to close ' // trim(filespec)
   call handle_ncerr( nf_close( ncid ), message )       

   end subroutine read_netcdf_file

   subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!	... netcdf error handling routine
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ret
   character(len=*), intent(in) :: mes

!---------------------------------------------------------------------
!  ... include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   if( ret /= nf_noerr ) then
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

   subroutine intrvl( t, ndx, nt, x, nx )
!---------------------------------------------------------------------
!  ... find interval enclosing point
!---------------------------------------------------------------------

   integer, intent(in)  :: nt
   real(8), intent(in)  :: t(nt)
   integer, intent(in)  :: nx
   real(8), intent(in)  :: x(nx)
   integer, intent(out) :: ndx(nt)

   integer :: k, kl, kh, l, nxm1
   real(8) :: tt

   l = 1
   nxm1 = nx - 1

target_loop : &
   do k = 1,nt
     tt = t(k)
     if( tt < x(l) ) then
       if( tt <= x(2) ) then
         l = 1
         ndx(k) = 1
         cycle target_loop
       else
         kl = 2; kh = l
       endif
     elseif( tt <= x(l+1) ) then
       ndx(k) = l
       cycle target_loop
     elseif( tt >= x(nxm1) ) then
       l = nxm1
       ndx(k) = l
       cycle target_loop
     else
       kl = l+1; kh = nxm1
     endif

     do
       l = (kl + kh)/2
       if( tt < x(l) ) then
         kh = l
       elseif( tt > x(l+1) ) then
         kl = l + 1
       else
         ndx(k) = l
         cycle target_loop
       endif
     end do
   end do target_loop
   
   end subroutine intrvl

   subroutine get_lct_fractions( filespec, real_fraction, int_fraction, do_average )
!---------------------------------------------------------------------
!  ... read and set the land type fraction
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  ... dummy arguments
!---------------------------------------------------------------------
   integer, pointer, optional :: int_fraction(:)
   real, pointer, optional    :: real_fraction(:)
   logical                    :: do_average
   character(len=*), intent(in) :: filespec

!---------------------------------------------------------------------
!  ... local variables
!---------------------------------------------------------------------
   integer                 :: nlons, nlats
   integer                 :: niu, nil
   integer                 :: nju, njl
   integer                 :: npartial
   integer(1), pointer     :: frac(:,:)
   integer(1), allocatable :: wrk_int1(:)
   integer(2)              :: div, average
   integer, allocatable    :: lon_ndx(:), lat_ndx(:)
   real(8), pointer        :: lons(:), lats(:)
   real(8), allocatable    :: elons(:), elats(:)

   call read_netcdf_file( filespec, nlons, nlats, lons, lats, frac )

   allocate( lon_ndx(nfires),lat_ndx(nfires), &
             elons(nlons+1),elats(nlats+1),stat=astat )
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate lon_ndx ... elats; error = ',astat
     stop 'Alloc err'
   end if

   elons(2:nlons) = .5_8*(lons(1:nlons-1) + lons(2:nlons))
   elons(1)       = lons(1) - (elons(2) - lons(1))
   elons(nlons+1) = lons(nlons) + (lons(nlons) - elons(nlons))

   call intrvl( real(lon(:),kind=8), lon_ndx, nfires, elons, nlons+1 )

!---------------------------------------------------------------------
!  ... if lats are monotonically decreasing then reorder
!---------------------------------------------------------------------
   if( lats(1) > lats(2) ) then
     allocate( wrk_int1(nlons),stat=astat )
     if( astat /= 0 ) then
       write(*,*) 'finn_fire: failed to allocate wrk_int1; error = ',astat
       stop 'Alloc err'
     end if
     do n = 1,nlats/2
       m        = nlats - n + 1
       wrk_real = lats(n)
       wrk_int1(:) = frac(:,n)
       lats(n)     = lats(m)
       frac(:,n)   = frac(:,m)
       lats(m)     = wrk_real
       frac(:,m)   = wrk_int1(:)
     end do
     deallocate( wrk_int1 )
   endif

   elats(2:nlats) = .5_8*(lats(1:nlats-1) + lats(2:nlats))
   elats(1)       = lats(1) - (elats(2) - lats(1))
   elats(nlats+1) = lats(nlats) + (lats(nlats) - elats(nlats))

   call intrvl( real(lat(:),kind=8), lat_ndx, nfires, elats, nlats+1 )

   write(*,*) 'min,max lon indices = ',minval(lon_ndx(:)),maxval(lon_ndx(:))
   write(*,*) 'min,max lat indices = ',minval(lat_ndx(:)),maxval(lat_ndx(:))

!---------------------------------------------------------------------
!  ... replace missing data with zero
!---------------------------------------------------------------------
   where( frac(:,:) == -1_1 )
     frac(:,:) = 0_1
   endwhere

   if( present( real_fraction ) ) then
     allocate( real_fraction(nfires),stat=astat )
   elseif( present( int_fraction ) ) then
     allocate( int_fraction(nfires),stat=astat )
   end if
   if( astat /= 0 ) then
     write(*,*) 'finn_fire: failed to allocate fraction; error = ',astat
     stop 'Alloc err'
   end if

   niu = 0 ; nil = 0 ; njl = 0 ; nju = 0 ; npartial = 0
   do n = 1,nfires
     i = lon_ndx(n) ; j = lat_ndx(n)
     if( do_average ) then
       if( lon(n) >= lons(i) ) then
         il = max( 1,i )
         niu = niu + 1 
       else
         il = max( 1,i-1 )
         nil = nil + 1 
       endif
       iu = min( nlons,il+1 )
       if( lat(n) >= lats(j) ) then
         jl = max( 1,j )
         nju = nju + 1 
       else
         jl = max( 1,j-1 )
         njl = njl + 1 
       endif
       ju = min( nlats,jl+1 )
       div = iu - il + ju - jl + 2
       if( div < 4_2 ) then
         div = div - 1_2
         npartial = npartial + 1
       endif
!      il = max( 1,i-1 ) ; iu = min( nlons,i+1 )
!      jl = max( 1,j-1 ) ; ju = min( nlats,j+1 )
!      average = sum( int(frac(il:iu,j),kind=2) ) + sum( int(frac(i,jl:ju),kind=2) )
!      average = average - int(frac(i,j),kind=2)
       average = sum( int(frac(il:iu,jl),kind=2) ) + sum( int(frac(il:iu,ju),kind=2) )
!      average = average - int(frac(i,j),kind=2)
       if( present( real_fraction ) ) then
!        real_fraction(n) = real(average/5_2)
         real_fraction(n) = real(average/div)
       elseif( present( int_fraction ) ) then
!        int_fraction(n) = int(average/5_2)
         int_fraction(n) = int(average/div)
       endif
       if( n == 10421 ) then
         write(*,*) '###################################################'
         write(*,'(''get_lct_frac: i,j         = '',2i6)') i,j
         write(*,'(''get_lct_frac: il,iu,jl,ju = '',4i6)') il,iu,jl,ju
         write(*,'(''get_lct_frac: fire lon,lat = '',1p2g15.7)') lon(n),lat(n)
         write(*,'(''get_lct_frac: grid lons    = '',1p2g15.7)') lons(il),lons(iu)
         write(*,'(''get_lct_frac: grid lats    = '',1p2g15.7)') lats(jl),lats(ju)
         write(*,*) 'get_lct_frac: frac'
         write(*,*) frac(il:iu,jl:ju)
         write(*,*) 'get_lct_frac: div,average = ',div,average
         write(*,*) '###################################################'
       endif
     else
       if( present( real_fraction ) ) then
         real_fraction(n) = real(frac(i,j))
       elseif( present( int_fraction ) ) then
         int_fraction(n) = int(frac(i,j))
       endif
     endif
   end do

   if( do_average ) then
     write(*,'(''get_lct_frac: nil,niu,njl,nju = '',4i6)') nil,niu,njl,nju
     write(*,'(''get_lct_frac: npartial        = '',i6)') npartial
   endif

   deallocate( lon_ndx, lat_ndx, lons, lats, elons, elats, frac )

   end subroutine get_lct_fractions

   subroutine reverse_array( n, arr )
!-----------------------------------------------------------------------------
!  ... reverse order in input array
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!  ... dummy arguments
!-----------------------------------------------------------------------------
   integer, intent(in)    :: n
   integer, intent(inout) :: arr(n)

!-----------------------------------------------------------------------------
!  ... local variables
!-----------------------------------------------------------------------------
   integer :: i, wrk
   integer :: nn, nh

   nh = n/2

   do i = 1,nh
     nn      = n - i + 1
     wrk     = arr(i)
     arr(i)  = arr(nn)
     arr(nn) = wrk
   end do

   end subroutine reverse_array

   end program finn_fire_tst
