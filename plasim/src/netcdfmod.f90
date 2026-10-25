!     =================
!     MODULE NETCDFMOD
!     =================

      module netcdfmod
      use pumamod
      use netcdf
      implicit none

      ! NetCDF dataset and dimension IDs
      integer :: ncid = -1
      integer :: dim_lon, dim_lat, dim_lev, dim_time
      integer :: vid_lon, vid_lat, vid_lev, vid_time

      ! 3D variable IDs (dimensions: lon, lat, lev, time)
      integer :: vid_ua, vid_va, vid_ta, vid_hus, vid_clw, vid_cl

      ! 2D variable IDs (dimensions: lon, lat, time)
      integer :: vid_ps, vid_ts, vid_tas, vid_tasmax, vid_tasmin
      integer :: vid_mrso, vid_snd, vid_prl, vid_prc, vid_pr, vid_prsn
      integer :: vid_evap, vid_mrro, vid_hfss, vid_hfls
      integer :: vid_rss, vid_rls, vid_rst, vid_rlut
      integer :: vid_rsut, vid_ssru, vid_stru
      integer :: vid_tauu, vid_tauv, vid_clt, vid_u3
      integer :: vid_tso, vid_tso2, vid_tso3, vid_tso4, vid_tsod
      integer :: vid_sic, vid_sit, vid_mld, vid_snm, vid_sndc
      integer :: vid_prw, vid_lsm, vid_z0, vid_as, vid_vegf
      integer :: vid_glac, vid_dwmax, vid_sg

      ! Ocean dataset and dimension IDs
      integer :: ncid_ocean = -1
      integer :: dim_lon_oce, dim_lat_oce, dim_time_oce
      integer :: vid_lon_oce, vid_lat_oce, vid_time_oce
      integer :: vid_oce_heata, vid_oce_ifluxa, vid_oce_fssta, vid_oce_dssta
      integer :: vid_oce_qhda, vid_oce_fldoa, vid_oce_icec, vid_oce_sst
      integer :: vid_oce_ls, vid_oce_clsst
      integer :: ocean_nc_time_step = 0

      ! Sea ice dataset and dimension IDs
      integer :: ncid_ice = -1
      integer :: dim_lon_ice, dim_lat_ice, dim_time_ice
      integer :: vid_lon_ice, vid_lat_ice, vid_time_ice
      integer :: vid_ice_heata, vid_ice_ofluxa, vid_ice_tsfluxa, vid_ice_smelta
      integer :: vid_ice_imelta, vid_ice_cfluxa, vid_ice_fluxca, vid_ice_qmelta
      integer :: vid_ice_xflxicea, vid_ice_icec, vid_ice_iced, vid_ice_scflxa
      integer :: vid_ice_cfluxra, vid_ice_cfluxna, vid_ice_ts, vid_ice_zsnow
      integer :: vid_ice_sst, vid_ice_ls, vid_ice_clicec2, vid_ice_cliced2
      integer :: vid_ice_icecc, vid_ice_cpmea, vid_ice_croffa, vid_ice_stoia
      integer :: ice_nc_time_step = 0

      ! Tracking counters
      integer :: nc_time_step = 0
      integer :: nc_acc_count = 0

      ! Static gridpoint surface geopotential (orography)
      real, allocatable :: oro_gp(:)

      ! Accumulation arrays for time-averaging
      real, allocatable :: acc_u(:,:), acc_v(:,:), acc_t(:,:)
      real, allocatable :: acc_q(:,:), acc_clw(:,:), acc_cl(:,:)
      real, allocatable :: acc_ps(:), acc_ts(:), acc_tas(:)
      real, allocatable :: acc_tasmax(:), acc_tasmin(:)
      real, allocatable :: acc_mrso(:), acc_snd(:), acc_prl(:), acc_prc(:)
      real, allocatable :: acc_prsn(:), acc_evap(:), acc_mrro(:)
      real, allocatable :: acc_hfss(:), acc_hfls(:)
      real, allocatable :: acc_rss(:), acc_rls(:), acc_rst(:), acc_rlut(:)
      real, allocatable :: acc_rsut(:), acc_ssru(:), acc_stru(:)
      real, allocatable :: acc_tauu(:), acc_tauv(:), acc_clt(:), acc_u3(:)
      real, allocatable :: acc_tso(:), acc_tso2(:), acc_tso3(:), acc_tso4(:), acc_tsod(:)
      real, allocatable :: acc_sic(:), acc_sit(:), acc_mld(:), acc_snm(:), acc_sndc(:)
      real, allocatable :: acc_prw(:)

      contains

      ! ------------------------------------------------------------------
      !  nc_check: error handler
      ! ------------------------------------------------------------------
      subroutine nc_check(status, context)
         integer,          intent(in) :: status
         character(len=*), intent(in) :: context
         if (status /= NF90_NOERR) then
            write(*,*) 'NetCDF error in ', trim(context), ': ', trim(nf90_strerror(status))
            call mpabort('NetCDF error occurred')
         endif
      end subroutine nc_check

      ! ------------------------------------------------------------------
      !  def_var_2d_file: define 2D float variable with deflate and chunking in target file
      ! ------------------------------------------------------------------
      subroutine def_var_2d_file(target_ncid, d_lon, d_lat, d_time, name, long_name, units, code, vid)
         integer,          intent(in)  :: target_ncid, d_lon, d_lat, d_time
         character(len=*), intent(in)  :: name, long_name, units
         integer,          intent(in)  :: code
         integer,          intent(out) :: vid
         integer                       :: dims(3), chunk(3), ret

         dims = (/ d_lon, d_lat, d_time /)
         chunk = (/ NLON, NLAT, 1 /)

         ret = nf90_def_var(target_ncid, trim(name), NF90_FLOAT, dims, vid)
         call nc_check(ret, 'def_var ' // trim(name))

         if (ndeflate > 0) then
            ret = nf90_def_var_deflate(target_ncid, vid, shuffle=1, deflate=1, deflate_level=ndeflate)
            call nc_check(ret, 'deflate ' // trim(name))
         endif

         ret = nf90_def_var_chunking(target_ncid, vid, NF90_CHUNKED, chunk)
         call nc_check(ret, 'chunking ' // trim(name))

         ret = nf90_put_att(target_ncid, vid, 'long_name', trim(long_name))
         call nc_check(ret, 'put_att long_name ' // trim(name))

         ret = nf90_put_att(target_ncid, vid, 'units', trim(units))
         call nc_check(ret, 'put_att units ' // trim(name))

         ret = nf90_put_att(target_ncid, vid, 'code', code)
         call nc_check(ret, 'put_att code ' // trim(name))
      end subroutine def_var_2d_file

      ! ------------------------------------------------------------------
      !  def_var_2d: define 2D float variable with deflate and chunking in primary dataset
      ! ------------------------------------------------------------------
      subroutine def_var_2d(name, long_name, units, code, vid)
         character(len=*), intent(in)  :: name, long_name, units
         integer,          intent(in)  :: code
         integer,          intent(out) :: vid
         call def_var_2d_file(ncid, dim_lon, dim_lat, dim_time, name, long_name, units, code, vid)
      end subroutine def_var_2d

      ! ------------------------------------------------------------------
      !  def_var_3d: define 3D float variable with deflate and chunking
      ! ------------------------------------------------------------------
      subroutine def_var_3d(name, long_name, units, code, vid)
         character(len=*), intent(in)  :: name, long_name, units
         integer,          intent(in)  :: code
         integer,          intent(out) :: vid
         integer                       :: dims(4), chunk(4), ret

         dims = (/ dim_lon, dim_lat, dim_lev, dim_time /)
         chunk = (/ NLON, NLAT, 1, 1 /)

         ret = nf90_def_var(ncid, trim(name), NF90_FLOAT, dims, vid)
         call nc_check(ret, 'def_var ' // trim(name))

         if (ndeflate > 0) then
            ret = nf90_def_var_deflate(ncid, vid, shuffle=1, deflate=1, deflate_level=ndeflate)
            call nc_check(ret, 'deflate ' // trim(name))
         endif

         ret = nf90_def_var_chunking(ncid, vid, NF90_CHUNKED, chunk)
         call nc_check(ret, 'chunking ' // trim(name))

         ret = nf90_put_att(ncid, vid, 'long_name', trim(long_name))
         call nc_check(ret, 'put_att long_name ' // trim(name))

         ret = nf90_put_att(ncid, vid, 'units', trim(units))
         call nc_check(ret, 'put_att units ' // trim(name))

         ret = nf90_put_att(ncid, vid, 'code', code)
         call nc_check(ret, 'put_att code ' // trim(name))
      end subroutine def_var_3d

      ! ------------------------------------------------------------------
      !  reset_accumulators: zeroes accumulation fields
      ! ------------------------------------------------------------------
      subroutine reset_accumulators
         acc_u(:,:)   = 0.0
         acc_v(:,:)   = 0.0
         acc_t(:,:)   = 0.0
         acc_q(:,:)   = 0.0
         acc_clw(:,:) = 0.0
         acc_cl(:,:)  = 0.0

         acc_ps(:)    = 0.0
         acc_ts(:)    = 0.0
         acc_tas(:)   = 0.0
         acc_tasmax(:)= -1.0e10
         acc_tasmin(:)=  1.0e10
         acc_mrso(:)  = 0.0
         acc_snd(:)   = 0.0
         acc_prl(:)   = 0.0
         acc_prc(:)   = 0.0
         acc_prsn(:)  = 0.0
         acc_evap(:)  = 0.0
         acc_mrro(:)  = 0.0
         acc_hfss(:)  = 0.0
         acc_hfls(:)  = 0.0
         acc_rss(:)   = 0.0
         acc_rls(:)   = 0.0
         acc_rst(:)   = 0.0
         acc_rlut(:)  = 0.0
         acc_rsut(:)  = 0.0
         acc_ssru(:)  = 0.0
         acc_stru(:)  = 0.0
         acc_tauu(:)  = 0.0
         acc_tauv(:)  = 0.0
         acc_clt(:)   = 0.0
         acc_u3(:)    = 0.0
         acc_tso(:)   = 0.0
         acc_tso2(:)  = 0.0
         acc_tso3(:)  = 0.0
         acc_tso4(:)  = 0.0
         acc_tsod(:)  = 0.0
         acc_sic(:)   = 0.0
         acc_sit(:)   = 0.0
         acc_mld(:)   = 0.0
         acc_snm(:)   = 0.0
         acc_sndc(:)  = 0.0
         acc_prw(:)   = 0.0

         nc_acc_count = 0
      end subroutine reset_accumulators

      ! ------------------------------------------------------------------
      !  inq_var_ids: inquire existing variable IDs when reopening
      ! ------------------------------------------------------------------
      subroutine inq_var_ids
         integer :: ret
         ret = nf90_inq_varid(ncid, 'time',   vid_time)
         ret = nf90_inq_varid(ncid, 'ua',     vid_ua)
         ret = nf90_inq_varid(ncid, 'va',     vid_va)
         ret = nf90_inq_varid(ncid, 'ta',     vid_ta)
         ret = nf90_inq_varid(ncid, 'hus',    vid_hus)
         ret = nf90_inq_varid(ncid, 'clw',    vid_clw)
         ret = nf90_inq_varid(ncid, 'cl',     vid_cl)
         ret = nf90_inq_varid(ncid, 'ps',     vid_ps)
         ret = nf90_inq_varid(ncid, 'ts',     vid_ts)
         ret = nf90_inq_varid(ncid, 'tas',    vid_tas)
         ret = nf90_inq_varid(ncid, 'tasmax', vid_tasmax)
         ret = nf90_inq_varid(ncid, 'tasmin', vid_tasmin)
         ret = nf90_inq_varid(ncid, 'mrso',   vid_mrso)
         ret = nf90_inq_varid(ncid, 'snd',    vid_snd)
         ret = nf90_inq_varid(ncid, 'prl',    vid_prl)
         ret = nf90_inq_varid(ncid, 'prc',    vid_prc)
         ret = nf90_inq_varid(ncid, 'pr',     vid_pr)
         ret = nf90_inq_varid(ncid, 'prsn',   vid_prsn)
         ret = nf90_inq_varid(ncid, 'evap',   vid_evap)
         ret = nf90_inq_varid(ncid, 'mrro',   vid_mrro)
         ret = nf90_inq_varid(ncid, 'hfss',   vid_hfss)
         ret = nf90_inq_varid(ncid, 'hfls',   vid_hfls)
         ret = nf90_inq_varid(ncid, 'rss',    vid_rss)
         ret = nf90_inq_varid(ncid, 'rls',    vid_rls)
         ret = nf90_inq_varid(ncid, 'rst',    vid_rst)
         ret = nf90_inq_varid(ncid, 'rlut',   vid_rlut)
         ret = nf90_inq_varid(ncid, 'rsut',   vid_rsut)
         ret = nf90_inq_varid(ncid, 'ssru',   vid_ssru)
         ret = nf90_inq_varid(ncid, 'stru',   vid_stru)
         ret = nf90_inq_varid(ncid, 'tauu',   vid_tauu)
         ret = nf90_inq_varid(ncid, 'tauv',   vid_tauv)
         ret = nf90_inq_varid(ncid, 'clt',    vid_clt)
         ret = nf90_inq_varid(ncid, 'u3',     vid_u3)
         ret = nf90_inq_varid(ncid, 'tso',    vid_tso)
         ret = nf90_inq_varid(ncid, 'tso2',   vid_tso2)
         ret = nf90_inq_varid(ncid, 'tso3',   vid_tso3)
         ret = nf90_inq_varid(ncid, 'tso4',   vid_tso4)
         ret = nf90_inq_varid(ncid, 'tsod',   vid_tsod)
         ret = nf90_inq_varid(ncid, 'sic',    vid_sic)
         ret = nf90_inq_varid(ncid, 'sit',    vid_sit)
         ret = nf90_inq_varid(ncid, 'mld',    vid_mld)
         ret = nf90_inq_varid(ncid, 'snm',    vid_snm)
         ret = nf90_inq_varid(ncid, 'sndc',   vid_sndc)
         ret = nf90_inq_varid(ncid, 'prw',    vid_prw)
         ret = nf90_inq_varid(ncid, 'lsm',    vid_lsm)
         ret = nf90_inq_varid(ncid, 'z0',     vid_z0)
         ret = nf90_inq_varid(ncid, 'as',     vid_as)
         ret = nf90_inq_varid(ncid, 'vegf',   vid_vegf)
         ret = nf90_inq_varid(ncid, 'glac',   vid_glac)
         ret = nf90_inq_varid(ncid, 'dwmax',  vid_dwmax)
         ret = nf90_inq_varid(ncid, 'sg',     vid_sg)
      end subroutine inq_var_ids

      ! ------------------------------------------------------------------
      !  netcdf_init: create or open NetCDF file, allocate accumulators
      ! ------------------------------------------------------------------
      subroutine netcdf_init
         character(len=256) :: fname, time_units
         real(kind=4)       :: coord_lon(NLON), coord_lat(NLAT), coord_lev(NLEV)
         real               :: zfc(NHOR)
         integer            :: jlon, jlat, jlev, ret, did_time, existing_len

         if (nnetcdf == 0 .or. noutput == 0) return

         ! Allocate local accumulation fields
         if (.not. allocated(acc_u)) then
            allocate(acc_u(NHOR, NLEV), acc_v(NHOR, NLEV), acc_t(NHOR, NLEV))
            allocate(acc_q(NHOR, NLEV), acc_clw(NHOR, NLEV), acc_cl(NHOR, NLEV))
            allocate(acc_ps(NHOR), acc_ts(NHOR), acc_tas(NHOR))
            allocate(acc_tasmax(NHOR), acc_tasmin(NHOR))
            allocate(acc_mrso(NHOR), acc_snd(NHOR), acc_prl(NHOR), acc_prc(NHOR))
            allocate(acc_prsn(NHOR), acc_evap(NHOR), acc_mrro(NHOR))
            allocate(acc_hfss(NHOR), acc_hfls(NHOR))
            allocate(acc_rss(NHOR), acc_rls(NHOR), acc_rst(NHOR), acc_rlut(NHOR))
            allocate(acc_rsut(NHOR), acc_ssru(NHOR), acc_stru(NHOR))
            allocate(acc_tauu(NHOR), acc_tauv(NHOR), acc_clt(NHOR), acc_u3(NHOR))
            allocate(acc_tso(NHOR), acc_tso2(NHOR), acc_tso3(NHOR), acc_tso4(NHOR), acc_tsod(NHOR))
            allocate(acc_sic(NHOR), acc_sit(NHOR), acc_mld(NHOR), acc_snm(NHOR), acc_sndc(NHOR))
            allocate(acc_prw(NHOR))
            allocate(oro_gp(NHOR))
         endif

         call reset_accumulators

         ! Compute gridpoint orography from spectral so
         zfc(:) = 0.0
         call sp2fc(so, zfc)
         call fc2gp(zfc, NLON, NLPP)
         oro_gp(:) = zfc(:) * (cv * cv)

         if (mypid /= NROOT) return

         ! Determine output filename
         fname = trim(plasim_output)
         if (index(fname, '.nc') == 0) then
            fname = trim(fname) // '.nc'
         endif

         ! On restart, try to append to existing file
         if (nrestart > 0) then
            ret = nf90_open(trim(fname), NF90_WRITE, ncid)
            if (ret == NF90_NOERR) then
               ret = nf90_inq_dimid(ncid, 'time', did_time)
               call nc_check(ret, 'inq_dimid time')
               ret = nf90_inquire_dimension(ncid, did_time, len=existing_len)
               call nc_check(ret, 'inquire_dim time')
               nc_time_step = existing_len
               call inq_var_ids
               return
            endif
         endif

         ! Create new NetCDF-4 file
         ret = nf90_create(trim(fname), IOR(NF90_CLOBBER, NF90_NETCDF4), ncid)
         call nc_check(ret, 'create ' // trim(fname))

         ! Global attributes
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'title', 'PlaSim-LSG TNG Direct NetCDF Output')
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.6')
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'source', 'PlaSim-LSG TNG')
         if (noutput_avg == 2) then
            ret = nf90_put_att(ncid, NF90_GLOBAL, 'frequency', 'mon')
         else
            ret = nf90_put_att(ncid, NF90_GLOBAL, 'frequency', 'day')
         endif

         ! Define dimensions
         ret = nf90_def_dim(ncid, 'lon', NLON, dim_lon)
         call nc_check(ret, 'def_dim lon')
         ret = nf90_def_dim(ncid, 'lat', NLAT, dim_lat)
         call nc_check(ret, 'def_dim lat')
         ret = nf90_def_dim(ncid, 'lev', NLEV, dim_lev)
         call nc_check(ret, 'def_dim lev')
         ret = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dim_time)
         call nc_check(ret, 'def_dim time')

         ! Define coordinate variables
         ret = nf90_def_var(ncid, 'lon', NF90_FLOAT, (/ dim_lon /), vid_lon)
         call nc_check(ret, 'def_var lon')
         ret = nf90_put_att(ncid, vid_lon, 'standard_name', 'longitude')
         ret = nf90_put_att(ncid, vid_lon, 'long_name', 'longitude')
         ret = nf90_put_att(ncid, vid_lon, 'units', 'degrees_east')
         ret = nf90_put_att(ncid, vid_lon, 'axis', 'X')

         ret = nf90_def_var(ncid, 'lat', NF90_FLOAT, (/ dim_lat /), vid_lat)
         call nc_check(ret, 'def_var lat')
         ret = nf90_put_att(ncid, vid_lat, 'standard_name', 'latitude')
         ret = nf90_put_att(ncid, vid_lat, 'long_name', 'latitude')
         ret = nf90_put_att(ncid, vid_lat, 'units', 'degrees_north')
         ret = nf90_put_att(ncid, vid_lat, 'axis', 'Y')

         ret = nf90_def_var(ncid, 'lev', NF90_FLOAT, (/ dim_lev /), vid_lev)
         call nc_check(ret, 'def_var lev')
         ret = nf90_put_att(ncid, vid_lev, 'standard_name', 'atmosphere_sigma_coordinate')
         ret = nf90_put_att(ncid, vid_lev, 'long_name', 'sigma at layer midpoints')
         ret = nf90_put_att(ncid, vid_lev, 'units', '1')
         ret = nf90_put_att(ncid, vid_lev, 'positive', 'down')
         ret = nf90_put_att(ncid, vid_lev, 'axis', 'Z')

         ret = nf90_def_var(ncid, 'time', NF90_DOUBLE, (/ dim_time /), vid_time)
         call nc_check(ret, 'def_var time')
         ret = nf90_put_att(ncid, vid_time, 'standard_name', 'time')
         ret = nf90_put_att(ncid, vid_time, 'long_name', 'time')
         write(time_units, '("days since ",I4.4,"-",I2.2,"-01 00:00:00")') n_start_year, n_start_month
         ret = nf90_put_att(ncid, vid_time, 'units', trim(time_units))
         if (n_days_per_year == 360) then
            ret = nf90_put_att(ncid, vid_time, 'calendar', '360_day')
         else if (n_days_per_year == 365) then
            ret = nf90_put_att(ncid, vid_time, 'calendar', '365_day')
         else
            ret = nf90_put_att(ncid, vid_time, 'calendar', 'proleptic_gregorian')
         endif
         ret = nf90_put_att(ncid, vid_time, 'axis', 'T')

         ! Define 3D variables
         call def_var_3d('ua',     'Zonal Wind',                      'm/s',      131, vid_ua)
         call def_var_3d('va',     'Meridional Wind',                 'm/s',      132, vid_va)
         call def_var_3d('ta',     'Temperature',                     'K',        130, vid_ta)
         call def_var_3d('hus',    'Specific Humidity',               'kg/kg',    133, vid_hus)
         call def_var_3d('clw',    'Liquid Water Content',            'kg/kg',    161, vid_clw)
         call def_var_3d('cl',     'Cloud Cover',                     'fraction', 162, vid_cl)

         ! Define 2D variables
         call def_var_2d('ps',     'Surface Pressure',                'Pa',       134, vid_ps)
         call def_var_2d('ts',     'Surface Temperature',             'K',        139, vid_ts)
         call def_var_2d('tas',    'Temperature at 2m',               'K',        167, vid_tas)
         call def_var_2d('tasmax', 'Maximum 2m Temperature',          'K',        201, vid_tasmax)
         call def_var_2d('tasmin', 'Minimum 2m Temperature',          'K',        202, vid_tasmin)
         call def_var_2d('mrso',   'Soil Wetness',                    'm',        140, vid_mrso)
         call def_var_2d('snd',    'Snow Depth',                      'm',        141, vid_snd)
         call def_var_2d('prl',    'Large Scale Precipitation',       'm/s',      142, vid_prl)
         call def_var_2d('prc',    'Convective Precipitation',        'm/s',      143, vid_prc)
         call def_var_2d('pr',     'Total Precipitation',             'm/s',      260, vid_pr)
         call def_var_2d('prsn',   'Snow Fall',                       'm/s',      144, vid_prsn)
         call def_var_2d('evap',   'Evaporation',                     'm/s',      182, vid_evap)
         call def_var_2d('mrro',   'Surface Runoff',                  'm/s',      160, vid_mrro)
         call def_var_2d('hfss',   'Surface Sensible Heat Flux',      'W/m2',     146, vid_hfss)
         call def_var_2d('hfls',   'Surface Latent Heat Flux',        'W/m2',     147, vid_hfls)
         call def_var_2d('rss',    'Surface Solar Radiation',         'W/m2',     176, vid_rss)
         call def_var_2d('rls',    'Surface Thermal Radiation',       'W/m2',     177, vid_rls)
         call def_var_2d('rst',    'Top Solar Radiation',             'W/m2',     178, vid_rst)
         call def_var_2d('rlut',   'Top Thermal Radiation',           'W/m2',     179, vid_rlut)
         call def_var_2d('rsut',   'Top Solar Radiation Upward',      'W/m2',     203, vid_rsut)
         call def_var_2d('ssru',   'Surface Solar Radiation Upward',  'W/m2',     204, vid_ssru)
         call def_var_2d('stru',   'Surface Thermal Radiation Upward','W/m2',     205, vid_stru)
         call def_var_2d('tauu',   'U-Stress',                        'Pa',       180, vid_tauu)
         call def_var_2d('tauv',   'V-Stress',                        'Pa',       181, vid_tauv)
         call def_var_2d('clt',    'Total Cloud Cover',               'fraction', 164, vid_clt)
         call def_var_2d('u3',     'ustar **3',                       'm3/s3',    159, vid_u3)
         call def_var_2d('tso',    'Soil Temperature',                'K',        183, vid_tso)
         call def_var_2d('tso2',   'Soil Temperature Level 2',        'K',        207, vid_tso2)
         call def_var_2d('tso3',   'Soil Temperature Level 3',        'K',        208, vid_tso3)
         call def_var_2d('tso4',   'Soil Temperature Level 4',        'K',        209, vid_tso4)
         call def_var_2d('tsod',   'Deep Soil Temperature',           'K',        170, vid_tsod)
         call def_var_2d('sic',    'Sea Ice Cover',                   'fraction', 210, vid_sic)
         call def_var_2d('sit',    'Sea Ice Thickness',               'm',        211, vid_sit)
         call def_var_2d('mld',    'Mixed Layer Depth',               'm',        110, vid_mld)
         call def_var_2d('snm',    'Snow Melt',                       'm/s',      218, vid_snm)
         call def_var_2d('sndc',   'Snow Depth Change',               'm/s',      221, vid_sndc)
         call def_var_2d('prw',    'Vert. Integrated Spec. Hum.',     'kg/m2',    230, vid_prw)
         call def_var_2d('lsm',    'Land Sea Mask',                   'fraction', 172, vid_lsm)
         call def_var_2d('z0',     'Surface Roughness',               'm',        173, vid_z0)
         call def_var_2d('as',     'Surface Albedo',                  'fraction', 175, vid_as)
         call def_var_2d('vegf',   'Forest Cover',                    'fraction', 212, vid_vegf)
         call def_var_2d('glac',   'Glacier Cover',                   'fraction', 232, vid_glac)
         call def_var_2d('dwmax',  'Field capacity',                  'm',        229, vid_dwmax)
         call def_var_2d('sg',     'Surface Geopotential Orography',  'm2/s2',    129, vid_sg)

         ! End definition mode
         ret = nf90_enddef(ncid)
         call nc_check(ret, 'enddef')

         ! Write static coordinates
         do jlon = 1, NLON
            coord_lon(jlon) = real((jlon - 1) * 360.0 / NLON, kind=4)
         enddo
         ret = nf90_put_var(ncid, vid_lon, coord_lon)
         call nc_check(ret, 'put_var lon')

         do jlat = 1, NLAT
            coord_lat(jlat) = real(deglat(jlat), kind=4)
         enddo
         ret = nf90_put_var(ncid, vid_lat, coord_lat)
         call nc_check(ret, 'put_var lat')

         do jlev = 1, NLEV
            coord_lev(jlev) = real(sigma(jlev), kind=4)
         enddo
         ret = nf90_put_var(ncid, vid_lev, coord_lev)
         call nc_check(ret, 'put_var lev')

         ret = nf90_sync(ncid)
         call nc_check(ret, 'sync')
      end subroutine netcdf_init

      ! ------------------------------------------------------------------
      !  write_field_2d_file: gather and write 2D variable at specified record in target file
      ! ------------------------------------------------------------------
      subroutine write_field_2d_file(target_ncid, rec_step, vid, loc_field)
         integer, intent(in) :: target_ncid, rec_step, vid
         real,    intent(in) :: loc_field(NHOR)
         real                :: glob(NLON*NLAT)
         real(kind=4)        :: out_grid(NLON, NLAT)
         integer             :: jlat, jlon, idx, ret

         call mpgagp(glob, loc_field, 1)

         if (mypid == NROOT) then
            do jlat = 1, NLAT
               do jlon = 1, NLON
                  idx = (jlat - 1) * NLON + jlon
                  out_grid(jlon, jlat) = real(glob(idx), kind=4)
               enddo
            enddo
            ret = nf90_put_var(target_ncid, vid, out_grid, start=(/ 1, 1, rec_step /), &
                               count=(/ NLON, NLAT, 1 /))
            call nc_check(ret, 'put_var 2D')
         endif
      end subroutine write_field_2d_file

      ! ------------------------------------------------------------------
      !  write_field_2d: gather and write 2D variable at current record
      ! ------------------------------------------------------------------
      subroutine write_field_2d(vid, loc_field)
         integer, intent(in) :: vid
         real,    intent(in) :: loc_field(NHOR)
         call write_field_2d_file(ncid, nc_time_step, vid, loc_field)
      end subroutine write_field_2d

      ! ------------------------------------------------------------------
      !  write_field_3d: gather and write 3D variable at current record
      ! ------------------------------------------------------------------
      subroutine write_field_3d(vid, loc_field)
         integer, intent(in) :: vid
         real,    intent(in) :: loc_field(NHOR, NLEV)
         real                :: glob(NLON*NLAT, NLEV)
         real(kind=4)        :: out_grid(NLON, NLAT, NLEV)
         integer             :: jlev, jlat, jlon, idx, ret

         call mpgagp(glob, loc_field, NLEV)

         if (mypid == NROOT) then
            do jlev = 1, NLEV
               do jlat = 1, NLAT
                  do jlon = 1, NLON
                     idx = (jlat - 1) * NLON + jlon
                     out_grid(jlon, jlat, jlev) = real(glob(idx, jlev), kind=4)
                  enddo
               enddo
            enddo
            ret = nf90_put_var(ncid, vid, out_grid, start=(/ 1, 1, 1, nc_time_step /), &
                               count=(/ NLON, NLAT, NLEV, 1 /))
            call nc_check(ret, 'put_var 3D')
         endif
      end subroutine write_field_3d

      ! ------------------------------------------------------------------
      !  netcdf_write_record: scale accumulators and write record
      ! ------------------------------------------------------------------
      subroutine netcdf_write_record
         real(kind=8) :: time_val
         real         :: zscale
         real         :: tmp_2d(NHOR)
         integer      :: ret

         if (nc_acc_count <= 0) return

         nc_time_step = nc_time_step + 1

         ! Write time coordinate
         if (mypid == NROOT) then
            time_val = (nstep - nstep1 + 1) / real(ntspd, kind=8)
            ret = nf90_put_var(ncid, vid_time, (/ time_val /), &
                               start=(/ nc_time_step /), count=(/ 1 /))
            call nc_check(ret, 'put_var time')
         endif

         if (noutput_avg == 0) then
            ! Instantaneous step mode: use current state variables
            call write_field_3d(vid_ua,     du(:,:))
            call write_field_3d(vid_va,     dv(:,:))
            call write_field_3d(vid_ta,     dt(:,1:NLEV))
            call write_field_3d(vid_hus,    dq(:,1:NLEV))
            call write_field_3d(vid_clw,    dql(:,1:NLEV))
            call write_field_3d(vid_cl,     dcc(:,1:NLEV))

            call write_field_2d(vid_ps,     dp(:))
            call write_field_2d(vid_ts,     dt(:,NLEP))
            call write_field_2d(vid_tas,    dtsa(:))
            call write_field_2d(vid_tasmax, acc_tasmax(:))
            call write_field_2d(vid_tasmin, acc_tasmin(:))
            call write_field_2d(vid_mrso,   dwatc(:))
            call write_field_2d(vid_snd,    dsnow(:))
            call write_field_2d(vid_tso,    dtsoil(:))
            call write_field_2d(vid_tso2,   dtd2(:))
            call write_field_2d(vid_tso3,   dtd3(:))
            call write_field_2d(vid_tso4,   dtd4(:))
            call write_field_2d(vid_tsod,   dtd5(:))
            call write_field_2d(vid_sic,    dicec(:))
            call write_field_2d(vid_sit,    diced(:))
            call write_field_2d(vid_mld,    dmld(:))

            ! Fluxes averaged over interval
            zscale = 1.0 / real(nc_acc_count)
            call write_field_2d(vid_prl,    acc_prl * zscale)
            call write_field_2d(vid_prc,    acc_prc * zscale)
            tmp_2d(:) = (acc_prl(:) + acc_prc(:)) * zscale
            call write_field_2d(vid_pr,     tmp_2d)
            call write_field_2d(vid_prsn,   acc_prsn * zscale)
            call write_field_2d(vid_evap,   acc_evap * zscale)
            call write_field_2d(vid_mrro,   acc_mrro * zscale)
            call write_field_2d(vid_hfss,   acc_hfss * zscale)
            call write_field_2d(vid_hfls,   acc_hfls * zscale)
            call write_field_2d(vid_rss,    acc_rss * zscale)
            call write_field_2d(vid_rls,    acc_rls * zscale)
            call write_field_2d(vid_rst,    acc_rst * zscale)
            call write_field_2d(vid_rlut,   acc_rlut * zscale)
            call write_field_2d(vid_rsut,   acc_rsut * zscale)
            call write_field_2d(vid_ssru,   acc_ssru * zscale)
            call write_field_2d(vid_stru,   acc_stru * zscale)
            call write_field_2d(vid_tauu,   acc_tauu * zscale)
            call write_field_2d(vid_tauv,   acc_tauv * zscale)
            call write_field_2d(vid_clt,    acc_clt * zscale)
            call write_field_2d(vid_u3,     acc_u3 * zscale)
            call write_field_2d(vid_snm,    acc_snm * zscale)
            call write_field_2d(vid_sndc,   acc_sndc * zscale)
            call write_field_2d(vid_prw,    acc_prw * zscale)
         else
            ! Averaged mode: scale all accumulated fields by 1/count
            zscale = 1.0 / real(nc_acc_count)

            call write_field_3d(vid_ua,     acc_u * zscale)
            call write_field_3d(vid_va,     acc_v * zscale)
            call write_field_3d(vid_ta,     acc_t * zscale)
            call write_field_3d(vid_hus,    acc_q * zscale)
            call write_field_3d(vid_clw,    acc_clw * zscale)
            call write_field_3d(vid_cl,     acc_cl * zscale)

            call write_field_2d(vid_ps,     acc_ps * zscale)
            call write_field_2d(vid_ts,     acc_ts * zscale)
            call write_field_2d(vid_tas,    acc_tas * zscale)
            call write_field_2d(vid_tasmax, acc_tasmax)
            call write_field_2d(vid_tasmin, acc_tasmin)
            call write_field_2d(vid_mrso,   acc_mrso * zscale)
            call write_field_2d(vid_snd,    acc_snd * zscale)
            call write_field_2d(vid_prl,    acc_prl * zscale)
            call write_field_2d(vid_prc,    acc_prc * zscale)
            tmp_2d(:) = (acc_prl(:) + acc_prc(:)) * zscale
            call write_field_2d(vid_pr,     tmp_2d)
            call write_field_2d(vid_prsn,   acc_prsn * zscale)
            call write_field_2d(vid_evap,   acc_evap * zscale)
            call write_field_2d(vid_mrro,   acc_mrro * zscale)
            call write_field_2d(vid_hfss,   acc_hfss * zscale)
            call write_field_2d(vid_hfls,   acc_hfls * zscale)
            call write_field_2d(vid_rss,    acc_rss * zscale)
            call write_field_2d(vid_rls,    acc_rls * zscale)
            call write_field_2d(vid_rst,    acc_rst * zscale)
            call write_field_2d(vid_rlut,   acc_rlut * zscale)
            call write_field_2d(vid_rsut,   acc_rsut * zscale)
            call write_field_2d(vid_ssru,   acc_ssru * zscale)
            call write_field_2d(vid_stru,   acc_stru * zscale)
            call write_field_2d(vid_tauu,   acc_tauu * zscale)
            call write_field_2d(vid_tauv,   acc_tauv * zscale)
            call write_field_2d(vid_clt,    acc_clt * zscale)
            call write_field_2d(vid_u3,     acc_u3 * zscale)
            call write_field_2d(vid_tso,    acc_tso * zscale)
            call write_field_2d(vid_tso2,   acc_tso2 * zscale)
            call write_field_2d(vid_tso3,   acc_tso3 * zscale)
            call write_field_2d(vid_tso4,   acc_tso4 * zscale)
            call write_field_2d(vid_tsod,   acc_tsod * zscale)
            call write_field_2d(vid_sic,    acc_sic * zscale)
            call write_field_2d(vid_sit,    acc_sit * zscale)
            call write_field_2d(vid_mld,    acc_mld * zscale)
            call write_field_2d(vid_snm,    acc_snm * zscale)
            call write_field_2d(vid_sndc,   acc_sndc * zscale)
            call write_field_2d(vid_prw,    acc_prw * zscale)
         endif

         ! Static fields (invariant)
         call write_field_2d(vid_lsm,    dls(:))
         call write_field_2d(vid_z0,     dz0(:))
         call write_field_2d(vid_as,     dalb(:))
         call write_field_2d(vid_vegf,   dforest(:))
         call write_field_2d(vid_glac,   dglac(:))
         call write_field_2d(vid_dwmax,  dwmax(:))
         call write_field_2d(vid_sg,     oro_gp(:))

         if (mypid == NROOT) then
            ret = nf90_sync(ncid)
            call nc_check(ret, 'sync')
         endif

         call reset_accumulators
      end subroutine netcdf_write_record

      ! ------------------------------------------------------------------
      !  netcdf_step: accumulate and conditionally write record
      ! ------------------------------------------------------------------
      subroutine netcdf_step
         integer :: nmin, nhour, nday, next_mon, next_year
         logical :: lwrite

         if (nnetcdf == 0 .or. noutput == 0) return

         ! Accumulate 3D prognostic/diagnostic fields
         acc_u(:,:)   = acc_u(:,:)   + du(:,:)
         acc_v(:,:)   = acc_v(:,:)   + dv(:,:)
         acc_t(:,:)   = acc_t(:,:)   + dt(:,1:NLEV)
         acc_q(:,:)   = acc_q(:,:)   + dq(:,1:NLEV)
         acc_clw(:,:) = acc_clw(:,:) + dql(:,1:NLEV)
         acc_cl(:,:)  = acc_cl(:,:)  + dcc(:,1:NLEV)

         ! Accumulate 2D state fields
         acc_ps(:)    = acc_ps(:)    + dp(:)
         acc_ts(:)    = acc_ts(:)    + dt(:,NLEP)
         acc_tas(:)   = acc_tas(:)   + dtsa(:)
         acc_tasmax(:)= max(acc_tasmax(:), dtsa(:))
         acc_tasmin(:)= min(acc_tasmin(:), dtsa(:))
         acc_mrso(:)  = acc_mrso(:)  + dwatc(:)
         acc_snd(:)   = acc_snd(:)   + dsnow(:)
         acc_tso(:)   = acc_tso(:)   + dtsoil(:)
         acc_tso2(:)  = acc_tso2(:)  + dtd2(:)
         acc_tso3(:)  = acc_tso3(:)  + dtd3(:)
         acc_tso4(:)  = acc_tso4(:)  + dtd4(:)
         acc_tsod(:)  = acc_tsod(:)  + dtd5(:)
         acc_sic(:)   = acc_sic(:)   + dicec(:)
         acc_sit(:)   = acc_sit(:)   + diced(:)
         acc_mld(:)   = acc_mld(:)   + dmld(:)

         ! Accumulate fluxes
         acc_prl(:)   = acc_prl(:)   + dprl(:)
         acc_prc(:)   = acc_prc(:)   + dprc(:)
         acc_prsn(:)  = acc_prsn(:)  + dprs(:)
         acc_evap(:)  = acc_evap(:)  + devap(:)
         acc_mrro(:)  = acc_mrro(:)  + drunoff(:)
         acc_hfss(:)  = acc_hfss(:)  + dshfl(:)
         acc_hfls(:)  = acc_hfls(:)  + dlhfl(:)
         acc_rss(:)   = acc_rss(:)   + dswfl(:,NLEP)
         acc_rls(:)   = acc_rls(:)   + dlwfl(:,NLEP)
         acc_rst(:)   = acc_rst(:)   + dswfl(:,1)
         acc_rlut(:)  = acc_rlut(:)  + dlwfl(:,1)
         acc_rsut(:)  = acc_rsut(:)  + dfu(:,1)
         acc_ssru(:)  = acc_ssru(:)  + dfu(:,NLEP)
         acc_stru(:)  = acc_stru(:)  + dftu(:,NLEP)
         acc_tauu(:)  = acc_tauu(:)  + dtaux(:)
         acc_tauv(:)  = acc_tauv(:)  + dtauy(:)
         acc_clt(:)   = acc_clt(:)   + dcc(:,NLEP)
         acc_u3(:)    = acc_u3(:)    + dust3(:)
         acc_snm(:)   = acc_snm(:)   + dsmelt(:)
         acc_sndc(:)  = acc_sndc(:)  + dsndch(:)
         acc_prw(:)   = acc_prw(:)   + dqvi(:)

         nc_acc_count = nc_acc_count + 1

         ! Check if it is time to write record
         lwrite = .false.
         if (noutput_avg == 2) then
            ! Monthly average: write on the last step of the month
            call ntomin(nstep + 1, nmin, nhour, nday, next_mon, next_year)
            if (next_mon /= ndatim(2)) lwrite = .true.
         else if (noutput_avg == 1) then
            ! Daily average: write on the last step of the day
            if (mod(nstep - nstep1 + 1, ntspd) == 0) lwrite = .true.
         else
            ! Step / instantaneous: write every nafter steps
            if (mod(nstep - nstep1 + 1, nafter) == 0) lwrite = .true.
         endif

         if (lwrite) call netcdf_write_record
      end subroutine netcdf_step

      ! ------------------------------------------------------------------
      !  netcdf_close: write any pending accumulated record and close file
      ! ------------------------------------------------------------------
      subroutine netcdf_close
         integer :: ret

         if (nnetcdf == 0 .or. noutput == 0) return

         ! Write any partial average if simulation ends mid-period
         if (nc_acc_count > 0) call netcdf_write_record

         if (mypid == NROOT .and. ncid >= 0) then
            ret = nf90_close(ncid)
            call nc_check(ret, 'close')
            ncid = -1
         endif

         if (allocated(acc_u)) then
            deallocate(acc_u, acc_v, acc_t, acc_q, acc_clw, acc_cl)
            deallocate(acc_ps, acc_ts, acc_tas, acc_tasmax, acc_tasmin)
            deallocate(acc_mrso, acc_snd, acc_prl, acc_prc, acc_prsn)
            deallocate(acc_evap, acc_mrro, acc_hfss, acc_hfls)
            deallocate(acc_rss, acc_rls, acc_rst, acc_rlut)
            deallocate(acc_rsut, acc_ssru, acc_stru)
            deallocate(acc_tauu, acc_tauv, acc_clt, acc_u3)
            deallocate(acc_tso, acc_tso2, acc_tso3, acc_tso4, acc_tsod)
            deallocate(acc_sic, acc_sit, acc_mld, acc_snm, acc_sndc)
            deallocate(acc_prw, oro_gp)
         endif

         call netcdf_ocean_close
         call netcdf_ice_close
      end subroutine netcdf_close

      ! ==================================================================
      !  Mixed-Layer Ocean NetCDF Output Subroutines
      ! ==================================================================

      subroutine inq_ocean_var_ids
         integer :: ret
         ret = nf90_inq_varid(ncid_ocean, 'time',   vid_time_oce)
         ret = nf90_inq_varid(ncid_ocean, 'heata',  vid_oce_heata)
         ret = nf90_inq_varid(ncid_ocean, 'ifluxa', vid_oce_ifluxa)
         ret = nf90_inq_varid(ncid_ocean, 'fssta',  vid_oce_fssta)
         ret = nf90_inq_varid(ncid_ocean, 'dssta',  vid_oce_dssta)
         ret = nf90_inq_varid(ncid_ocean, 'qhda',   vid_oce_qhda)
         ret = nf90_inq_varid(ncid_ocean, 'fldoa',  vid_oce_fldoa)
         ret = nf90_inq_varid(ncid_ocean, 'icec',   vid_oce_icec)
         ret = nf90_inq_varid(ncid_ocean, 'sst',    vid_oce_sst)
         ret = nf90_inq_varid(ncid_ocean, 'ls',     vid_oce_ls)
         ret = nf90_inq_varid(ncid_ocean, 'clsst',  vid_oce_clsst)
      end subroutine inq_ocean_var_ids

      subroutine netcdf_ocean_init(youtput)
         character(len=*), intent(in) :: youtput
         character(len=256)           :: fname, time_units
         real(kind=4)                 :: coord_lon(NLON), coord_lat(NLAT)
         integer                      :: jlon, jlat, ret, did_time, existing_len

         if (nnetcdf == 0 .or. noutput == 0) return

         ocean_nc_time_step = 0

         if (mypid /= NROOT) return

         fname = trim(youtput)
         if (len_trim(fname) == 0) fname = 'ocean_output'
         if (index(fname, '.nc') == 0) then
            fname = trim(fname) // '.nc'
         endif

         if (nrestart > 0) then
            ret = nf90_open(trim(fname), NF90_WRITE, ncid_ocean)
            if (ret == NF90_NOERR) then
               ret = nf90_inq_dimid(ncid_ocean, 'time', did_time)
               call nc_check(ret, 'inq_dimid time ocean')
               ret = nf90_inquire_dimension(ncid_ocean, did_time, len=existing_len)
               call nc_check(ret, 'inquire_dim time ocean')
               ocean_nc_time_step = existing_len
               call inq_ocean_var_ids
               return
            endif
         endif

         ret = nf90_create(trim(fname), IOR(NF90_CLOBBER, NF90_NETCDF4), ncid_ocean)
         call nc_check(ret, 'create ' // trim(fname))

         ret = nf90_put_att(ncid_ocean, NF90_GLOBAL, 'title', 'PlaSim-LSG TNG Mixed Layer Ocean NetCDF Output')
         ret = nf90_put_att(ncid_ocean, NF90_GLOBAL, 'Conventions', 'CF-1.6')
         ret = nf90_put_att(ncid_ocean, NF90_GLOBAL, 'source', 'PlaSim-LSG TNG')

         ret = nf90_def_dim(ncid_ocean, 'lon', NLON, dim_lon_oce)
         call nc_check(ret, 'def_dim lon ocean')
         ret = nf90_def_dim(ncid_ocean, 'lat', NLAT, dim_lat_oce)
         call nc_check(ret, 'def_dim lat ocean')
         ret = nf90_def_dim(ncid_ocean, 'time', NF90_UNLIMITED, dim_time_oce)
         call nc_check(ret, 'def_dim time ocean')

         ret = nf90_def_var(ncid_ocean, 'lon', NF90_FLOAT, (/ dim_lon_oce /), vid_lon_oce)
         call nc_check(ret, 'def_var lon ocean')
         ret = nf90_put_att(ncid_ocean, vid_lon_oce, 'standard_name', 'longitude')
         ret = nf90_put_att(ncid_ocean, vid_lon_oce, 'long_name', 'longitude')
         ret = nf90_put_att(ncid_ocean, vid_lon_oce, 'units', 'degrees_east')
         ret = nf90_put_att(ncid_ocean, vid_lon_oce, 'axis', 'X')

         ret = nf90_def_var(ncid_ocean, 'lat', NF90_FLOAT, (/ dim_lat_oce /), vid_lat_oce)
         call nc_check(ret, 'def_var lat ocean')
         ret = nf90_put_att(ncid_ocean, vid_lat_oce, 'standard_name', 'latitude')
         ret = nf90_put_att(ncid_ocean, vid_lat_oce, 'long_name', 'latitude')
         ret = nf90_put_att(ncid_ocean, vid_lat_oce, 'units', 'degrees_north')
         ret = nf90_put_att(ncid_ocean, vid_lat_oce, 'axis', 'Y')

         ret = nf90_def_var(ncid_ocean, 'time', NF90_DOUBLE, (/ dim_time_oce /), vid_time_oce)
         call nc_check(ret, 'def_var time ocean')
         ret = nf90_put_att(ncid_ocean, vid_time_oce, 'standard_name', 'time')
         ret = nf90_put_att(ncid_ocean, vid_time_oce, 'long_name', 'time')
         write(time_units, '("days since ",I4.4,"-",I2.2,"-01 00:00:00")') n_start_year, n_start_month
         ret = nf90_put_att(ncid_ocean, vid_time_oce, 'units', trim(time_units))
         if (n_days_per_year == 360) then
            ret = nf90_put_att(ncid_ocean, vid_time_oce, 'calendar', '360_day')
         else if (n_days_per_year == 365) then
            ret = nf90_put_att(ncid_ocean, vid_time_oce, 'calendar', '365_day')
         else
            ret = nf90_put_att(ncid_ocean, vid_time_oce, 'calendar', 'proleptic_gregorian')
         endif
         ret = nf90_put_att(ncid_ocean, vid_time_oce, 'axis', 'T')

         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'heata', 'Heat flux from atm/ice', 'W/m2', 901, vid_oce_heata)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'ifluxa', 'Heat flux into ice', 'W/m2', 902, vid_oce_ifluxa)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'fssta', 'Flux correction', 'W/m2', 903, vid_oce_fssta)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'dssta', 'Heat flux from vertical diffusion', 'W/m2', 904, vid_oce_dssta)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'qhda', 'Heat flux from horizontal diffusion', 'W/m2', 905, vid_oce_qhda)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'fldoa', 'Flux from deep ocean (LSG)', 'W/m2', 906, vid_oce_fldoa)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'icec', 'Sea ice cover', 'fraction', 910, vid_oce_icec)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'sst', 'SST (ML ocean temperature)', 'K', 939, vid_oce_sst)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'ls', 'Land-sea mask', 'fraction', 972, vid_oce_ls)
         call def_var_2d_file(ncid_ocean, dim_lon_oce, dim_lat_oce, dim_time_oce, &
                              'clsst', 'Climatological SST', 'K', 990, vid_oce_clsst)

         ret = nf90_enddef(ncid_ocean)
         call nc_check(ret, 'enddef ocean')

         do jlon = 1, NLON
            coord_lon(jlon) = real((jlon - 1) * 360.0 / NLON, kind=4)
         enddo
         ret = nf90_put_var(ncid_ocean, vid_lon_oce, coord_lon)
         call nc_check(ret, 'put_var lon ocean')

         do jlat = 1, NLAT
            coord_lat(jlat) = real(deglat(jlat), kind=4)
         enddo
         ret = nf90_put_var(ncid_ocean, vid_lat_oce, coord_lat)
         call nc_check(ret, 'put_var lat ocean')

         ret = nf90_sync(ncid_ocean)
         call nc_check(ret, 'sync ocean')
      end subroutine netcdf_ocean_init

      subroutine netcdf_ocean_write(pheata, pifluxa, pfssta, pdssta, &
                                    pqhda, pfldoa, picec, psst, pls, pclsst, kstep, ktspd)
         real,    intent(in) :: pheata(NHOR), pifluxa(NHOR), pfssta(NHOR), pdssta(NHOR)
         real,    intent(in) :: pqhda(NHOR), pfldoa(NHOR), picec(NHOR), psst(NHOR)
         real,    intent(in) :: pls(NHOR), pclsst(NHOR)
         integer, intent(in) :: kstep, ktspd
         real(kind=8)        :: time_val
         integer             :: ret

         if (nnetcdf == 0 .or. noutput == 0) return

         ocean_nc_time_step = ocean_nc_time_step + 1

         if (mypid == NROOT) then
            time_val = kstep / real(ktspd, kind=8)
            ret = nf90_put_var(ncid_ocean, vid_time_oce, (/ time_val /), &
                               start=(/ ocean_nc_time_step /), count=(/ 1 /))
            call nc_check(ret, 'put_var time ocean')
         endif

         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_heata,  pheata)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_ifluxa, pifluxa)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_fssta,  pfssta)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_dssta,  pdssta)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_qhda,   pqhda)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_fldoa,  pfldoa)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_icec,   picec)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_sst,    psst)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_ls,     pls)
         call write_field_2d_file(ncid_ocean, ocean_nc_time_step, vid_oce_clsst,  pclsst)

         if (mypid == NROOT) then
            ret = nf90_sync(ncid_ocean)
            call nc_check(ret, 'sync ocean')
         endif
      end subroutine netcdf_ocean_write

      subroutine netcdf_ocean_close
         integer :: ret
         if (nnetcdf == 0 .or. noutput == 0) return
         if (mypid == NROOT .and. ncid_ocean >= 0) then
            ret = nf90_close(ncid_ocean)
            call nc_check(ret, 'close ocean')
            ncid_ocean = -1
         endif
         ocean_nc_time_step = 0
      end subroutine netcdf_ocean_close

      ! ==================================================================
      !  Sea Ice NetCDF Output Subroutines
      ! ==================================================================

      subroutine inq_ice_var_ids
         integer :: ret
         ret = nf90_inq_varid(ncid_ice, 'time',     vid_time_ice)
         ret = nf90_inq_varid(ncid_ice, 'heata',    vid_ice_heata)
         ret = nf90_inq_varid(ncid_ice, 'ofluxa',   vid_ice_ofluxa)
         ret = nf90_inq_varid(ncid_ice, 'tsfluxa',  vid_ice_tsfluxa)
         ret = nf90_inq_varid(ncid_ice, 'smelta',   vid_ice_smelta)
         ret = nf90_inq_varid(ncid_ice, 'imelta',   vid_ice_imelta)
         ret = nf90_inq_varid(ncid_ice, 'cfluxa',   vid_ice_cfluxa)
         ret = nf90_inq_varid(ncid_ice, 'fluxca',   vid_ice_fluxca)
         ret = nf90_inq_varid(ncid_ice, 'qmelta',   vid_ice_qmelta)
         ret = nf90_inq_varid(ncid_ice, 'xflxicea', vid_ice_xflxicea)
         ret = nf90_inq_varid(ncid_ice, 'icec',     vid_ice_icec)
         ret = nf90_inq_varid(ncid_ice, 'iced',     vid_ice_iced)
         ret = nf90_inq_varid(ncid_ice, 'scflxa',   vid_ice_scflxa)
         ret = nf90_inq_varid(ncid_ice, 'cfluxra',  vid_ice_cfluxra)
         ret = nf90_inq_varid(ncid_ice, 'cfluxna',  vid_ice_cfluxna)
         ret = nf90_inq_varid(ncid_ice, 'ts',       vid_ice_ts)
         ret = nf90_inq_varid(ncid_ice, 'zsnow',    vid_ice_zsnow)
         ret = nf90_inq_varid(ncid_ice, 'sst',      vid_ice_sst)
         ret = nf90_inq_varid(ncid_ice, 'ls',       vid_ice_ls)
         ret = nf90_inq_varid(ncid_ice, 'clicec2',  vid_ice_clicec2)
         ret = nf90_inq_varid(ncid_ice, 'cliced2',  vid_ice_cliced2)
         ret = nf90_inq_varid(ncid_ice, 'icecc',    vid_ice_icecc)
         ret = nf90_inq_varid(ncid_ice, 'cpmea',    vid_ice_cpmea)
         ret = nf90_inq_varid(ncid_ice, 'croffa',   vid_ice_croffa)
         ret = nf90_inq_varid(ncid_ice, 'stoia',    vid_ice_stoia)
      end subroutine inq_ice_var_ids

      subroutine netcdf_ice_init(youtput)
         character(len=*), intent(in) :: youtput
         character(len=256)           :: fname, time_units
         real(kind=4)                 :: coord_lon(NLON), coord_lat(NLAT)
         integer                      :: jlon, jlat, ret, did_time, existing_len

         if (nnetcdf == 0 .or. noutput == 0) return

         ice_nc_time_step = 0

         if (mypid /= NROOT) return

         fname = trim(youtput)
         if (len_trim(fname) == 0) fname = 'ice_output'
         if (index(fname, '.nc') == 0) then
            fname = trim(fname) // '.nc'
         endif

         if (nrestart > 0) then
            ret = nf90_open(trim(fname), NF90_WRITE, ncid_ice)
            if (ret == NF90_NOERR) then
               ret = nf90_inq_dimid(ncid_ice, 'time', did_time)
               call nc_check(ret, 'inq_dimid time ice')
               ret = nf90_inquire_dimension(ncid_ice, did_time, len=existing_len)
               call nc_check(ret, 'inquire_dim time ice')
               ice_nc_time_step = existing_len
               call inq_ice_var_ids
               return
            endif
         endif

         ret = nf90_create(trim(fname), IOR(NF90_CLOBBER, NF90_NETCDF4), ncid_ice)
         call nc_check(ret, 'create ' // trim(fname))

         ret = nf90_put_att(ncid_ice, NF90_GLOBAL, 'title', 'PlaSim-LSG TNG Sea Ice NetCDF Output')
         ret = nf90_put_att(ncid_ice, NF90_GLOBAL, 'Conventions', 'CF-1.6')
         ret = nf90_put_att(ncid_ice, NF90_GLOBAL, 'source', 'PlaSim-LSG TNG')

         ret = nf90_def_dim(ncid_ice, 'lon', NLON, dim_lon_ice)
         call nc_check(ret, 'def_dim lon ice')
         ret = nf90_def_dim(ncid_ice, 'lat', NLAT, dim_lat_ice)
         call nc_check(ret, 'def_dim lat ice')
         ret = nf90_def_dim(ncid_ice, 'time', NF90_UNLIMITED, dim_time_ice)
         call nc_check(ret, 'def_dim time ice')

         ret = nf90_def_var(ncid_ice, 'lon', NF90_FLOAT, (/ dim_lon_ice /), vid_lon_ice)
         call nc_check(ret, 'def_var lon ice')
         ret = nf90_put_att(ncid_ice, vid_lon_ice, 'standard_name', 'longitude')
         ret = nf90_put_att(ncid_ice, vid_lon_ice, 'long_name', 'longitude')
         ret = nf90_put_att(ncid_ice, vid_lon_ice, 'units', 'degrees_east')
         ret = nf90_put_att(ncid_ice, vid_lon_ice, 'axis', 'X')

         ret = nf90_def_var(ncid_ice, 'lat', NF90_FLOAT, (/ dim_lat_ice /), vid_lat_ice)
         call nc_check(ret, 'def_var lat ice')
         ret = nf90_put_att(ncid_ice, vid_lat_ice, 'standard_name', 'latitude')
         ret = nf90_put_att(ncid_ice, vid_lat_ice, 'long_name', 'latitude')
         ret = nf90_put_att(ncid_ice, vid_lat_ice, 'units', 'degrees_north')
         ret = nf90_put_att(ncid_ice, vid_lat_ice, 'axis', 'Y')

         ret = nf90_def_var(ncid_ice, 'time', NF90_DOUBLE, (/ dim_time_ice /), vid_time_ice)
         call nc_check(ret, 'def_var time ice')
         ret = nf90_put_att(ncid_ice, vid_time_ice, 'standard_name', 'time')
         ret = nf90_put_att(ncid_ice, vid_time_ice, 'long_name', 'time')
         write(time_units, '("days since ",I4.4,"-",I2.2,"-01 00:00:00")') n_start_year, n_start_month
         ret = nf90_put_att(ncid_ice, vid_time_ice, 'units', trim(time_units))
         if (n_days_per_year == 360) then
            ret = nf90_put_att(ncid_ice, vid_time_ice, 'calendar', '360_day')
         else if (n_days_per_year == 365) then
            ret = nf90_put_att(ncid_ice, vid_time_ice, 'calendar', '365_day')
         else
            ret = nf90_put_att(ncid_ice, vid_time_ice, 'calendar', 'proleptic_gregorian')
         endif
         ret = nf90_put_att(ncid_ice, vid_time_ice, 'axis', 'T')

         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'heata', 'Flux from the atmosphere', 'W/m2', 701, vid_ice_heata)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'ofluxa', 'Flux from the ocean', 'W/m2', 702, vid_ice_ofluxa)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'tsfluxa', 'Flux to warm/cool ice/snow', 'W/m2', 703, vid_ice_tsfluxa)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'smelta', 'Flux for snow melt', 'W/m2', 704, vid_ice_smelta)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'imelta', 'Flux used for ice melt', 'W/m2', 705, vid_ice_imelta)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'cfluxa', 'Flux to the ocean', 'W/m2', 706, vid_ice_cfluxa)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'fluxca', 'Cond. heatflux', 'W/m2', 707, vid_ice_fluxca)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'qmelta', 'Res flux to ice', 'W/m2', 708, vid_ice_qmelta)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'xflxicea', 'Flux correction', 'W/m2', 709, vid_ice_xflxicea)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'icec', 'Ice cover', 'fraction', 710, vid_ice_icec)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'iced', 'Ice thickness', 'm', 711, vid_ice_iced)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'scflxa', 'Flux from snow -> ice conversion', 'W/m2', 712, vid_ice_scflxa)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'cfluxra', 'Flux for limiting ice to xmaxd', 'W/m2', 713, vid_ice_cfluxra)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'cfluxna', 'Flux due to neg. ice', 'W/m2', 714, vid_ice_cfluxna)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'ts', 'Surface temperature', 'K', 739, vid_ice_ts)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'zsnow', 'Snow depth', 'm', 741, vid_ice_zsnow)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'sst', 'Sea surface temperature', 'K', 769, vid_ice_sst)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'ls', 'Land sea mask', 'fraction', 772, vid_ice_ls)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'clicec2', 'Climatological ice cover', 'fraction', 790, vid_ice_clicec2)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'cliced2', 'Climatological ice thickness', 'm', 791, vid_ice_cliced2)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'icecc', 'Ice cover computed prognostically', 'fraction', 792, vid_ice_icecc)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'cpmea', 'Fresh water (p-e) for LSG', 'm/s', 794, vid_ice_cpmea)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'croffa', 'Runoff for LSG', 'm/s', 795, vid_ice_croffa)
         call def_var_2d_file(ncid_ice, dim_lon_ice, dim_lat_ice, dim_time_ice, &
                              'stoia', 'Snow converted to ice', 'm', 796, vid_ice_stoia)

         ret = nf90_enddef(ncid_ice)
         call nc_check(ret, 'enddef ice')

         do jlon = 1, NLON
            coord_lon(jlon) = real((jlon - 1) * 360.0 / NLON, kind=4)
         enddo
         ret = nf90_put_var(ncid_ice, vid_lon_ice, coord_lon)
         call nc_check(ret, 'put_var lon ice')

         do jlat = 1, NLAT
            coord_lat(jlat) = real(deglat(jlat), kind=4)
         enddo
         ret = nf90_put_var(ncid_ice, vid_lat_ice, coord_lat)
         call nc_check(ret, 'put_var lat ice')

         ret = nf90_sync(ncid_ice)
         call nc_check(ret, 'sync ice')
      end subroutine netcdf_ice_init

      subroutine netcdf_ice_write(pxheata, pxofluxa, pxtsfluxa, pxsmelta, pximelta, &
                                  pxcfluxa, pxfluxca, pxqmelta, pxflxicea, pxicec,  &
                                  pxiced, pxscflxa, pxcfluxra, pxcfluxna, pxts,     &
                                  pzsnow, pxsst, pxls, pxclicec2, pxcliced2,        &
                                  pxicecc, pxcpmea, pxcroffa, pxstoia, kstep, ktspd)
         real,    intent(in) :: pxheata(NHOR), pxofluxa(NHOR), pxtsfluxa(NHOR)
         real,    intent(in) :: pxsmelta(NHOR), pximelta(NHOR), pxcfluxa(NHOR)
         real,    intent(in) :: pxfluxca(NHOR), pxqmelta(NHOR), pxflxicea(NHOR)
         real,    intent(in) :: pxicec(NHOR), pxiced(NHOR), pxscflxa(NHOR)
         real,    intent(in) :: pxcfluxra(NHOR), pxcfluxna(NHOR), pxts(NHOR)
         real,    intent(in) :: pzsnow(NHOR), pxsst(NHOR), pxls(NHOR)
         real,    intent(in) :: pxclicec2(NHOR), pxcliced2(NHOR), pxicecc(NHOR)
         real,    intent(in) :: pxcpmea(NHOR), pxcroffa(NHOR), pxstoia(NHOR)
         integer, intent(in) :: kstep, ktspd
         real(kind=8)        :: time_val
         integer             :: ret

         if (nnetcdf == 0 .or. noutput == 0) return

         ice_nc_time_step = ice_nc_time_step + 1

         if (mypid == NROOT) then
            time_val = kstep / real(ktspd, kind=8)
            ret = nf90_put_var(ncid_ice, vid_time_ice, (/ time_val /), &
                               start=(/ ice_nc_time_step /), count=(/ 1 /))
            call nc_check(ret, 'put_var time ice')
         endif

         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_heata,    pxheata)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_ofluxa,   pxofluxa)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_tsfluxa,  pxtsfluxa)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_smelta,   pxsmelta)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_imelta,   pximelta)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_cfluxa,   pxcfluxa)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_fluxca,   pxfluxca)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_qmelta,   pxqmelta)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_xflxicea, pxflxicea)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_icec,     pxicec)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_iced,     pxiced)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_scflxa,   pxscflxa)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_cfluxra,  pxcfluxra)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_cfluxna,  pxcfluxna)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_ts,       pxts)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_zsnow,    pzsnow)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_sst,      pxsst)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_ls,       pxls)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_clicec2,  pxclicec2)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_cliced2,  pxcliced2)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_icecc,    pxicecc)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_cpmea,    pxcpmea)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_croffa,   pxcroffa)
         call write_field_2d_file(ncid_ice, ice_nc_time_step, vid_ice_stoia,    pxstoia)

         if (mypid == NROOT) then
            ret = nf90_sync(ncid_ice)
            call nc_check(ret, 'sync ice')
         endif
      end subroutine netcdf_ice_write

      subroutine netcdf_ice_close
         integer :: ret
         if (nnetcdf == 0 .or. noutput == 0) return
         if (mypid == NROOT .and. ncid_ice >= 0) then
            ret = nf90_close(ncid_ice)
            call nc_check(ret, 'close ice')
            ncid_ice = -1
         endif
         ice_nc_time_step = 0
      end subroutine netcdf_ice_close

      end module netcdfmod
