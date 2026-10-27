!     ======================
!     MODULE NETCDFMOD (STUB)
!     ======================

      module netcdfmod
      use pumamod
      implicit none

      contains

      subroutine netcdf_init
         if (nnetcdf > 0) then
            if (mypid == NROOT) then
               write(*,*) 'Error: Direct NetCDF output requested (nnetcdf > 0)'
               write(*,*) '       but model was compiled with netcdfmod_stub.'
               write(*,*) '       Please recompile with NetCDF-Fortran support.'
            endif
            call mpabort('NetCDF output not available in stub build')
         endif
      end subroutine netcdf_init

      subroutine netcdf_step
      end subroutine netcdf_step

      subroutine netcdf_close
      end subroutine netcdf_close

      subroutine netcdf_ocean_init(youtput)
         character(len=*), intent(in) :: youtput
      end subroutine netcdf_ocean_init

      subroutine netcdf_ocean_write(pheata, pifluxa, pfssta, pdssta, &
                                    pqhda, pfldoa, picec, psst, pls, pclsst, kstep, ktspd)
         real,    intent(in) :: pheata(*), pifluxa(*), pfssta(*), pdssta(*)
         real,    intent(in) :: pqhda(*), pfldoa(*), picec(*), psst(*)
         real,    intent(in) :: pls(*), pclsst(*)
         integer, intent(in) :: kstep, ktspd
      end subroutine netcdf_ocean_write

      subroutine netcdf_ocean_close
      end subroutine netcdf_ocean_close

      subroutine netcdf_ice_init(youtput)
         character(len=*), intent(in) :: youtput
      end subroutine netcdf_ice_init

      subroutine netcdf_ice_write(pxheata, pxofluxa, pxtsfluxa, pxsmelta, pximelta, &
                                  pxcfluxa, pxfluxca, pxqmelta, pxflxicea, pxicec,  &
                                  pxiced, pxscflxa, pxcfluxra, pxcfluxna, pxts,     &
                                  pzsnow, pxsst, pxls, pxclicec2, pxcliced2,        &
                                  pxicecc, pxcpmea, pxcroffa, pxstoia, kstep, ktspd)
         real,    intent(in) :: pxheata(*), pxofluxa(*), pxtsfluxa(*)
         real,    intent(in) :: pxsmelta(*), pximelta(*), pxcfluxa(*)
         real,    intent(in) :: pxfluxca(*), pxqmelta(*), pxflxicea(*)
         real,    intent(in) :: pxicec(*), pxiced(*), pxscflxa(*)
         real,    intent(in) :: pxcfluxra(*), pxcfluxna(*), pxts(*)
         real,    intent(in) :: pzsnow(*), pxsst(*), pxls(*)
         real,    intent(in) :: pxclicec2(*), pxcliced2(*), pxicecc(*)
         real,    intent(in) :: pxcpmea(*), pxcroffa(*), pxstoia(*)
         integer, intent(in) :: kstep, ktspd
      end subroutine netcdf_ice_write

      subroutine netcdf_ice_close
      end subroutine netcdf_ice_close

      end module netcdfmod
