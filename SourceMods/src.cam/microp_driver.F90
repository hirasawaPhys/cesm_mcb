module microp_driver

!-------------------------------------------------------------------------------------------------------
!
! Driver for CAM microphysics parameterizations
!
!-------------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pver
use physics_types,  only: physics_state, physics_ptend, physics_tend,  &
                          physics_ptend_copy, physics_ptend_sum
use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
use phys_control,   only: phys_getopts

use micro_mg_cam,   only: micro_mg_cam_readnl, micro_mg_cam_register, &
                          micro_mg_cam_implements_cnst, micro_mg_cam_init_cnst, &
                          micro_mg_cam_init, micro_mg_cam_tend
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use perf_mod,       only: t_startf, t_stopf

implicit none
private
save

public :: &
   microp_driver_readnl,          &
   microp_driver_register,        &
   microp_driver_init_cnst,       &
   microp_driver_implements_cnst, &
   microp_driver_init,            &
   microp_driver_tend,             &
   set_cdnc

character(len=16)  :: microp_scheme   ! Microphysics scheme

!===============================================================================
contains
!===============================================================================

subroutine microp_driver_readnl(nlfile)

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Read in namelist for microphysics scheme
   !-----------------------------------------------------------------------

   call phys_getopts(microp_scheme_out=microp_scheme)

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_readnl(nlfile)
   case ('NONE', 'RK', 'SPCAM_sam1mom', 'SPCAM_m2005')
      continue
   case default
      call endrun('microp_driver_readnl:: unrecognized microp_scheme, "'//trim(microp_scheme)//'"')
   end select

end subroutine microp_driver_readnl

subroutine microp_driver_register

   ! Register microphysics constituents and fields in the physics buffer.
   !-----------------------------------------------------------------------


   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_register()
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_register:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_register

!===============================================================================

function microp_driver_implements_cnst(name)

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: microp_driver_implements_cnst    ! return value

   ! Local workspace
   integer :: m
   !-----------------------------------------------------------------------

   microp_driver_implements_cnst = .false.

   select case (microp_scheme)
   case ('MG')
      microp_driver_implements_cnst = micro_mg_cam_implements_cnst(name)
   case ('NONE', 'RK', 'SPCAM_sam1mom', 'SPCAM_m2005')
      continue
   case default
      call endrun('microp_driver_implements_cnst:: unrecognized microp_scheme, '//trim(microp_scheme))
   end select

end function microp_driver_implements_cnst

!===============================================================================

subroutine microp_driver_init_cnst(name, latvals, lonvals, mask, q)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in)  :: name       ! constituent name
   real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
   real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
   logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev
   !-----------------------------------------------------------------------

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init_cnst(name, latvals, lonvals, mask, q)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case ('SPCAM_m2005')
      ! microp_driver doesn't handle this one
      continue
   case ('SPCAM_sam1mom')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_init_cnst:: unrecognized microp_scheme'//trim(microp_scheme))
   end select

end subroutine microp_driver_init_cnst

!===============================================================================

subroutine microp_driver_init(pbuf2d)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! Initialize the microphysics parameterizations
   !-----------------------------------------------------------------------

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init(pbuf2d)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_init:: unrecognized microp_scheme'//trim(microp_scheme))
   end select


end subroutine microp_driver_init

!===============================================================================

subroutine microp_driver_tend(state, ptend, dtime, pbuf)

   ! Call the microphysics parameterization run methods.

   ! Input arguments

   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend       ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep

   ! Local variables

   integer :: lchnk
   integer :: ncol

   !======================================================================

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Call MG Microphysics

   select case (microp_scheme)
   case ('MG')
      call t_startf('microp_mg_tend')
      call micro_mg_cam_tend(state, ptend, dtime, pbuf)
      call t_stopf('microp_mg_tend')
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_tend:: unrecognized microp_scheme'//trim(microp_scheme))
   end select

end subroutine microp_driver_tend

subroutine set_cdnc(state, ptend, pbuf, nc_reg1, nc_reg2, nc_reg3)

   use constituents,     only: cnst_get_ind
   use physics_buffer,   only: pbuf_old_tim_idx, pbuf_get_index
   use ppgrid,           only: pcols, pver

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state) :: state
   type(physics_ptend) :: ptend
   real(r8)            :: nc_reg1 ! CDNC target for region 1 (NEPac)
   real(r8)            :: nc_reg2 ! CDNC target for region 2 (SEPac)
   real(r8)            :: nc_reg3 ! CDNC target for region 3 (SEAtl)
   integer             :: ixnumliq
   integer             :: itim_old
   integer             :: alst_idx = -1
   real(r8), pointer, dimension(:,:) :: alst        ! cloud fraction
   integer             :: i, k, ncol

   ncol  = state%ncol

   ! Get liquid stratiform cloud fraction
   itim_old = pbuf_old_tim_idx()
   alst_idx  = pbuf_get_index('ALST')
   call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   ! Get index of cloud liquid number concentration
   call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)

   ! set CDNC for columns in specified regions
   do i = 1, ncol  
      !! North East Pacific [0, 30N], [110W, 150W]
      if(state%lat(i) > 0._r8 .and. state%lat(i) < 0.523599_r8 .and. &
         state%lon(i) > 3.665191_r8 .and. state%lon(i) < 4.363323_r8) then
            ! debug write, remove later
            write(iulog,*) 'testing region 1 def lat=',state%lat(i), ', lon=',state%lon(i)
            ! modify over the entire column
            do k = ptend%top_level, ptend%bot_level
               ! q(i,k,ixnumliq) is the grid box average value, so we need to multiply 
               !  by the liquid cloud fraction (alst)
               state%q(i,k,ixnumliq) = nc_reg1 * alst(i,k)
            end do
      end if

      !! South East Pacific [30S, 0], [70W, 110W]
      if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. &
         state%lon(i) < 5.061455_r8 .and. state%lon(i) > 4.363323_r8) then
            write(iulog,*) 'testing region 2 def lat=',state%lat(i), ', lon=',state%lon(i)
            do k = ptend%top_level, ptend%bot_level
               state%q(i,k,ixnumliq) = nc_reg2 * alst(i,k)
            end do
      end if

      !! South East Atlantic [30S, 0], [15E, 25W]
      ! East of 25W
      if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. &
         state%lon(i) > 5.846853_r8) then
            write(iulog,*) 'testing region 3 def lat=',state%lat(i), ', lon=',state%lon(i)
            do k = ptend%top_level, ptend%bot_level
               state%q(i,k,ixnumliq) = nc_reg3 * alst(i,k)
            end do
      end if
      ! note need to split because we're spanning over the 0 meridian - this is for west of 15E
      if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. &
         state%lon(i) < 0.261799_r8) then
            write(iulog,*) 'testing region 3 def lat=',state%lat(i), ', lon=',state%lon(i)
            do k = ptend%top_level, ptend%bot_level
               state%q(i,k,ixnumliq) = nc_reg3 * alst(i,k)
            end do
      end if

   end do

end subroutine set_cdnc

end module microp_driver
