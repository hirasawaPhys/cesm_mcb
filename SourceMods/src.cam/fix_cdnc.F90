module fix_cdnc
    use physics_types,    only: physics_state, physics_ptend
    use constituents,     only: cnst_get_ind
    use physics_buffer,   only: pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
    use ppgrid,           only: pcols, pver

    implicit none
    private

    public set_cdnc

subroutine set_cdnc(state, ptend, pbuf, ocnfrac, nc_reg1, nc_reg2, nc_reg3)

   use constituents,     only: cnst_get_ind
   use physics_buffer,   only: pbuf_old_tim_idx, pbuf_get_index
   use ppgrid,           only: pcols, pver

   type(physics_buffer_desc), pointer  :: pbuf(:)
   type(physics_state), intent(inout)  :: state
   type(physics_ptend), intent(in)     :: ptend
   real(r8),            intent(in)     :: ocnfrac(pcols)    ! Land fraction
   real(r8),            intent(in)     :: nc_reg1 ! CDNC target for region 1 (NEPac)
   real(r8),            intent(in)     :: nc_reg2 ! CDNC target for region 2 (SEPac)
   real(r8),            intent(in)     :: nc_reg3 ! CDNC target for region 3 (SEAtl)

   integer             :: ixnumliq
   integer             :: itim_old
   integer             :: alst_idx = -1
   integer             :: i, k, ncol
   real(r8), pointer, dimension(:,:) :: alst        ! cloud fraction

   ncol  = state%ncol

   itim_old = pbuf_old_tim_idx()
   alst_idx  = pbuf_get_index('ALST')
   call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)

   do i = 1, ncol  
      !! North East Pacific
      if(state%lat(i) > 0._r8 .and. state%lat(i) < 0.523599_r8 .and. &
         state%lon(i) > 3.665191_r8 .and. state%lon(i) < 4.363323_r8) then
            write(iulog,*) 'testing region 1 def lat=',state%lat(i), ', lon=',state%lon(i)
            do k = ptend%top_level, ptend%bot_level
               state%q(i,k,ixnumliq) = nc_reg1 * alst(i,k) * ocnfrac(i)
            end do
      end if

      !! South East Pacific
      if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. &
         state%lon(i) < 5.061455_r8 .and. state%lon(i) > 4.363323_r8) then
            write(iulog,*) 'testing region 2 def lat=',state%lat(i), ', lon=',state%lon(i)
            do k = ptend%top_level, ptend%bot_level
               state%q(i,k,ixnumliq) = nc_reg2 * alst(i,k) * ocnfrac(i)
            end do
      end if

      !! South East Atlantic
      if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. &
         state%lon(i) > 5.846853_r8) then
            write(iulog,*) 'testing region 3 def lat=',state%lat(i), ', lon=',state%lon(i)
            do k = ptend%top_level, ptend%bot_level
               state%q(i,k,ixnumliq) = nc_reg3 * alst(i,k)
            end do
      end if

      if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. &
         state%lon(i) < 0.261799_r8) then
            write(iulog,*) 'testing region 3 def lat=',state%lat(i), ', lon=',state%lon(i)
            do k = ptend%top_level, ptend%bot_level
               state%q(i,k,ixnumliq) = nc_reg3 * alst(i,k) * ocnfrac(i)
            end do
      end if

   end do

end subroutine set_cdnc


end module fix_cdnc