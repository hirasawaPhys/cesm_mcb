module fix_cdnc

    use physics_types,    only: physics_state, physics_ptend
    use constituents,     only: cnst_get_ind

    public set_cdnc

    subroutine set_cdnc(state, ptend, nc_reg1, nc_reg2, nc_reg3)
        type(physics_ptend) :: ptend
        type(physics_state) :: state
        real(r8)            :: nc_reg1
        real(r8)            :: nc_reg2
        real(r8)            :: nc_reg3

        call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)
        do i = 1, ncol  
            ! North East Pacific
            if(state%lat(i) > 0._r8 .and. state%lat(i) < 0.523599_r8 .and. state%lon(i) > 3.665191_r8 .and. state%lon(i) < 4.363323_r8) then
                write(iulog,*) 'testing region 1 def lat=',state%lat(i), ', lon=',state%lon(i)
                do k = ptend%top_level, ptend%bot_level
                    state%q(i,k,ixnumliq) = nc_reg1 * cld(i,k)
                end do
            end if
            ! South East Pacific
            if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. state%lon(i) < 5.061455_r8 .and. state%lon(i) > 4.363323_r8) then
                write(iulog,*) 'testing region 2 def lat=',state%lat(i), ', lon=',state%lon(i)
                do k = ptend%top_level, ptend%bot_level
                    state%q(i,k,ixnumliq) = nc_reg2 * cld(i,k)
                end do
            end if
            ! South East Atlantic
            if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. state%lon(i) > 5.846853_r8) then
                write(iulog,*) 'testing region 3 def lat=',state%lat(i), ', lon=',state%lon(i)
                do k = ptend%top_level, ptend%bot_level
                    state%q(i,k,ixnumliq) = nc_reg3 * cld(i,k)
                end do
            end if
            if(state%lat(i) < 0._r8 .and. state%lat(i) > -0.523599_r8 .and. state%lon(i) < 0.261799_r8) then
                write(iulog,*) 'testing region 3 def lat=',state%lat(i), ', lon=',state%lon(i)
                do k = ptend%top_level, ptend%bot_level
                    state%q(i,k,ixnumliq) = nc_reg3 * cld(i,k)
                end do
            end if

        end do
    end subroutine set_cdnc