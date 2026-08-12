! This file contains miscellaneous subroutines needed by METISSE
! to work in stand alone or otherwise
! Ideally these should be packed in a module
! But Fortran 77 does not know how to use modules
! So here we are

subroutine initialize_front_end(front_end_name)
    use track_support
    character (len=*), intent(in) :: front_end_name

    if (verbose)print*, 'Setting front end to ',trim(front_end_name)

    if (ANY((/'TEST','test'/)== trim(front_end_name))) then
        ! for running unit tests
        front_end = test
    elseif (ANY((/'MAIN','main'/)== trim(front_end_name))) then
        ! METISSE's main code as described in Agrawal et al. 2020
        ! Can be used to evolve single stars and/or debugging purposes.
        front_end = main

    elseif (ANY((/'SSE','sse','BSE','bse'/)== trim(front_end_name))) then
        ! SSE (Single Star Evolution) from Hurley et al. 2000
        ! BSE (Binary Star Evolution) from Hurley et al. 2002
        front_end = BSE
        
    elseif (ANY((/'COSMIC','cosmic'/)== trim(front_end_name))) then
        ! COSMIC (Compact Object Synthesis and Monte Carlo Investigation Code)
        ! Binary evolution code from Breivik et al. 2020
        front_end = COSMIC
    elseif (ANY((/'AMUSE','amuse'/)== trim(front_end_name))) then
        ! Astrophysical Multipurpose Software Environment
        ! Pelupessy et al. 2010
        front_end = AMUSE
        
    else
        print*, "METISSE error: Unrecongnized front_end_name for METISSE"
        print*, "Choose from 'MAIN', 'SSE', 'BSE', 'COSMIC' "
    endif
    
end subroutine initialize_front_end

subroutine set_file_mode(i)
    use track_support, only: mode
    integer, intent(in) :: i
    mode = i
end subroutine

subroutine set_cmc_windowed_interp(flag)
    use track_support, only: cmc_windowed_interp
    logical, intent(in) :: flag
    cmc_windowed_interp = flag
end subroutine

subroutine allocate_track(n,mass)
    use track_support
    implicit none

    integer, intent(in):: n
    real(dp), intent(in), optional :: mass(:)
        
!    print*,"I am in alloc_track with ", n,mass

    allocate(tarr(n))
    tarr% star_type = unknown
    tarr% pars% age = 0.d0
    tarr% pars% extra = 0
    tarr% pars% bhspin = 0.d0
    tarr% ierr = 0
    tarr% pars% dms = 0.d0
    tarr% pars% delta = 0.d0
    tarr% reju = .false.
    
end subroutine allocate_track


subroutine dealloc_track()
    use track_support
    implicit none
    integer:: n,i

    n = size(tarr)

    do i = 1,n
        if (allocated(tarr(i)% eep)) then
            deallocate(tarr(i)% eep)
            deallocate(tarr(i)% tr)
            deallocate(tarr(i)% cols)
            deallocate(tarr(i)% bounds)
            if((tarr(i)% ierr/=0).and.verbose) write(UNIT=err_unit,fmt=*)'METISSE: error in evolving the system',i
        endif
        code_error = .false.
    end do
    
    deallocate(tarr)
end subroutine dealloc_track


subroutine set_star_type(id)

! set star type to rejuvenated before calling star
! id is the binary-component index (1 or 2), called directly from
! evolv2.f/comenv.f/mix.f (bypassing star.f/hrdiag.f); remap it to the
! persistent track-pool slot exactly as star.f/hrdiag.f do, so this
! doesn't corrupt (or fall short of) tarr under using_cmc
    use track_support
        implicit none
        integer, intent(in) :: id
        integer :: metisse_id

        real(dp) :: merger
        integer(8) :: id1_pass, id2_pass, using_cmc
        COMMON /CMCPASS/ merger, id1_pass, id2_pass, using_cmc
        integer :: track_id(2)
        COMMON /TRACKIDMAP/ track_id

!        print*, 'setting star to reju',tarr(id)% pars% age,id

        metisse_id = id
        if (using_cmc == 1) metisse_id = track_id(id)

        tarr(metisse_id)% star_type = rejuvenated
        tarr(metisse_id)% reju = .true.
end subroutine set_star_type


! ---------------------------------------------------------------------
! Persistent per-star track pool (CMC only). Lets tarr persist across
! repeated evolv2 calls for the same star, keyed by a CMC-supplied
! external id (permanent, globally unique, never reused). See the
! plan/design note in track_support.f90 next to the pool_* declarations.
! ---------------------------------------------------------------------

subroutine init_track_slot(slot)
    ! resets one tarr slot exactly as allocate_track resets a fresh array
    use track_support
    implicit none
    integer, intent(in) :: slot

    tarr(slot)% star_type = unknown
    tarr(slot)% pars% age = 0.d0
    tarr(slot)% pars% extra = 0
    tarr(slot)% pars% bhspin = 0.d0
    tarr(slot)% ierr = 0
    tarr(slot)% pars% dms = 0.d0
    tarr(slot)% pars% delta = 0.d0
    tarr(slot)% reju = .false.
end subroutine init_track_slot

subroutine pool_release_slot_arrays(slot)
    ! per-slot equivalent of dealloc_track's per-element cleanup. Does not
    ! touch the shared code_error flag (dealloc_track's reset of it makes
    ! sense when tearing down every star at once; here it would incorrectly
    ! clear a flag that may belong to a different star's error this timestep)
    use track_support
    implicit none
    integer, intent(in) :: slot

    if (allocated(tarr(slot)% eep)) then
        if ((tarr(slot)% ierr/=0) .and. verbose) &
            write(UNIT=err_unit,fmt=*)'METISSE: error in evolving the system',slot
        deallocate(tarr(slot)% eep)
        deallocate(tarr(slot)% tr)
        deallocate(tarr(slot)% cols)
        deallocate(tarr(slot)% bounds)
    endif
end subroutine pool_release_slot_arrays

subroutine pool_grow_tarr(min_size)
    ! grows tarr (preserving existing slots/indices) to hold at least min_size
    use track_support
    implicit none
    integer, intent(in) :: min_size
    type(track), allocatable :: tmp(:)
    integer :: old_size, new_size

    old_size = 0
    if (allocated(tarr)) old_size = size(tarr)
    if (old_size >= min_size) return

    new_size = max(min_size, max(8, old_size*2))
    allocate(tmp(new_size))
    if (old_size > 0) tmp(1:old_size) = tarr(1:old_size)
    call move_alloc(tmp, tarr)
end subroutine pool_grow_tarr

subroutine pool_grow_free_slots()
    use track_support
    implicit none
    integer, allocatable :: tmp(:)
    integer :: old_size, new_size

    old_size = 0
    if (allocated(pool_free_slots)) old_size = size(pool_free_slots)
    new_size = max(16, old_size*2)
    allocate(tmp(new_size))
    if (old_size > 0) tmp(1:old_size) = pool_free_slots(1:old_size)
    call move_alloc(tmp, pool_free_slots)
end subroutine pool_grow_free_slots

integer function pool_hash_func(key, table_size) result(h)
    use track_support, only: dp
    implicit none
    integer(8), intent(in) :: key
    integer, intent(in) :: table_size
    integer(8) :: k

    k = key
    if (k < 0_8) k = -k
    k = ieor(k, ishft(k, -17))
    k = k * 2654435761_8
    k = ieor(k, ishft(k, -13))
    h = int(iand(k, int(table_size-1, 8))) + 1
end function pool_hash_func

subroutine pool_hash_init(initial_size)
    use track_support
    implicit none
    integer, intent(in) :: initial_size

    if (allocated(pool_hash_key)) deallocate(pool_hash_key, pool_hash_slot)
    allocate(pool_hash_key(initial_size))
    allocate(pool_hash_slot(initial_size))
    pool_hash_key = pool_hash_empty
    pool_hash_count = 0
end subroutine pool_hash_init

subroutine pool_hash_find(key, idx, found)
    ! returns the slot in the hash table holding key (found=.true.), or
    ! the first empty/tombstone slot where key could be inserted (found=.false.)
    use track_support
    implicit none
    integer(8), intent(in) :: key
    integer, intent(out) :: idx
    logical, intent(out) :: found
    integer :: pool_hash_func
    integer :: h, i, first_tomb

    found = .false.
    first_tomb = -1
    h = pool_hash_func(key, size(pool_hash_key))
    i = h
    do
        if (pool_hash_key(i) == pool_hash_empty) then
            if (first_tomb > 0) then
                idx = first_tomb
            else
                idx = i
            endif
            return
        elseif (pool_hash_key(i) == pool_hash_tomb) then
            if (first_tomb < 0) first_tomb = i
        elseif (pool_hash_key(i) == key) then
            found = .true.
            idx = i
            return
        endif
        i = i + 1
        if (i > size(pool_hash_key)) i = 1
        if (i == h) then
            idx = -1
            return
        endif
    end do
end subroutine pool_hash_find

subroutine pool_hash_grow()
    use track_support
    implicit none
    integer(8), allocatable :: old_key(:)
    integer, allocatable :: old_slot(:)
    integer :: old_size, new_size, i, idx
    logical :: found

    old_size = size(pool_hash_key)
    old_key = pool_hash_key
    old_slot = pool_hash_slot
    new_size = old_size * 2

    call pool_hash_init(new_size)
    do i = 1, old_size
        if (old_key(i) >= 0_8) then
            call pool_hash_find(old_key(i), idx, found)
            pool_hash_key(idx) = old_key(i)
            pool_hash_slot(idx) = old_slot(i)
            pool_hash_count = pool_hash_count + 1
        endif
    end do
end subroutine pool_hash_grow

subroutine get_track_slot(id, slot)
    ! returns the persistent tarr slot for external id, creating one (with a
    ! freshly reset state, same as allocate_track gives a new star) if this
    ! id hasn't been seen before or was previously released. id==0 is CMC's
    ! "no second star" sentinel and always maps to one shared inert slot.
    use track_support
    implicit none
    integer(8), intent(in) :: id
    integer, intent(out) :: slot
    integer :: idx
    logical :: found

    if (id == 0_8) then
        if (pool_zero_slot == 0) then
            call pool_grow_tarr(pool_next_slot+1)
            pool_next_slot = pool_next_slot+1
            pool_zero_slot = pool_next_slot
            call init_track_slot(pool_zero_slot)
        endif
        slot = pool_zero_slot
        return
    endif

    if (.not. allocated(pool_hash_key)) call pool_hash_init(64)

    call pool_hash_find(id, idx, found)
    if (found) then
        slot = pool_hash_slot(idx)
        return
    endif

    ! keep load factor <= 0.6; re-probe in the grown table afterward
    if ((pool_hash_count+1) * 10 > size(pool_hash_key) * 6) then
        call pool_hash_grow()
        call pool_hash_find(id, idx, found)
    endif

    if (pool_n_free > 0) then
        slot = pool_free_slots(pool_n_free)
        pool_n_free = pool_n_free - 1
    else
        call pool_grow_tarr(pool_next_slot+1)
        pool_next_slot = pool_next_slot+1
        slot = pool_next_slot
    endif

    pool_hash_key(idx) = id
    pool_hash_slot(idx) = slot
    pool_hash_count = pool_hash_count + 1

    call init_track_slot(slot)
end subroutine get_track_slot

subroutine release_track_slot(id)
    ! frees the persistent tarr slot owned by external id, if any. Called by
    ! CMC when an id will never be evolved again (star ejected, merged away,
    ! etc). No-op for id==0 or an id with no current slot.
    use track_support
    implicit none
    integer(8), intent(in) :: id
    integer :: idx
    logical :: found

    if (id == 0_8) return
    if (.not. allocated(pool_hash_key)) return

    call pool_hash_find(id, idx, found)
    if (.not. found) return

    call pool_release_slot_arrays(pool_hash_slot(idx))

    if (.not. allocated(pool_free_slots)) call pool_grow_free_slots()
    if (pool_n_free >= size(pool_free_slots)) call pool_grow_free_slots()
    pool_n_free = pool_n_free + 1
    pool_free_slots(pool_n_free) = pool_hash_slot(idx)

    pool_hash_key(idx) = pool_hash_tomb
    pool_hash_count = pool_hash_count - 1
end subroutine release_track_slot

