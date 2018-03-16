! Written in 2015 by Shun Sakuraba
! 
! To the extent possible under law, the author has dedicated all copyright
! and related and neighboring rights to this software to the public domain
! worldwide. This software is distributed without any warranty.
! See <http://creativecommons.org/publicdomain/zero/1.0/>.
! 
! Vigna's xorshift128+ pseudorandom generator.
! (Sebastiano Vigna. Further scramblings of Marsaglia's xorshift generators. CoRR, abs/1402.6246, 2014.)
! xorshift128+ is known as a fast pseudorandom generator with reasonably good resilience to randomness tests.
! Since its primary imporance is its speed, do not forget to add inlining directives depending your compiler.
! 
! TODO FIXME: initialization

module xorshift128plus
  implicit none
  
  ! random number state
  public :: xorshift128plus_state
  ! state initialization functions
  public :: state_init_full, state_init, state_init_self
  ! global state initialization functions
  public :: rand_init, rand_init_self

  ! Draw integer from range, or uniform random number
  public :: draw_integer, draw_uniform
  ! Draw integer from range, or uniform random number, from global rng state
  public :: rand_integer, rand_uniform

  type xorshift128plus_state
     integer(8) :: s0, s1
  end type xorshift128plus_state

  type(xorshift128plus_state) :: global_state
  data global_state%s0 / 451509857683038208_8 /, global_state%s1 / -7168371501937287167_8 /

contains

  subroutine state_init_full(state, s0, s1)
    implicit none
    type(xorshift128plus_state), intent(out) :: state
    integer(8), intent(in) :: s0, s1
    integer(8) :: s0_copy

    s0_copy = s0
    ! xorshift generators does not allow 0/0 state!
    if(s0_copy == 0_8 .and. s1 == 0_8) then
       s0_copy = -1_8
    end if
    
    state%s0 = s0_copy
    state%s1 = s1
  end subroutine state_init_full

  integer(8) function draw_integer8(state)
    implicit none
    type(xorshift128plus_state), intent(inout) :: state
    
    integer(8) :: s0_new, s1_new
    
    ! swap state variables
    s1_new = state%s0
    s0_new = state%s1
    state%s0 = s0_new

    ! xorshift
    s1_new = ieor(ishft(s1_new, 23), s1_new)

    ! scramble
    s1_new = ieor(ieor(s1_new, s0_new), ieor(ishft(s1_new, -17), ishft(s0_new, -23)))
    state%s1 = s1_new

    draw_integer8 = s1_new + s0_new
  end function draw_integer8
  
  ! Following are utility functions
  
  subroutine state_init(state, seed)
    implicit none
    type(xorshift128plus_state), intent(out) :: state
    integer, intent(in) :: seed
    integer(8) :: tmp
    integer(8) :: x0
    integer(8), parameter :: spreader = 2685821657736338717_8

    ! do something similar to xorshift64*
    tmp = int(seed, 8)
    tmp = tmp * spreader
    x0 = tmp
    x0 = ieor(x0, ishft(x0, -12))
    x0 = ieor(x0, ishft(x0, 25))
    x0 = ieor(x0, ishft(x0, -27))
    x0 = x0 * spreader

    call state_init_full(state, x0, tmp)
  end subroutine state_init

  subroutine rand_init(seed)
    implicit none
    integer, intent(in) :: seed
    call state_init(global_state, seed)
  end subroutine rand_init

  subroutine state_init_self(state)
    implicit none
    type(xorshift128plus_state), intent(out) :: state
    integer :: current_time(8)
    call date_and_time(VALUES=current_time)
    if (current_time(4) == -huge(0)) then
       stop "xorshift128plus:state_init_self: too exotic system! (UNIX time not available)"
    endif
    call state_init(state, current_time(4))
  end subroutine state_init_self

  subroutine rand_init_self()
    call state_init_self(global_state)
  end subroutine rand_init_self
  
  ! rmax must be positive
  integer function draw_integer(state, rmax)
    implicit none
    type(xorshift128plus_state), intent(inout) :: state
    integer, intent(in) :: rmax
    integer(8) :: rmax8
    integer(8), parameter :: rmask = 9223372036854775807_8 ! (1 << 63) - 1
    integer(8) :: rnd, qmax, q

    rmax8 = int(rmax, 8)
    ! real maximum is (rmask + 1) / rmax8, but it is impossible to do in 64-bit arithmetic.
    ! Instead we evaluate conservatively. Even if the compiler uses integer(8) as a default integer, 
    ! rmax is maxed at ((1 << 63) - 1), thus qmax >= 1 is ensured.
    qmax = rmask / rmax8
    do
       ! mask to convert to positive integer
       rnd = iand(draw_integer8(state), rmask)
       ! Now both are positive, divide to check whether it is ok to go
       q = rnd / rmax8
       if(q < qmax) then
          draw_integer = int(mod(rnd, rmax8), kind=kind(rmax))
          exit
       endif
       ! otherwise repeat the last step to ensure equidistribution
    end do
  end function draw_integer

  integer(8) function rand_integer(rmax)
    implicit none
    integer, intent(in) :: rmax
    rand_integer = draw_integer(global_state, rmax)
  end function rand_integer

  real(8) function draw_uniform(state)
    implicit none
    type(xorshift128plus_state), intent(inout) :: state
    integer(8) :: rnd 

    ! 1.0 / (1 << 53)
    real(8), parameter :: multiplier = 1.0d0 / 9007199254740992d0

    rnd = draw_integer8(state)
    
    ! 53-bit, divided by 2^53
    draw_uniform = real(ishft(rnd, -11), kind=8) * multiplier
  end function draw_uniform

  real(8) function rand_uniform()
    implicit none
    rand_uniform = draw_uniform(global_state)
  end function rand_uniform

end module xorshift128plus
