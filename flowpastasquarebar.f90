module my_precision
  implicit none
  ! 有効桁数から種別番号を求めるためにselected_real_kindを用いる
  ! integer,parameter :: i_r = selected_real_kind(p=6)   !　p=6  : single precision
  integer,parameter :: i_r = selected_real_kind(p=12)  !　p=12 : double precision
  ! integer,parameter :: i_r = selected_real_kind(p=24)  !　p=24 : quadruple precision
  ! integer,parameter :: i_r = selected_real_kind(p=18)  !　p=18 : real(10)
end module my_precision

program potential_flow
  use my_precision
  implicit none
  ! two-dimensional potential flow described by complex velocity potential
  !   W   = phi + i . psi = U0 . z + i . zeta_j . log ( z - z_j )
  ! dW/dz = u - i . v = U0 + i . zeta_j / ( z - z_j )
  ! phi : velocity potential
  ! psi : stream function
  !  u  : x component of velocity vector
  !  v  : y component of velocity vector
  real(i_r) :: U0 ! uniform velocity
  real(i_r) :: a ! length of square
  real(i_r),dimension(:),allocatable :: zeta ! strength of vortex filament
  complex(i_r),dimension(:),allocatable :: zp ! positions located vortex filaments
  complex(i_r),dimension(:),allocatable :: bp ! positions where b. c. are evaluated
  complex(i_r),dimension(:,:),allocatable :: mat_a
  complex(i_r),dimension(:),allocatable :: vec_x, vec_y
  integer :: np ! number of points having vortex filament
  integer :: i,j,n

  U0 = 1.0_i_r
  a = 1.0_i_r
  np = 4

  allocate ( zeta (1:np) ) ; zeta = 0.0_i_r
  allocate ( zp (1:np) ) ; zp = 0.0_i_r
  allocate ( bp (1:np) ) ; bp = 0.0_i_r
  allocate ( vec_x (1:np) ) ; vec_x = 0.0_i_r
  allocate ( vec_y (1:np) ) ; vec_y = 0.0_i_r
  allocate ( mat_a (1:np,1:np) ) ; mat_a = 0.0_i_r

  zp (1) = cmplx ( -a/2., -a/2., i_r )
  zp (2) = cmplx ( +a/2., -a/2., i_r )
  zp (3) = cmplx ( +a/2., +a/2., i_r )
  zp (4) = cmplx ( -a/2., +a/2., i_r )

  bp (1:np-1) = ( zp(1:np-1) + zp(2:np) ) / 2.
  bp (np) = ( zp(np) + zp(1) ) / 2.

  do i = 1, np
    write (*,*) i, zp ( i )
  end do
  do i = 1, np
    write (*,*) i, bp ( i )
  end do
  
  stop
end program potential_flow

