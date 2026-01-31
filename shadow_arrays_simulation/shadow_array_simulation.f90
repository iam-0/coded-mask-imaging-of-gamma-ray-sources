program shadow_pattern_simulation

  implicit none
  integer, parameter :: N = 16          ! length of the masked grid
  real, parameter    :: side = 40       ! total side length (mm)
  real, parameter    :: H = 600         ! distance from mask to detector (mm)
  real, parameter    :: pixel = side/N  ! size of each pixel/square
  real, parameter    :: pi = acos(-1.0)
  real, parameter    :: deg = 180/pi    ! for radian to degree conversion 

  integer :: mask(N, N), shadow(N, N)
  integer :: jx, jy, r, c, ios, i !index no. for each shadow pattern 
  real    :: dx, dy, theta_x, theta_y

  print*, "Side length of each pixel is ",pixel,"mm."

! -------------------------------
! 1. Read the mask pattern
! -------------------------------
  call read_mask(mask, N, ios)
  if (ios /= 0) then
     print *, "Error: could not read mask pattern file. ios = ", ios
     stop
  end if

  print *, "Mask pattern (1=blocked, 0=open):"
  do r = 1, N
     print *, (mask(r, c), c = 1, N)
  end do

! -------------------------------
! 2. Generate shadow patterns
! -------------------------------
  open(unit=20, action='write', file='shadow_arrays.txt', status='replace', iostat=ios)   
  if (ios /=0) then
          print*, "Error in opening the shadow arrays text file ", ios
          stop
  end if

  write(20, '(A)') "# Shadow pattern simulation"
  write(20, '(A)') "# Shadow array (1=light, 0=shadowed)"
  write(20, '(A)') "# S_i (jx,jy) (dx,dy) (theta_x,theta_y)"

  do jy = -15, 15           !jy goes from top to bottom as its value increases
     do jx = -15, 15        !jx goes from left to right as its value increases
        i = jx + 31*jy + 480

        ! Calculate shadow displacement and angles
        dx = jx * pixel
        dy = -jy * pixel
        theta_x = atan(dx / H) * deg
        theta_y = atan(dy / H) * deg

        ! Compute shadow pattern
        call make_shadow(mask, shadow, N, jx, jy)

        ! Write results
        write(20, '(A,I3.3, 3X,"(",I3,",",I3,")", 3X,"(",F10.4,",",F10.4,")", 3X,"(",F10.6,",",F10.6,")")')&
                "S_", i, jx, jy, dx, dy, theta_x, theta_y
        write(20, '(256I2)') ((shadow(r,c), c=1,N), r=1,N) 

        
     end do
  end do

  close(20)

  print *, "Simulation complete. Total patterns generated:", i + 1

contains
! ---------------------------------------------------------------
! Subroutine: read_mask
! Reads the mask pattern from file "asymmetric_mask_pattern.txt"
! ---------------------------------------------------------------
  subroutine read_mask(mask, N, ios)
    implicit none
    integer, intent(out) :: mask(N, N)
    integer, intent(in)  :: N
    integer, intent(out) :: ios
    integer :: r, c
    ios = 0

    open(unit=10, file='asymmetric_mask_pattern.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) return

    do r = 1, N
       read(10, *, iostat=ios) (mask(r, c), c = 1, N)
       if (ios /= 0) exit
    end do

    close(10)
  end subroutine read_mask

! -------------------------------------------------------------------
! Subroutine: make_shadow
! Creates the shifted shadow pattern based on jx, jy
! -------------------------------------------------------------------
  subroutine make_shadow(mask, shadow, N, jx, jy)
    implicit none
    integer, intent(in) :: N, jx, jy
    integer, intent(in) :: mask(N, N)
    integer, intent(out) :: shadow(N, N)
    integer :: r, c

    ! Initialize all pixels, count to zero (0)
    shadow = 0

    do r = 1, N
       do c = 1, N
          if (mask(r, c) == 0) then
             if (r+jy >= 1 .and. r+jy <= N .and. c+jx >= 1 .and. c+jx <= N) then
                shadow(r+jy, c+jx) = 1
             end if
          end if
       end do
    end do
  end subroutine make_shadow

end program shadow_pattern_simulation

