module utils
implicit none
private
public dp, pi, imagunit

! Double precision kind
integer, parameter :: dp = kind(0.d0)

! Pi
real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

! Imaginary unit
complex(dp), parameter :: imagunit = (0.0_dp, 1.0_dp)

end module utils