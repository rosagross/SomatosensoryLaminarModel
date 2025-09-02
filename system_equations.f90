module system_equations

double precision :: PI = 4.0*atan(1.0)
complex :: I = (0.0, 1.0)

contains


subroutine vector_field(t,y,dy,u,tau,h,tau_v1,h_v1,tau_v2,h_v2,tau_v3,&
     & h_v3,v_thr,s,m_max,weight_v2,v_thr_v1,s_v1,m_max_v1,weight_v3,&
     & v_thr_v2,s_v2,m_max_v2,weight,weight_v1)

implicit none

double precision, intent(in) :: t
double precision, intent(in) :: y(8)
double precision :: v
double precision :: x
double precision :: v_v1
double precision :: x_v1
double precision :: v_v2
double precision :: x_v2
double precision :: v_v3
double precision :: x_v3
double precision :: m_v2
double precision :: m_in_v2
double precision :: m
double precision :: m_in_v3
double precision :: m_v1
double precision :: m_in
double precision :: m_in_v1
double precision, intent(inout) :: dy(8)
double precision, intent(in) :: u
double precision, intent(in) :: tau
double precision, intent(in) :: h
double precision, intent(in) :: tau_v1
double precision, intent(in) :: h_v1
double precision, intent(in) :: tau_v2
double precision, intent(in) :: h_v2
double precision, intent(in) :: tau_v3
double precision, intent(in) :: h_v3
double precision, intent(in) :: v_thr
double precision, intent(in) :: s
double precision, intent(in) :: m_max
double precision, intent(in) :: weight_v2
double precision, intent(in) :: v_thr_v1
double precision, intent(in) :: s_v1
double precision, intent(in) :: m_max_v1
double precision, intent(in) :: weight_v3
double precision, intent(in) :: v_thr_v2
double precision, intent(in) :: s_v2
double precision, intent(in) :: m_max_v2
double precision, intent(in) :: weight
double precision, intent(in) :: weight_v1


v = y(1)
x = y(2)
v_v1 = y(3)
x_v1 = y(4)
v_v2 = y(5)
x_v2 = y(6)
v_v3 = y(7)
x_v3 = y(8)
m_v2 = m_max/(exp(s*(-v + v_thr - v_v1)) + 1)
m_in_v2 = m_v2*weight_v2
m = m_max_v1/(exp(s_v1*(v_thr_v1 - v_v2)) + 1)
m_in_v3 = m_v2*weight_v3
m_v1 = m_max_v2/(exp(s_v2*(v_thr_v2 - v_v3)) + 1)
m_in = m*weight
m_in_v1 = m_v1*weight_v1

dy(1) = x
dy(2) = h*(m_in + u)/tau - 2*x/tau - v/tau**2
dy(3) = x_v1
dy(4) = h_v1*m_in_v1/tau_v1 - 2*x_v1/tau_v1 - v_v1/tau_v1**2
dy(5) = x_v2
dy(6) = h_v2*m_in_v2/tau_v2 - 2*x_v2/tau_v2 - v_v2/tau_v2**2
dy(7) = x_v3
dy(8) = h_v3*m_in_v3/tau_v3 - 2*x_v3/tau_v3 - v_v3/tau_v3**2

end subroutine


end module


subroutine func(ndim,y,icp,args,ijac,dy,dfdu,dfdp)

use system_equations
implicit none
integer, intent(in) :: ndim, icp(*), ijac
double precision, intent(in) :: y(ndim), args(*)
double precision, intent(out) :: dy(ndim)
double precision, intent(inout) :: dfdu(ndim,ndim), dfdp(ndim,*)

call vector_field(args(14), y, dy, args(1), args(2), args(3), args(4), &
     & args(5), args(6), args(7), args(8), args(9), args(15), args(16),&
     &  args(17), args(18), args(19), args(20), args(21), args(22), &
     & args(23), args(24), args(25), args(26), args(27))

end subroutine func


subroutine stpnt(ndim, y, args, t)

implicit None
integer, intent(in) :: ndim
double precision, intent(inout) :: y(ndim), args(*)
double precision, intent(in) :: t

args(1) = 400.0  ! u
args(2) = 0.01  ! tau
args(3) = 0.00325  ! h
args(4) = 0.02  ! tau_v1
args(5) = -0.022  ! h_v1
args(6) = 0.01  ! tau_v2
args(7) = 0.00325  ! h_v2
args(8) = 0.01  ! tau_v3
args(9) = 0.00325  ! h_v3
args(15) = 0.006  ! v_thr
args(16) = 560.0  ! s
args(17) = 5.0  ! m_max
args(18) = 135.0  ! weight_v2
args(19) = 0.006  ! v_thr_v1
args(20) = 560.0  ! s_v1
args(21) = 5.0  ! m_max_v1
args(22) = 33.75  ! weight_v3
args(23) = 0.006  ! v_thr_v2
args(24) = 560.0  ! s_v2
args(25) = 5.0  ! m_max_v2
args(26) = 108.0  ! weight
args(27) = 33.75  ! weight_v1
y(1) = 0.03052549  ! v
y(2) = 0.0  ! x
y(3) = -0.02195986  ! v_v1
y(4) = 0.0  ! x_v1
y(5) = 0.01762243  ! v_v2
y(6) = 0.0  ! x_v2
y(7) = 0.004405607  ! v_v3
y(8) = 0.0  ! x_v3

end subroutine stpnt



subroutine bcnd
end subroutine bcnd


subroutine icnd
end subroutine icnd


subroutine fopt
end subroutine fopt


subroutine pvls
end subroutine pvls
