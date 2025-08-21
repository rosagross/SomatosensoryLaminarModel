module system_equations

double precision :: PI = 4.0*atan(1.0)
complex :: I = (0.0, 1.0)

contains


subroutine vector_field(t,y,dy,tau,m_it,c_it,H,tau_v1,H_v1,tau_v2,H_v2,&
     & tau_v3,m_it_v1,c_it_v1,H_v3,tau_v4,H_v4,tau_v5,m_it_v2,c_it_v2,&
     & H_v5,tau_v6,H_v6,V_thr,r,m_max,weight_v2,V_thr_v1,r_v1,m_max_v1,&
     & weight_v3,V_thr_v2,r_v2,m_max_v2,weight,weight_v1)

implicit none

double precision, intent(in) :: t
double precision, intent(in) :: y(14)
double precision :: v
double precision :: i
double precision :: v_v1
double precision :: i_v1
double precision :: v_v2
double precision :: i_v2
double precision :: v_v3
double precision :: i_v3
double precision :: v_v4
double precision :: i_v4
double precision :: v_v5
double precision :: i_v5
double precision :: v_v6
double precision :: i_v6
double precision :: m_out_v2
double precision :: m_in_v2
double precision :: m_out
double precision :: m_in_v3
double precision :: m_out_v1
double precision :: m_in
double precision :: m_in_v1
double precision, intent(inout) :: dy(14)
double precision, intent(in) :: tau
double precision, intent(in) :: m_it
double precision, intent(in) :: c_it
double precision, intent(in) :: H
double precision, intent(in) :: tau_v1
double precision, intent(in) :: H_v1
double precision, intent(in) :: tau_v2
double precision, intent(in) :: H_v2
double precision, intent(in) :: tau_v3
double precision, intent(in) :: m_it_v1
double precision, intent(in) :: c_it_v1
double precision, intent(in) :: H_v3
double precision, intent(in) :: tau_v4
double precision, intent(in) :: H_v4
double precision, intent(in) :: tau_v5
double precision, intent(in) :: m_it_v2
double precision, intent(in) :: c_it_v2
double precision, intent(in) :: H_v5
double precision, intent(in) :: tau_v6
double precision, intent(in) :: H_v6
double precision, intent(in) :: V_thr
double precision, intent(in) :: r
double precision, intent(in) :: m_max
double precision, intent(in) :: weight_v2
double precision, intent(in) :: V_thr_v1
double precision, intent(in) :: r_v1
double precision, intent(in) :: m_max_v1
double precision, intent(in) :: weight_v3
double precision, intent(in) :: V_thr_v2
double precision, intent(in) :: r_v2
double precision, intent(in) :: m_max_v2
double precision, intent(in) :: weight
double precision, intent(in) :: weight_v1


v = y(1)
i = y(2)
v_v1 = y(3)
i_v1 = y(4)
v_v2 = y(5)
i_v2 = y(6)
v_v3 = y(7)
i_v3 = y(8)
v_v4 = y(9)
i_v4 = y(10)
v_v5 = y(11)
i_v5 = y(12)
v_v6 = y(13)
i_v6 = y(14)
m_out_v2 = 2.0*m_max/(exp(r*(V_thr - v - v_v1 - v_v2)) + 1)
m_in_v2 = m_out_v2*weight_v2
m_out = 2.0*m_max_v1/(exp(r_v1*(V_thr_v1 - v_v3 - v_v4)) + 1)
m_in_v3 = m_out_v2*weight_v3
m_out_v1 = 2.0*m_max_v2/(exp(r_v2*(V_thr_v2 - v_v5 - v_v6)) + 1)
m_in = m_out*weight
m_in_v1 = m_out_v1*weight_v1

dy(1) = i
dy(2) = H*c_it*m_it/tau - 2*i/tau - v/tau**2
dy(3) = i_v1
dy(4) = H_v1*m_in/tau_v1 - 2*i_v1/tau_v1 - v_v1/tau_v1**2
dy(5) = i_v2
dy(6) = H_v2*m_in_v1/tau_v2 - 2*i_v2/tau_v2 - v_v2/tau_v2**2
dy(7) = i_v3
dy(8) = H_v3*c_it_v1*m_it_v1/tau_v3 - 2*i_v3/tau_v3 - v_v3/tau_v3**2
dy(9) = i_v4
dy(10) = H_v4*m_in_v2/tau_v4 - 2*i_v4/tau_v4 - v_v4/tau_v4**2
dy(11) = i_v5
dy(12) = H_v5*c_it_v2*m_it_v2/tau_v5 - 2*i_v5/tau_v5 - v_v5/tau_v5**2
dy(13) = i_v6
dy(14) = H_v6*m_in_v3/tau_v6 - 2*i_v6/tau_v6 - v_v6/tau_v6**2

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
     & args(23), args(24), args(25), args(26), args(27), args(28), &
     & args(29), args(30), args(31), args(32), args(33), args(34), &
     & args(35), args(36), args(37), args(38))

end subroutine func


subroutine stpnt(ndim, y, args, t)

implicit None
integer, intent(in) :: ndim
double precision, intent(inout) :: y(ndim), args(*)
double precision, intent(in) :: t

args(1) = 0.01  ! tau
args(2) = 0.0  ! m_it
args(3) = 1.0  ! c_it
args(4) = 0.00325  ! H
args(5) = 0.01  ! tau_v1
args(6) = 0.00325  ! H_v1
args(7) = 0.02  ! tau_v2
args(8) = -0.022  ! H_v2
args(9) = 0.01  ! tau_v3
args(15) = 0.0  ! m_it_v1
args(16) = 1.0  ! c_it_v1
args(17) = 0.00325  ! H_v3
args(18) = 0.01  ! tau_v4
args(19) = 0.00325  ! H_v4
args(20) = 0.01  ! tau_v5
args(21) = 0.0  ! m_it_v2
args(22) = 1.0  ! c_it_v2
args(23) = 0.00325  ! H_v5
args(24) = 0.01  ! tau_v6
args(25) = 0.00325  ! H_v6
args(26) = 0.006  ! V_thr
args(27) = 560.0  ! r
args(28) = 2.5  ! m_max
args(29) = 135.0  ! weight_v2
args(30) = 0.006  ! V_thr_v1
args(31) = 560.0  ! r_v1
args(32) = 2.5  ! m_max_v1
args(33) = 33.75  ! weight_v3
args(34) = 0.006  ! V_thr_v2
args(35) = 560.0  ! r_v2
args(36) = 2.5  ! m_max_v2
args(37) = 108.0  ! weight
args(38) = 33.75  ! weight_v1
y(1) = 0.0  ! v
y(2) = 0.0  ! i
y(3) = 0.0  ! v_v1
y(4) = 0.0  ! i_v1
y(5) = 0.0  ! v_v2
y(6) = 0.0  ! i_v2
y(7) = 0.0  ! v_v3
y(8) = 0.0  ! i_v3
y(9) = 0.0  ! v_v4
y(10) = 0.0  ! i_v4
y(11) = 0.0  ! v_v5
y(12) = 0.0  ! i_v5
y(13) = 0.0  ! v_v6
y(14) = 0.0  ! i_v6

end subroutine stpnt



subroutine bcnd
end subroutine bcnd


subroutine icnd
end subroutine icnd


subroutine fopt
end subroutine fopt


subroutine pvls
end subroutine pvls
