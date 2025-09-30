module system_equations

double precision :: PI = 4.0*atan(1.0)
complex :: I = (0.0, 1.0)

contains


subroutine vector_field(t,y,dy,tau,H,tau_v1,H_v1,tau_v2,H_v2,tau_v3,&
     & H_v3,tau_v4,H_v4,tau_v5,H_v5,tau_v6,H_v6,tau_v7,H_v7,tau_v8,&
     & H_v8,tau_v9,H_v9,tau_v10,H_v10,tau_v11,H_v11,tau_v12,H_v12,&
     & tau_v13,H_v13,tau_v14,H_v14,tau_v15,H_v15,tau_v16,H_v16,tau_v17,&
     & H_v17,tau_v18,H_v18,tau_v19,H_v19,tau_v20,H_v20,tau_v21,H_v21,&
     & tau_v22,H_v22,tau_v23,H_v23,tau_v24,H_v24,tau_v25,H_v25,tau_v26,&
     & H_v26,tau_v27,H_v27,tau_v28,H_v28,tau_v29,H_v29,tau_v30,H_v30,&
     & tau_v31,H_v31,tau_v32,H_v32,tau_v33,H_v33,tau_v34,H_v34,tau_v35,&
     & H_v35,tau_v36,H_v36,tau_v37,bI,H_v37,V_thr,r,m_max,V_thr_v1,&
     & v_v12,r_v1,m_max_v1,V_thr_v2,v_v19,r_v2,m_max_v2,V_thr_v3,v_v26,&
     & r_v3,m_max_v3,onset,dur,A,V_thr_v4,r_v4,m_max_v4,V_thr_v5,r_v5,&
     & m_max_v5,g_input,g_thal_input,bEI_input,bEI_thal_input,&
     & connect_reverse_factor,weight_v6,connect_reverse_factor_v1,&
     & weight_v12,weight_v13,connect_reverse_factor_v2,weight_v18,&
     & weight_v19,weight_v20,connect_reverse_factor_v3,weight_v24,&
     & weight_v25,weight_v26,weight_v27,connect_reverse_factor_thal,&
     & weight_v30,weight_v31,weight_v32,weight_v33,weight_v34,&
     & connect_reverse_factor_thal_v1,weight,weight_v1,weight_v2,&
     & weight_v3,weight_v4,weight_v5,weight_v7,weight_v8,weight_v9,&
     & weight_v10,weight_v11,weight_v14,weight_v15,weight_v16,&
     & weight_v17,weight_v21,weight_v22,weight_v23,weight_v28,&
     & weight_v29,weight_v35)

implicit none

double precision, intent(in) :: t
double precision, intent(in) :: y(76)
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
double precision :: v_v7
double precision :: i_v7
double precision :: v_v8
double precision :: i_v8
double precision :: v_v9
double precision :: i_v9
double precision :: v_v10
double precision :: i_v10
double precision :: v_v11
double precision :: i_v11
double precision :: v_v13
double precision :: i_v12
double precision :: v_v14
double precision :: i_v13
double precision :: v_v15
double precision :: i_v14
double precision :: v_v16
double precision :: i_v15
double precision :: v_v17
double precision :: i_v16
double precision :: v_v18
double precision :: i_v17
double precision :: v_v20
double precision :: i_v18
double precision :: v_v21
double precision :: i_v19
double precision :: v_v22
double precision :: i_v20
double precision :: v_v23
double precision :: i_v21
double precision :: v_v24
double precision :: i_v22
double precision :: v_v25
double precision :: i_v23
double precision :: v_v27
double precision :: i_v24
double precision :: v_v28
double precision :: i_v25
double precision :: v_v29
double precision :: i_v26
double precision :: v_v30
double precision :: i_v27
double precision :: v_v31
double precision :: i_v28
double precision :: v_v32
double precision :: i_v29
double precision :: v_v33
double precision :: i_v30
double precision :: v_v34
double precision :: i_v31
double precision :: v_v35
double precision :: i_v32
double precision :: v_v36
double precision :: i_v33
double precision :: v_v37
double precision :: i_v34
double precision :: v_v38
double precision :: i_v35
double precision :: v_v39
double precision :: i_v36
double precision :: v_bI
double precision :: i_v37
double precision :: v_bIn
double precision :: m_outC
double precision :: v_bIn_v1
double precision :: m_outC_v1
double precision :: v_bIn_v2
double precision :: m_outC_v2
double precision :: v_bIn_v3
double precision :: m_outC_v3
double precision :: Iext
double precision :: m_outC_v4
double precision :: m_outC_v5
double precision :: gC
double precision :: g_thalC
double precision :: bEIC
double precision :: bEI_thalC
double precision :: g
double precision :: bEI
double precision :: m_out
double precision :: g_v1
double precision :: bEI_v1
double precision :: m_in_v6
double precision :: m_out_v1
double precision :: g_v2
double precision :: bEI_v2
double precision :: m_in_v12
double precision :: m_in_v13
double precision :: m_out_v2
double precision :: g_v3
double precision :: bEI_v3
double precision :: m_in_v18
double precision :: m_in_v19
double precision :: m_in_v20
double precision :: m_out_v3
double precision :: g_thal
double precision :: bEI_thal
double precision :: m_in_v24
double precision :: m_in_v25
double precision :: m_in_v26
double precision :: m_in_v27
double precision :: m_out_v4
double precision :: g_thal_v1
double precision :: bEI_thal_v1
double precision :: m_in_v30
double precision :: m_in_v31
double precision :: m_in_v32
double precision :: m_in_v33
double precision :: m_in_v34
double precision :: m_out_v5
double precision :: m_in
double precision :: m_in_v1
double precision :: m_in_v2
double precision :: m_in_v3
double precision :: m_in_v4
double precision :: m_in_v5
double precision :: m_in_v7
double precision :: m_in_v8
double precision :: m_in_v9
double precision :: m_in_v10
double precision :: m_in_v11
double precision :: m_in_v14
double precision :: m_in_v15
double precision :: m_in_v16
double precision :: m_in_v17
double precision :: m_in_v21
double precision :: m_in_v22
double precision :: m_in_v23
double precision :: m_in_v28
double precision :: m_in_v29
double precision :: m_in_v35
double precision, intent(inout) :: dy(76)
double precision, intent(in) :: tau
double precision, intent(in) :: H
double precision, intent(in) :: tau_v1
double precision, intent(in) :: H_v1
double precision, intent(in) :: tau_v2
double precision, intent(in) :: H_v2
double precision, intent(in) :: tau_v3
double precision, intent(in) :: H_v3
double precision, intent(in) :: tau_v4
double precision, intent(in) :: H_v4
double precision, intent(in) :: tau_v5
double precision, intent(in) :: H_v5
double precision, intent(in) :: tau_v6
double precision, intent(in) :: H_v6
double precision, intent(in) :: tau_v7
double precision, intent(in) :: H_v7
double precision, intent(in) :: tau_v8
double precision, intent(in) :: H_v8
double precision, intent(in) :: tau_v9
double precision, intent(in) :: H_v9
double precision, intent(in) :: tau_v10
double precision, intent(in) :: H_v10
double precision, intent(in) :: tau_v11
double precision, intent(in) :: H_v11
double precision, intent(in) :: tau_v12
double precision, intent(in) :: H_v12
double precision, intent(in) :: tau_v13
double precision, intent(in) :: H_v13
double precision, intent(in) :: tau_v14
double precision, intent(in) :: H_v14
double precision, intent(in) :: tau_v15
double precision, intent(in) :: H_v15
double precision, intent(in) :: tau_v16
double precision, intent(in) :: H_v16
double precision, intent(in) :: tau_v17
double precision, intent(in) :: H_v17
double precision, intent(in) :: tau_v18
double precision, intent(in) :: H_v18
double precision, intent(in) :: tau_v19
double precision, intent(in) :: H_v19
double precision, intent(in) :: tau_v20
double precision, intent(in) :: H_v20
double precision, intent(in) :: tau_v21
double precision, intent(in) :: H_v21
double precision, intent(in) :: tau_v22
double precision, intent(in) :: H_v22
double precision, intent(in) :: tau_v23
double precision, intent(in) :: H_v23
double precision, intent(in) :: tau_v24
double precision, intent(in) :: H_v24
double precision, intent(in) :: tau_v25
double precision, intent(in) :: H_v25
double precision, intent(in) :: tau_v26
double precision, intent(in) :: H_v26
double precision, intent(in) :: tau_v27
double precision, intent(in) :: H_v27
double precision, intent(in) :: tau_v28
double precision, intent(in) :: H_v28
double precision, intent(in) :: tau_v29
double precision, intent(in) :: H_v29
double precision, intent(in) :: tau_v30
double precision, intent(in) :: H_v30
double precision, intent(in) :: tau_v31
double precision, intent(in) :: H_v31
double precision, intent(in) :: tau_v32
double precision, intent(in) :: H_v32
double precision, intent(in) :: tau_v33
double precision, intent(in) :: H_v33
double precision, intent(in) :: tau_v34
double precision, intent(in) :: H_v34
double precision, intent(in) :: tau_v35
double precision, intent(in) :: H_v35
double precision, intent(in) :: tau_v36
double precision, intent(in) :: H_v36
double precision, intent(in) :: tau_v37
double precision, intent(in) :: bI
double precision, intent(in) :: H_v37
double precision, intent(in) :: V_thr
double precision, intent(in) :: r
double precision, intent(in) :: m_max
double precision, intent(in) :: V_thr_v1
double precision, intent(in) :: v_v12
double precision, intent(in) :: r_v1
double precision, intent(in) :: m_max_v1
double precision, intent(in) :: V_thr_v2
double precision, intent(in) :: v_v19
double precision, intent(in) :: r_v2
double precision, intent(in) :: m_max_v2
double precision, intent(in) :: V_thr_v3
double precision, intent(in) :: v_v26
double precision, intent(in) :: r_v3
double precision, intent(in) :: m_max_v3
double precision, intent(in) :: onset
double precision, intent(in) :: dur
double precision, intent(in) :: A
double precision, intent(in) :: V_thr_v4
double precision, intent(in) :: r_v4
double precision, intent(in) :: m_max_v4
double precision, intent(in) :: V_thr_v5
double precision, intent(in) :: r_v5
double precision, intent(in) :: m_max_v5
double precision, intent(in) :: g_input
double precision, intent(in) :: g_thal_input
double precision, intent(in) :: bEI_input
double precision, intent(in) :: bEI_thal_input
integer, intent(in) :: connect_reverse_factor(1)
double precision, intent(in) :: weight_v6
integer, intent(in) :: connect_reverse_factor_v1(1)
double precision, intent(in) :: weight_v12
double precision, intent(in) :: weight_v13
integer, intent(in) :: connect_reverse_factor_v2(1)
double precision, intent(in) :: weight_v18
double precision, intent(in) :: weight_v19
double precision, intent(in) :: weight_v20
integer, intent(in) :: connect_reverse_factor_v3(1)
double precision, intent(in) :: weight_v24
double precision, intent(in) :: weight_v25
double precision, intent(in) :: weight_v26
double precision, intent(in) :: weight_v27
integer, intent(in) :: connect_reverse_factor_thal(1)
double precision, intent(in) :: weight_v30
double precision, intent(in) :: weight_v31
double precision, intent(in) :: weight_v32
double precision, intent(in) :: weight_v33
double precision, intent(in) :: weight_v34
integer, intent(in) :: connect_reverse_factor_thal_v1(1)
double precision, intent(in) :: weight
double precision, intent(in) :: weight_v1
double precision, intent(in) :: weight_v2
double precision, intent(in) :: weight_v3
double precision, intent(in) :: weight_v4
double precision, intent(in) :: weight_v5
double precision, intent(in) :: weight_v7
double precision, intent(in) :: weight_v8
double precision, intent(in) :: weight_v9
double precision, intent(in) :: weight_v10
double precision, intent(in) :: weight_v11
double precision, intent(in) :: weight_v14
double precision, intent(in) :: weight_v15
double precision, intent(in) :: weight_v16
double precision, intent(in) :: weight_v17
double precision, intent(in) :: weight_v21
double precision, intent(in) :: weight_v22
double precision, intent(in) :: weight_v23
double precision, intent(in) :: weight_v28
double precision, intent(in) :: weight_v29
double precision, intent(in) :: weight_v35


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
v_v7 = y(15)
i_v7 = y(16)
v_v8 = y(17)
i_v8 = y(18)
v_v9 = y(19)
i_v9 = y(20)
v_v10 = y(21)
i_v10 = y(22)
v_v11 = y(23)
i_v11 = y(24)
v_v13 = y(25)
i_v12 = y(26)
v_v14 = y(27)
i_v13 = y(28)
v_v15 = y(29)
i_v14 = y(30)
v_v16 = y(31)
i_v15 = y(32)
v_v17 = y(33)
i_v16 = y(34)
v_v18 = y(35)
i_v17 = y(36)
v_v20 = y(37)
i_v18 = y(38)
v_v21 = y(39)
i_v19 = y(40)
v_v22 = y(41)
i_v20 = y(42)
v_v23 = y(43)
i_v21 = y(44)
v_v24 = y(45)
i_v22 = y(46)
v_v25 = y(47)
i_v23 = y(48)
v_v27 = y(49)
i_v24 = y(50)
v_v28 = y(51)
i_v25 = y(52)
v_v29 = y(53)
i_v26 = y(54)
v_v30 = y(55)
i_v27 = y(56)
v_v31 = y(57)
i_v28 = y(58)
v_v32 = y(59)
i_v29 = y(60)
v_v33 = y(61)
i_v30 = y(62)
v_v34 = y(63)
i_v31 = y(64)
v_v35 = y(65)
i_v32 = y(66)
v_v36 = y(67)
i_v33 = y(68)
v_v37 = y(69)
i_v34 = y(70)
v_v38 = y(71)
i_v35 = y(72)
v_v39 = y(73)
i_v36 = y(74)
v_bI = y(75)
i_v37 = y(76)
v_bIn = v_bI
m_outC = m_max/(exp(r*(V_thr &
     & - 2*v - v_v1 - v_v2 - v_v3 - v_v4 - v_v5)) + 1)
v_bIn_v1 = v_bI
m_outC_v1 = m_max_v1/(exp(r_v1*(V_thr_v1 &
     & - v_v10 - v_v11 - v_v12 - v_v6 - v_v7 - v_v8 - v_v9)) + 1)
v_bIn_v2 = v_bI
m_outC_v2 = m_max_v2/(exp(r_v2*(V_thr_v2 &
     & - v_v13 - v_v14 - v_v15 - v_v16 - v_v17 - v_v18 - v_v19)) + 1)
v_bIn_v3 = v_bI
m_outC_v3 = m_max_v3/(exp(r_v3*(V_thr_v3 &
     & - v_v20 - v_v21 - v_v22 - v_v23 - v_v24 - v_v25 - v_v26)) + 1)
Iext = A*(fsign_1(-onset + t) + 1)&
     & /2 - A*(fsign_1(-dur - onset + t) + 1)/2
m_outC_v4 = m_max_v4/(exp(r_v4*(V_thr_v4 &
     & - v_v27 - v_v28 - v_v29 - v_v30 - v_v31 - v_v32 - v_v33)) + 1)
m_outC_v5 = m_max_v5/(exp(r_v5*(V_thr_v5 &
     & - v_v34 - v_v35 - v_v36 - v_v37 - v_v38 - v_v39)) + 1)
gC = g_input
g_thalC = g_thal_input
bEIC = bEI_input
bEI_thalC = bEI_thal_input
g = gC
bEI = bEIC
m_out = bEI*g*m_outC/connect_reverse_factor
g_v1 = gC
bEI_v1 = bEIC
m_in_v6 = m_out*weight_v6
m_out_v1 = -g_v1*m_outC_v1*(1 - bEI_v1)/connect_reverse_factor_v1
g_v2 = gC
bEI_v2 = bEIC
m_in_v12 = m_out*weight_v12
m_in_v13 = m_out_v1*weight_v13
m_out_v2 = -g_v2*m_outC_v2*(1 - bEI_v2)/connect_reverse_factor_v2
g_v3 = gC
bEI_v3 = bEIC
m_in_v18 = m_out*weight_v18
m_in_v19 = m_out_v1*weight_v19
m_in_v20 = m_out_v2*weight_v20
m_out_v3 = -g_v3*m_outC_v3*(1 - bEI_v3)/connect_reverse_factor_v3
g_thal = g_thalC
bEI_thal = bEI_thalC
m_in_v24 = m_out*weight_v24
m_in_v25 = m_out_v1*weight_v25
m_in_v26 = m_out_v2*weight_v26
m_in_v27 = m_out_v3*weight_v27
m_out_v4 = bEI_thal*g_thal*m_outC_v4/connect_reverse_factor_thal
g_thal_v1 = g_thalC
bEI_thal_v1 = bEI_thalC
m_in_v30 = m_out*weight_v30
m_in_v31 = m_out_v1*weight_v31
m_in_v32 = m_out_v2*weight_v32
m_in_v33 = m_out_v3*weight_v33
m_in_v34 = m_out_v4*weight_v34
m_out_v5 = -g_thal_v1*m_outC_v5*(1 - bEI_thal_v1)&
     & /connect_reverse_factor_thal_v1
m_in = m_out*weight
m_in_v1 = m_out_v1*weight_v1
m_in_v2 = m_out_v2*weight_v2
m_in_v3 = m_out_v3*weight_v3
m_in_v4 = m_out_v4*weight_v4
m_in_v5 = m_out_v5*weight_v5
m_in_v7 = m_out_v1*weight_v7
m_in_v8 = m_out_v2*weight_v8
m_in_v9 = m_out_v3*weight_v9
m_in_v10 = m_out_v4*weight_v10
m_in_v11 = m_out_v5*weight_v11
m_in_v14 = m_out_v2*weight_v14
m_in_v15 = m_out_v3*weight_v15
m_in_v16 = m_out_v4*weight_v16
m_in_v17 = m_out_v5*weight_v17
m_in_v21 = m_out_v3*weight_v21
m_in_v22 = m_out_v4*weight_v22
m_in_v23 = m_out_v5*weight_v23
m_in_v28 = m_out_v4*weight_v28
m_in_v29 = m_out_v5*weight_v29
m_in_v35 = m_out_v5*weight_v35

dy(1) = i
dy(2) = H*m_in/tau - 2*i/tau - v/tau**2
dy(3) = i_v1
dy(4) = H_v1*m_in_v1/tau_v1 - 2*i_v1/tau_v1 - v_v1/tau_v1**2
dy(5) = i_v2
dy(6) = H_v2*m_in_v2/tau_v2 - 2*i_v2/tau_v2 - v_v2/tau_v2**2
dy(7) = i_v3
dy(8) = H_v3*m_in_v3/tau_v3 - 2*i_v3/tau_v3 - v_v3/tau_v3**2
dy(9) = i_v4
dy(10) = H_v4*m_in_v4/tau_v4 - 2*i_v4/tau_v4 - v_v4/tau_v4**2
dy(11) = i_v5
dy(12) = H_v5*m_in_v5/tau_v5 - 2*i_v5/tau_v5 - v_v5/tau_v5**2
dy(13) = i_v6
dy(14) = H_v6*m_in_v6/tau_v6 - 2*i_v6/tau_v6 - v_v6/tau_v6**2
dy(15) = i_v7
dy(16) = H_v7*m_in_v7/tau_v7 - 2*i_v7/tau_v7 - v_v7/tau_v7**2
dy(17) = i_v8
dy(18) = H_v8*m_in_v8/tau_v8 - 2*i_v8/tau_v8 - v_v8/tau_v8**2
dy(19) = i_v9
dy(20) = H_v9*m_in_v9/tau_v9 - 2*i_v9/tau_v9 - v_v9/tau_v9**2
dy(21) = i_v10
dy(22) = H_v10*m_in_v10/tau_v10 - 2*i_v10/tau_v10 - v_v10/tau_v10**2
dy(23) = i_v11
dy(24) = H_v11*m_in_v11/tau_v11 - 2*i_v11/tau_v11 - v_v11/tau_v11**2
dy(25) = i_v12
dy(26) = H_v12*m_in_v12/tau_v12 - 2*i_v12/tau_v12 - v_v13/tau_v12**2
dy(27) = i_v13
dy(28) = H_v13*m_in_v13/tau_v13 - 2*i_v13/tau_v13 - v_v14/tau_v13**2
dy(29) = i_v14
dy(30) = H_v14*m_in_v14/tau_v14 - 2*i_v14/tau_v14 - v_v15/tau_v14**2
dy(31) = i_v15
dy(32) = H_v15*m_in_v15/tau_v15 - 2*i_v15/tau_v15 - v_v16/tau_v15**2
dy(33) = i_v16
dy(34) = H_v16*m_in_v16/tau_v16 - 2*i_v16/tau_v16 - v_v17/tau_v16**2
dy(35) = i_v17
dy(36) = H_v17*m_in_v17/tau_v17 - 2*i_v17/tau_v17 - v_v18/tau_v17**2
dy(37) = i_v18
dy(38) = H_v18*m_in_v18/tau_v18 - 2*i_v18/tau_v18 - v_v20/tau_v18**2
dy(39) = i_v19
dy(40) = H_v19*m_in_v19/tau_v19 - 2*i_v19/tau_v19 - v_v21/tau_v19**2
dy(41) = i_v20
dy(42) = H_v20*m_in_v20/tau_v20 - 2*i_v20/tau_v20 - v_v22/tau_v20**2
dy(43) = i_v21
dy(44) = H_v21*m_in_v21/tau_v21 - 2*i_v21/tau_v21 - v_v23/tau_v21**2
dy(45) = i_v22
dy(46) = H_v22*m_in_v22/tau_v22 - 2*i_v22/tau_v22 - v_v24/tau_v22**2
dy(47) = i_v23
dy(48) = H_v23*m_in_v23/tau_v23 - 2*i_v23/tau_v23 - v_v25/tau_v23**2
dy(49) = i_v24
dy(50) = H_v24*m_in_v24/tau_v24 - 2*i_v24/tau_v24 - v_v27/tau_v24**2
dy(51) = i_v25
dy(52) = H_v25*m_in_v25/tau_v25 - 2*i_v25/tau_v25 - v_v28/tau_v25**2
dy(53) = i_v26
dy(54) = H_v26*m_in_v26/tau_v26 - 2*i_v26/tau_v26 - v_v29/tau_v26**2
dy(55) = i_v27
dy(56) = H_v27*m_in_v27/tau_v27 - 2*i_v27/tau_v27 - v_v30/tau_v27**2
dy(57) = i_v28
dy(58) = H_v28*m_in_v28/tau_v28 - 2*i_v28/tau_v28 - v_v31/tau_v28**2
dy(59) = i_v29
dy(60) = H_v29*m_in_v29/tau_v29 - 2*i_v29/tau_v29 - v_v32/tau_v29**2
dy(61) = i_v30
dy(62) = H_v30*Iext/tau_v30 - 2*i_v30/tau_v30 - v_v33/tau_v30**2
dy(63) = i_v31
dy(64) = H_v31*m_in_v30/tau_v31 - 2*i_v31/tau_v31 - v_v34/tau_v31**2
dy(65) = i_v32
dy(66) = H_v32*m_in_v31/tau_v32 - 2*i_v32/tau_v32 - v_v35/tau_v32**2
dy(67) = i_v33
dy(68) = H_v33*m_in_v32/tau_v33 - 2*i_v33/tau_v33 - v_v36/tau_v33**2
dy(69) = i_v34
dy(70) = H_v34*m_in_v33/tau_v34 - 2*i_v34/tau_v34 - v_v37/tau_v34**2
dy(71) = i_v35
dy(72) = H_v35*m_in_v34/tau_v35 - 2*i_v35/tau_v35 - v_v38/tau_v35**2
dy(73) = i_v36
dy(74) = H_v36*m_in_v35/tau_v36 - 2*i_v36/tau_v36 - v_v39/tau_v36**2
dy(75) = i_v37
dy(76) = H_v37*bI/tau_v37 - 2*i_v37/tau_v37 - v_bI/tau_v37**2

end subroutine



function fsign_1(x)

implicit none

double precision :: fsign_1
double precision :: x
double precision :: a


a = 1.0
fsign_1 = sign(a,x)

end function fsign_1
    

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
     & args(35), args(36), args(37), args(38), args(39), args(40), &
     & args(41), args(42), args(43), args(44), args(45), args(46), &
     & args(47), args(48), args(49), args(50), args(51), args(52), &
     & args(53), args(54), args(55), args(56), args(57), args(58), &
     & args(59), args(60), args(61), args(62), args(63), args(64), &
     & args(65), args(66), args(67), args(68), args(69), args(70), &
     & args(71), args(72), args(73), args(74), args(75), args(76), &
     & args(77), args(78), args(79), args(80), args(81), args(82), &
     & args(83), args(84), args(85), args(86), args(87), args(88), &
     & args(89), args(90), args(91), args(92), args(93), args(94), &
     & args(95), args(96), args(97), args(98), args(99), args(100), &
     & args(101), args(102), args(103), args(104), args(105), args(106)&
     & , args(107), args(108), args(109), args(110), args(111), &
     & args(112), args(113), args(114), args(115), args(116), args(117)&
     & , args(118), args(119), args(120), args(121), args(122), &
     & args(123), args(124), args(125), args(126), args(127), args(128)&
     & , args(129), args(130), args(131), args(132), args(133), &
     & args(134), args(135), args(136), args(137), args(138), args(139)&
     & , args(140), args(141), args(142), args(143), args(144), &
     & args(145), args(146), args(147), args(148), args(149), args(150)&
     & , args(151), args(152))

end subroutine func


subroutine stpnt(ndim, y, args, t)

implicit None
integer, intent(in) :: ndim
double precision, intent(inout) :: y(ndim), args(*)
double precision, intent(in) :: t

args(1) = 0.006  ! tau
args(2) = 1.0  ! H
args(3) = 0.003  ! tau_v1
args(4) = 1.0  ! H_v1
args(5) = 0.02  ! tau_v2
args(6) = 1.0  ! H_v2
args(7) = 0.015  ! tau_v3
args(8) = 1.0  ! H_v3
args(9) = 0.003  ! tau_v4
args(15) = 1.0  ! H_v4
args(16) = 0.003  ! tau_v5
args(17) = 1.0  ! H_v5
args(18) = 0.006  ! tau_v6
args(19) = 1.0  ! H_v6
args(20) = 0.003  ! tau_v7
args(21) = 1.0  ! H_v7
args(22) = 0.02  ! tau_v8
args(23) = 1.0  ! H_v8
args(24) = 0.015  ! tau_v9
args(25) = 1.0  ! H_v9
args(26) = 0.003  ! tau_v10
args(27) = 1.0  ! H_v10
args(28) = 0.003  ! tau_v11
args(29) = 1.0  ! H_v11
args(30) = 0.006  ! tau_v12
args(31) = 1.0  ! H_v12
args(32) = 0.003  ! tau_v13
args(33) = 1.0  ! H_v13
args(34) = 0.02  ! tau_v14
args(35) = 1.0  ! H_v14
args(36) = 0.015  ! tau_v15
args(37) = 1.0  ! H_v15
args(38) = 0.003  ! tau_v16
args(39) = 1.0  ! H_v16
args(40) = 0.003  ! tau_v17
args(41) = 1.0  ! H_v17
args(42) = 0.006  ! tau_v18
args(43) = 1.0  ! H_v18
args(44) = 0.003  ! tau_v19
args(45) = 1.0  ! H_v19
args(46) = 0.02  ! tau_v20
args(47) = 1.0  ! H_v20
args(48) = 0.015  ! tau_v21
args(49) = 1.0  ! H_v21
args(50) = 0.003  ! tau_v22
args(51) = 1.0  ! H_v22
args(52) = 0.003  ! tau_v23
args(53) = 1.0  ! H_v23
args(54) = 0.006  ! tau_v24
args(55) = 1.0  ! H_v24
args(56) = 0.003  ! tau_v25
args(57) = 1.0  ! H_v25
args(58) = 0.02  ! tau_v26
args(59) = 1.0  ! H_v26
args(60) = 0.015  ! tau_v27
args(61) = 1.0  ! H_v27
args(62) = 0.003  ! tau_v28
args(63) = 1.0  ! H_v28
args(64) = 0.003  ! tau_v29
args(65) = 1.0  ! H_v29
args(66) = 0.003  ! tau_v30
args(67) = 1.0  ! H_v30
args(68) = 0.006  ! tau_v31
args(69) = 1.0  ! H_v31
args(70) = 0.003  ! tau_v32
args(71) = 1.0  ! H_v32
args(72) = 0.02  ! tau_v33
args(73) = 1.0  ! H_v33
args(74) = 0.015  ! tau_v34
args(75) = 1.0  ! H_v34
args(76) = 0.003  ! tau_v35
args(77) = 1.0  ! H_v35
args(78) = 0.003  ! tau_v36
args(79) = 1.0  ! H_v36
args(80) = 0.003  ! tau_v37
args(81) = 200.0  ! bI
args(82) = 1.0  ! H_v37
args(83) = 32.10540543  ! V_thr
args(84) = 0.12782346  ! r
args(85) = 31.39696397  ! m_max
args(86) = 40.03107351  ! V_thr_v1
args(87) = 0.0  ! v_v12
args(88) = 0.14218422  ! r_v1
args(89) = 166.82960408  ! m_max_v1
args(90) = 42.01276379  ! V_thr_v2
args(91) = 0.0  ! v_v19
args(92) = 0.07937015  ! r_v2
args(93) = 56.95305832  ! m_max_v2
args(94) = 37.86409387  ! V_thr_v3
args(95) = 0.0  ! v_v26
args(96) = 0.0704119  ! r_v3
args(97) = 38.52689646  ! m_max_v3
args(98) = 1.0  ! onset
args(99) = 1.0  ! dur
args(100) = 50.0  ! A
args(101) = 40.0  ! V_thr_v4
args(102) = 0.1  ! r_v4
args(103) = 30.0  ! m_max_v4
args(104) = 40.0  ! V_thr_v5
args(105) = 0.1  ! r_v5
args(106) = 30.0  ! m_max_v5
args(107) = 100.0  ! g_input
args(108) = 200.0  ! g_thal_input
args(109) = 0.5  ! bEI_input
args(110) = 0.5  ! bEI_thal_input
args(111) = [6448]  ! connect_reverse_factor
args(112) = 126.85331071899701  ! weight_v6
args(113) = [6448]  ! connect_reverse_factor_v1
args(114) = 34.73829288056395  ! weight_v12
args(115) = 8.711911694438868  ! weight_v13
args(116) = [6448]  ! connect_reverse_factor_v2
args(117) = 20.2590912478185  ! weight_v18
args(118) = 5.156935882352942  ! weight_v19
args(119) = 9.086982935153586  ! weight_v20
args(120) = [6448]  ! connect_reverse_factor_v3
args(121) = 0.0  ! weight_v24
args(122) = 0.0  ! weight_v25
args(123) = 0.0  ! weight_v26
args(124) = 0.0  ! weight_v27
args(125) = [6448]  ! connect_reverse_factor_thal
args(126) = 0.0  ! weight_v30
args(127) = 0.0  ! weight_v31
args(128) = 0.0  ! weight_v32
args(129) = 0.0  ! weight_v33
args(130) = 0.0  ! weight_v34
args(131) = [6448]  ! connect_reverse_factor_thal_v1
args(132) = 53.61853519164954  ! weight
args(133) = 9.920668080279233  ! weight_v1
args(134) = 5.541078936916653  ! weight_v2
args(135) = 3.958166666666667  ! weight_v3
args(136) = 571.71583  ! weight_v4
args(137) = 0.0  ! weight_v5
args(138) = 7.405289065743945  ! weight_v7
args(139) = 5.908537795623369  ! weight_v8
args(140) = 2.9094999999999995  ! weight_v9
args(141) = 35.69699  ! weight_v10
args(142) = 0.0  ! weight_v11
args(143) = 1.1365268618155135  ! weight_v14
args(144) = 7.097935153583618  ! weight_v15
args(145) = 2.3520000000000003  ! weight_v16
args(146) = 0.0  ! weight_v17
args(147) = 5.355  ! weight_v21
args(148) = 0.0  ! weight_v22
args(149) = 0.0  ! weight_v23
args(150) = 0.0  ! weight_v28
args(151) = 0.0  ! weight_v29
args(152) = 0.0  ! weight_v35
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
y(15) = 0.0  ! v_v7
y(16) = 0.0  ! i_v7
y(17) = 0.0  ! v_v8
y(18) = 0.0  ! i_v8
y(19) = 0.0  ! v_v9
y(20) = 0.0  ! i_v9
y(21) = 0.0  ! v_v10
y(22) = 0.0  ! i_v10
y(23) = 0.0  ! v_v11
y(24) = 0.0  ! i_v11
y(25) = 0.0  ! v_v13
y(26) = 0.0  ! i_v12
y(27) = 0.0  ! v_v14
y(28) = 0.0  ! i_v13
y(29) = 0.0  ! v_v15
y(30) = 0.0  ! i_v14
y(31) = 0.0  ! v_v16
y(32) = 0.0  ! i_v15
y(33) = 0.0  ! v_v17
y(34) = 0.0  ! i_v16
y(35) = 0.0  ! v_v18
y(36) = 0.0  ! i_v17
y(37) = 0.0  ! v_v20
y(38) = 0.0  ! i_v18
y(39) = 0.0  ! v_v21
y(40) = 0.0  ! i_v19
y(41) = 0.0  ! v_v22
y(42) = 0.0  ! i_v20
y(43) = 0.0  ! v_v23
y(44) = 0.0  ! i_v21
y(45) = 0.0  ! v_v24
y(46) = 0.0  ! i_v22
y(47) = 0.0  ! v_v25
y(48) = 0.0  ! i_v23
y(49) = 0.0  ! v_v27
y(50) = 0.0  ! i_v24
y(51) = 0.0  ! v_v28
y(52) = 0.0  ! i_v25
y(53) = 0.0  ! v_v29
y(54) = 0.0  ! i_v26
y(55) = 0.0  ! v_v30
y(56) = 0.0  ! i_v27
y(57) = 0.0  ! v_v31
y(58) = 0.0  ! i_v28
y(59) = 0.0  ! v_v32
y(60) = 0.0  ! i_v29
y(61) = 0.0  ! v_v33
y(62) = 0.0  ! i_v30
y(63) = 0.0  ! v_v34
y(64) = 0.0  ! i_v31
y(65) = 0.0  ! v_v35
y(66) = 0.0  ! i_v32
y(67) = 0.0  ! v_v36
y(68) = 0.0  ! i_v33
y(69) = 0.0  ! v_v37
y(70) = 0.0  ! i_v34
y(71) = 0.0  ! v_v38
y(72) = 0.0  ! i_v35
y(73) = 0.0  ! v_v39
y(74) = 0.0  ! i_v36
y(75) = 0.0  ! v_bI
y(76) = 0.0  ! i_v37

end subroutine stpnt



subroutine bcnd
end subroutine bcnd


subroutine icnd
end subroutine icnd


subroutine fopt
end subroutine fopt


subroutine pvls
end subroutine pvls
