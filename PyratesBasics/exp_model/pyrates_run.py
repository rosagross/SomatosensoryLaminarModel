from numpy import pi, sqrt
from numpy import dot
from numpy import exp
from numpy import sign


def vector_field(t,y,dy,tau,H,tau_v1,H_v1,tau_v2,H_v2,tau_v3,H_v3,tau_v4,H_v4,tau_v5,H_v5,tau_v6,H_v6,tau_v7,H_v7,tau_v8,H_v8,tau_v9,H_v9,tau_v10,H_v10,tau_v11,H_v11,tau_v12,H_v12,tau_v13,H_v13,tau_v14,H_v14,tau_v15,H_v15,tau_v16,H_v16,tau_v17,H_v17,tau_v18,H_v18,tau_v19,bI,H_v19,weight,V_thr,r,m_max,onset,dur,A,source_idx_v4,weight_v7,source_idx_v5,weight_v8,source_idx_v6,weight_v9,source_idx_v7,weight_v10,V_thr_v1,r_v1,m_max_v1,source_idx_v8,weight_v13,source_idx_v9,weight_v14,source_idx_v10,weight_v15,source_idx_v11,weight_v16,weight_v17,V_thr_v2,r_v2,m_max_v2,source_idx,weight_v1,source_idx_v1,weight_v2,source_idx_v2,weight_v3,source_idx_v3,weight_v4,weight_v5,weight_v6,weight_v11,weight_v12,weight_v18):


	v = y[0:4]
	i = y[4:8]
	v_v1 = y[8:12]
	i_v1 = y[12:16]
	v_v2 = y[16:20]
	i_v2 = y[20:24]
	v_v3 = y[24:28]
	i_v3 = y[28:32]
	v_v4 = y[32:36]
	i_v4 = y[36:40]
	v_v5 = y[40:44]
	i_v5 = y[44:48]
	v_v6 = y[48]
	i_v6 = y[49]
	v_v7 = y[50]
	i_v7 = y[51]
	v_v8 = y[52]
	i_v8 = y[53]
	v_v9 = y[54]
	i_v9 = y[55]
	v_v10 = y[56]
	i_v10 = y[57]
	v_v11 = y[58]
	i_v11 = y[59]
	v_v12 = y[60]
	i_v12 = y[61]
	v_v13 = y[62]
	i_v13 = y[63]
	v_v14 = y[64]
	i_v14 = y[65]
	v_v15 = y[66]
	i_v15 = y[67]
	v_v16 = y[68]
	i_v16 = y[69]
	v_v17 = y[70]
	i_v17 = y[71]
	v_v18 = y[72]
	i_v18 = y[73]
	v_bI = y[74]
	i_v19 = y[75]
	v_bIn = dot(weight, v_bI)
	m_out = m_max/(exp(r*(V_thr - 2*v - v_v1 - v_v2 - v_v3 - v_v4 - v_v5)) + 1)
	Iext = A*(sign(-onset + t) + 1)/2 - A*(sign(-dur - onset + t) + 1)/2
	m_in_v6 = weight_v7*m_out[source_idx_v4]
	m_in_v7 = weight_v8*m_out[source_idx_v5]
	m_in_v8 = weight_v9*m_out[source_idx_v6]
	m_in_v9 = weight_v10*m_out[source_idx_v7]
	m_out_v1 = m_max_v1/(exp(r_v1*(V_thr_v1 - v_v10 - v_v11 - v_v12 - v_v6 - v_v7 - v_v8 - v_v9)) + 1)
	m_in_v12 = weight_v13*m_out[source_idx_v8]
	m_in_v13 = weight_v14*m_out[source_idx_v9]
	m_in_v14 = weight_v15*m_out[source_idx_v10]
	m_in_v15 = weight_v16*m_out[source_idx_v11]
	m_in_v16 = m_out_v1*weight_v17
	m_out_v2 = m_max_v2/(exp(r_v2*(V_thr_v2 - v_v13 - v_v14 - v_v15 - v_v16 - v_v17 - v_v18)) + 1)
	m_in = dot(weight_v1, m_out[source_idx])
	m_in_v1 = dot(weight_v2, m_out[source_idx_v1])
	m_in_v2 = dot(weight_v3, m_out[source_idx_v2])
	m_in_v3 = dot(weight_v4, m_out[source_idx_v3])
	m_in_v4 = dot(weight_v5, m_out_v1)
	m_in_v5 = dot(weight_v6, m_out_v2)
	m_in_v10 = m_out_v1*weight_v11
	m_in_v11 = m_out_v2*weight_v12
	m_in_v17 = m_out_v2*weight_v18
	
	dy[0:4] = i
	dy[4:8] = H*m_in/tau - 2*i/tau - v/tau**2
	dy[8:12] = i_v1
	dy[12:16] = H_v1*m_in_v1/tau_v1 - 2*i_v1/tau_v1 - v_v1/tau_v1**2
	dy[16:20] = i_v2
	dy[20:24] = H_v2*m_in_v2/tau_v2 - 2*i_v2/tau_v2 - v_v2/tau_v2**2
	dy[24:28] = i_v3
	dy[28:32] = H_v3*m_in_v3/tau_v3 - 2*i_v3/tau_v3 - v_v3/tau_v3**2
	dy[32:36] = i_v4
	dy[36:40] = H_v4*m_in_v4/tau_v4 - 2*i_v4/tau_v4 - v_v4/tau_v4**2
	dy[40:44] = i_v5
	dy[44:48] = H_v5*m_in_v5/tau_v5 - 2*i_v5/tau_v5 - v_v5/tau_v5**2
	dy[48] = i_v6
	dy[49] = H_v6*m_in_v6/tau_v6 - 2*i_v6/tau_v6 - v_v6/tau_v6**2
	dy[50] = i_v7
	dy[51] = H_v7*m_in_v7/tau_v7 - 2*i_v7/tau_v7 - v_v7/tau_v7**2
	dy[52] = i_v8
	dy[53] = H_v8*m_in_v8/tau_v8 - 2*i_v8/tau_v8 - v_v8/tau_v8**2
	dy[54] = i_v9
	dy[55] = H_v9*m_in_v9/tau_v9 - 2*i_v9/tau_v9 - v_v9/tau_v9**2
	dy[56] = i_v10
	dy[57] = H_v10*m_in_v10/tau_v10 - 2*i_v10/tau_v10 - v_v10/tau_v10**2
	dy[58] = i_v11
	dy[59] = H_v11*m_in_v11/tau_v11 - 2*i_v11/tau_v11 - v_v11/tau_v11**2
	dy[60] = i_v12
	dy[61] = H_v12*Iext/tau_v12 - 2*i_v12/tau_v12 - v_v12/tau_v12**2
	dy[62] = i_v13
	dy[63] = H_v13*m_in_v12/tau_v13 - 2*i_v13/tau_v13 - v_v13/tau_v13**2
	dy[64] = i_v14
	dy[65] = H_v14*m_in_v13/tau_v14 - 2*i_v14/tau_v14 - v_v14/tau_v14**2
	dy[66] = i_v15
	dy[67] = H_v15*m_in_v14/tau_v15 - 2*i_v15/tau_v15 - v_v15/tau_v15**2
	dy[68] = i_v16
	dy[69] = H_v16*m_in_v15/tau_v16 - 2*i_v16/tau_v16 - v_v16/tau_v16**2
	dy[70] = i_v17
	dy[71] = H_v17*m_in_v16/tau_v17 - 2*i_v17/tau_v17 - v_v17/tau_v17**2
	dy[72] = i_v18
	dy[73] = H_v18*m_in_v17/tau_v18 - 2*i_v18/tau_v18 - v_v18/tau_v18**2
	dy[74] = i_v19
	dy[75] = H_v19*bI/tau_v19 - 2*i_v19/tau_v19 - v_bI/tau_v19**2

	return dy