from numpy import pi, sqrt
from numpy import exp


def vector_field(t,y,dy,bI,tau,H,bI_v1,tau_v1,H_v1,bI_v2,tau_v2,H_v2,bI_v3,tau_v3,H_v3,bI_v4,tau_v4,H_v4,bI_v5,tau_v5,H_v5,bI_v6,tau_v6,H_v6,bI_v7,tau_v7,H_v7,bI_v8,tau_v8,H_v8,bI_v9,tau_v9,H_v9,bI_v10,tau_v10,H_v10,bI_v11,tau_v11,H_v11,bI_v12,tau_v12,H_v12,V_thr,r,m_max,source_idx,weight,source_idx_v1,weight_v1,source_idx_v2,weight_v2,source_idx_v3,weight_v3,source_idx_v4,weight_v4,source_idx_v5,weight_v5,source_idx_v6,weight_v6,source_idx_v7,weight_v7,source_idx_v8,weight_v8,source_idx_v9,weight_v9,source_idx_v10,weight_v10,source_idx_v11,weight_v11,source_idx_v12,weight_v12):


	v = y[0:13]
	i = y[13:26]
	v_v1 = y[26:39]
	i_v1 = y[39:52]
	v_v2 = y[52:65]
	i_v2 = y[65:78]
	v_v3 = y[78:91]
	i_v3 = y[91:104]
	v_v4 = y[104:117]
	i_v4 = y[117:130]
	v_v5 = y[130:143]
	i_v5 = y[143:156]
	v_v6 = y[156:169]
	i_v6 = y[169:182]
	v_v7 = y[182:195]
	i_v7 = y[195:208]
	v_v8 = y[208:221]
	i_v8 = y[221:234]
	v_v9 = y[234:247]
	i_v9 = y[247:260]
	v_v10 = y[260:273]
	i_v10 = y[273:286]
	v_v11 = y[286:299]
	i_v11 = y[299:312]
	v_v12 = y[312:325]
	i_v12 = y[325:338]
	m_out = 2.0*m_max/(exp(r*(V_thr - v - v_v1 - v_v10 - v_v11 - v_v12 - v_v2 - v_v3 - v_v4 - v_v5 - v_v6 - v_v7 - v_v8 - v_v9)) + 1)
	m_in = weight*m_out[source_idx]
	m_in_v1 = weight_v1*m_out[source_idx_v1]
	m_in_v2 = weight_v2*m_out[source_idx_v2]
	m_in_v3 = weight_v3*m_out[source_idx_v3]
	m_in_v4 = weight_v4*m_out[source_idx_v4]
	m_in_v5 = weight_v5*m_out[source_idx_v5]
	m_in_v6 = weight_v6*m_out[source_idx_v6]
	m_in_v7 = weight_v7*m_out[source_idx_v7]
	m_in_v8 = weight_v8*m_out[source_idx_v8]
	m_in_v9 = weight_v9*m_out[source_idx_v9]
	m_in_v10 = weight_v10*m_out[source_idx_v10]
	m_in_v11 = weight_v11*m_out[source_idx_v11]
	m_in_v12 = weight_v12*m_out[source_idx_v12]
	
	dy[0:13] = i
	dy[13:26] = H*(bI + m_in)/tau - 2*i/tau - v/tau**2
	dy[26:39] = i_v1
	dy[39:52] = H_v1*(bI_v1 + m_in_v1)/tau_v1 - 2*i_v1/tau_v1 - v_v1/tau_v1**2
	dy[52:65] = i_v2
	dy[65:78] = H_v2*(bI_v2 + m_in_v2)/tau_v2 - 2*i_v2/tau_v2 - v_v2/tau_v2**2
	dy[78:91] = i_v3
	dy[91:104] = H_v3*(bI_v3 + m_in_v3)/tau_v3 - 2*i_v3/tau_v3 - v_v3/tau_v3**2
	dy[104:117] = i_v4
	dy[117:130] = H_v4*(bI_v4 + m_in_v4)/tau_v4 - 2*i_v4/tau_v4 - v_v4/tau_v4**2
	dy[130:143] = i_v5
	dy[143:156] = H_v5*(bI_v5 + m_in_v5)/tau_v5 - 2*i_v5/tau_v5 - v_v5/tau_v5**2
	dy[156:169] = i_v6
	dy[169:182] = H_v6*(bI_v6 + m_in_v6)/tau_v6 - 2*i_v6/tau_v6 - v_v6/tau_v6**2
	dy[182:195] = i_v7
	dy[195:208] = H_v7*(bI_v7 + m_in_v7)/tau_v7 - 2*i_v7/tau_v7 - v_v7/tau_v7**2
	dy[208:221] = i_v8
	dy[221:234] = H_v8*(bI_v8 + m_in_v8)/tau_v8 - 2*i_v8/tau_v8 - v_v8/tau_v8**2
	dy[234:247] = i_v9
	dy[247:260] = H_v9*(bI_v9 + m_in_v9)/tau_v9 - 2*i_v9/tau_v9 - v_v9/tau_v9**2
	dy[260:273] = i_v10
	dy[273:286] = H_v10*(bI_v10 + m_in_v10)/tau_v10 - 2*i_v10/tau_v10 - v_v10/tau_v10**2
	dy[286:299] = i_v11
	dy[299:312] = H_v11*(bI_v11 + m_in_v11)/tau_v11 - 2*i_v11/tau_v11 - v_v11/tau_v11**2
	dy[312:325] = i_v12
	dy[325:338] = H_v12*(bI_v12 + m_in_v12)/tau_v12 - 2*i_v12/tau_v12 - v_v12/tau_v12**2

	return dy