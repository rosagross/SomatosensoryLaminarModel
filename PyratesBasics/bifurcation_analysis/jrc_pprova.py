from numpy import pi, sqrt
from numpy import exp


def jrc_rhs(t,y,dy,u,tau,H,tau_v1,H_v1,tau_v2,H_v2,tau_v3,H_v3,V_thr,r,m_max,weight_v2,V_thr_v1,r_v1,m_max_v1,weight_v3,V_thr_v2,r_v2,m_max_v2,weight,weight_v1):


	v = y[0]
	i = y[1]
	v_v1 = y[2]
	i_v1 = y[3]
	v_v2 = y[4]
	i_v2 = y[5]
	v_v3 = y[6]
	i_v3 = y[7]
	m_out_v2 = 2.0*m_max/(exp(r*(V_thr - v - v_v1)) + 1)
	m_in_v2 = m_out_v2*weight_v2
	m_out = 2.0*m_max_v1/(exp(r_v1*(V_thr_v1 - v_v2)) + 1)
	m_in_v3 = m_out_v2*weight_v3
	m_out_v1 = 2.0*m_max_v2/(exp(r_v2*(V_thr_v2 - v_v3)) + 1)
	m_in = m_out*weight
	m_in_v1 = m_out_v1*weight_v1
	
	dy[0] = i
	dy[1] = H*(m_in + u)/tau - 2*i/tau - v/tau**2
	dy[2] = i_v1
	dy[3] = H_v1*m_in_v1/tau_v1 - 2*i_v1/tau_v1 - v_v1/tau_v1**2
	dy[4] = i_v2
	dy[5] = H_v2*m_in_v2/tau_v2 - 2*i_v2/tau_v2 - v_v2/tau_v2**2
	dy[6] = i_v3
	dy[7] = H_v3*m_in_v3/tau_v3 - 2*i_v3/tau_v3 - v_v3/tau_v3**2

	return dy