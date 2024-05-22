from numpy import pi, sqrt
from numpy import dot
from numpy import exp


def vector_field(t,y,dy,tau,H,weight,V_thr,r,m_max):


	v = y[0:11]
	i = y[11:22]
	v_in = dot(weight, v)
	m = 2.0*m_max/(exp(r*(V_thr - 1000*v_in)) + 1)
	
	dy[0:11] = i
	dy[11:22] = H*m/tau - 2*i/tau - v/tau**2

	return dy