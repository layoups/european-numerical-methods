function u_inf = u_p_inf_call(x,tau,k)
u_inf = exp( 0.25*((k+1)^2)*tau ) * ( exp( 0.5*(k+1)*x ) - exp( 0.5*(k-1)*x )*exp( -k*tau ) ); 











