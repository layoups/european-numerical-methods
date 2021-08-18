function [u,xgrid] = crank_fd_LU( payoff_fn, u_m_inf, u_p_inf, r, sigma, xLeft, xRight, Nx, tau_Max, M )

xgrid = linspace( xLeft, xRight, Nx+1 ); 
dx    = (xRight-xLeft)/Nx; 
dt    = tau_Max/M; 

a     = dt/(dx*dx); 
a2    = 0.5*a; 
k     = r/(0.5*sigma^2); 

tau  = 0.0; 
oldu = payoff_fn( xgrid, tau, k ); 
oldu = oldu(:).'; 

y = lu_find_y( zeros(1,Nx+1), a2 ); 

u      = zeros(M+1,Nx+1);   
u(1,:) = oldu; 
newu   = zeros(1,Nx+1); 
b      = zeros(1,Nx+1); 
for m=1:M
  tau = m*dt; 

  b( 2:(end-1) ) = (1-a)*oldu( 2:(end-1) ) + a2*( oldu( 3:end ) + oldu( 1:(end-2) ) );
  
  oldu(1)   = u_m_inf( xgrid(1),   tau, k ); 
  oldu(end) = u_p_inf( xgrid(end), tau, k ); 

  b(2)     = b(2)     + a2*oldu(1); 
  b(end-1) = b(end-1) + a2*oldu(end); 
  
  [newu] = lu_solver( oldu, b, y, a2 );
  
  oldu = newu;   

  u(m+1,:) = newu; 
end


