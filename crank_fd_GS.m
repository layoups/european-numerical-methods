function [u,xgrid] = crank_fd_GS( payoff_fn, u_m_inf, u_p_inf, r, sigma, xLeft, xRight, Nx, tau_Max, M )

xgrid = linspace( xLeft, xRight, Nx+1 ); 
dx    = (xRight-xLeft)/Nx; 
dt    = tau_Max/M; 

a        = dt/(dx*dx); 
a2       = 0.5*a; 
n_eps    = 1.e-8; 
k        = r/(0.5*sigma^2); 

tau  = 0.0; 
oldu = payoff_fn( xgrid, tau, k ); oldu = oldu(:).'; 

u      = zeros(M+1,Nx+1);   
u(1,:) = oldu; 
newu   = zeros(1,Nx+1); 
b      = zeros(1,Nx+1); 
for m=1:M
  tau = m*dt; 

  b( 2:(end-1) ) = (1-a)*oldu( 2:(end-1) ) + a2*( oldu( 3:end ) + oldu( 1:(end-2) ) );
  
  newu(1)   = u_m_inf( xgrid(1),   tau, k ); 
  newu(end) = u_p_inf( xgrid(end), tau, k ); 

  newu = GS_solver( newu, b, a2, n_eps );
  oldu = newu;   

  u(m+1,:) = newu; 
end


