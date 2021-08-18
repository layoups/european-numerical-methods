close all; drawnow; 
clear all; 
clear functions; 

% parameters
r     = 0.1;
sigma = 0.45; 
E     = 10.0; 
T     = 0.33;                
k     = r/(0.5*sigma^2); 

% x grid 
SRight = 16.0;                % +\infinity 
SLeft  = 1e-9;                % -\infinity 
xLeft  = log(SLeft/E); 
xRight = log(SRight/E); 
Nx = 4000;
 
for p=1:4
    dx = p*(xRight-xLeft)/Nx; 
    % alpha, delta t
    a  = 0.5; 
    dt = a*dx^2;

    % tau range
    tau_Max = (0.5*sigma^2)*T; 
    M       = ceil(tau_Max/dt);

    fprintf('running dx = %10.6f\n',dx); 
    [u,xgrid] = crank_fd_GS(@tran_payoff_call, @u_m_inf_call, @u_p_inf_call, r, sigma, xLeft, xRight, Nx, tau_Max, M );
    [u_1,xgrid_1] = crank_fd_LU(@tran_payoff_call, @u_m_inf_call, @u_p_inf_call, r, sigma, xLeft, xRight, Nx, tau_Max, M );


    % return to financial variables
    S   = E*exp( xgrid ); 
    t   = 0;                       
    tau = 0.5*(sigma^2)*(T-t);     

    Spow = (S.^(0.5*(1-k))); 
    Smat = repmat( Spow(:).', [M+1, 1] ); 
    V    = (E^(0.5*(1+k))) * Smat * exp( -(1/4)*((k+1)^2)*tau ).*u; 

    % explicit B-S solution
    [C, P] = blsprice(S,E,r,T,sigma);
    figure; 
    ns=plot( S, V(end,:), 'x' ); 
    grid on; 
    hold on; 
    xlabel( 'S' ); 
    as=plot( S, C, '-og' ); ylabel('C'); 
    title( ['The European Call Example.  dx=',num2str(dx)] );   
    legend( [ ns,as ], {'numerical solution', 'analytic solution'}, 'location', 'north' ); 

    S_1 = E*exp( xgrid_1 ); 
    Spow_1 = (S_1.^(0.5*(1-k))); 
    Smat_1 = repmat( Spow_1(:).', [M+1, 1] ); 
    V_1    = (E^(0.5*(1+k))) * Smat_1 * exp( -(1/4)*((k+1)^2)*tau ).*u_1; 

    % explicit B-S solution
    [C, P] = blsprice(S_1,E,r,T,sigma);
    figure; 
    ns=plot( S_1, V_1(end,:), 'x' ); 
    grid on; 
    hold on; 
    xlabel( 'S_1' ); 
    as=plot( S_1, C, '-og' ); ylabel('C'); 
    title( ['The European Call Example.  dx=',num2str(dx)] );   
    legend( [ ns,as ], {'numerical solution', 'analytic solution'}, 'location', 'north' ); 


    % for the sake of comparison with values in book
    Sgrid = [ 0.00, 2, 4, 6, 7, 8, 9, 10:16 ]; 
    Vgrid = interp1( S, V(end,:), Sgrid, 'linear', 'extrap' )
    Vgrid_1 = interp1( S_1, V_1(end,:), Sgrid, 'linear', 'extrap' )
end
