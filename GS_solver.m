function u = GS_solver( u, b, a, n_eps )

error = 1e6; 
nelts = length(u); 

while( error>n_eps )
  error = 0; 
 
  for n=2:(nelts-1)
    y     = ( b(n) + a * ( u(n-1) + u(n+1) ) )/(1+2*a); 
    error = error + (u(n)-y)^2;
    u(n)  = y;                                
  end
    

end

