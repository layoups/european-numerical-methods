function y = lu_find_y( y, a )
asq    = a*a; 
y(2)   = 1+2*a; 
Nminus = 1; 
Nplus  = length(y); 
for n=3:(Nplus-1)
  y(n) = 1+2*a - asq/y(n-1); 
  if( y(n)==0 ) error( 'failure / singular ...' ); end
end


