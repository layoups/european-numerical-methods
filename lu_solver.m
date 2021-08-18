function u = lu_solver( u, b, y, a )
q = zeros(size(u)); 

Nminus = 1; 
Nplus  = length(y); 

q(2) = b(2);
for n=3:(Nplus-1) 
  q(n) = b(n) + a*q(n-1)/y(n-1);
end

u(Nplus-1) = q(Nplus-1)/y(Nplus-1); 

for n=(Nplus-2):-1:(Nminus+1) 
  u(n) = (q(n)+a*u(n+1))/y(n); 
end



