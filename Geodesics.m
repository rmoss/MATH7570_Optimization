function [phi ,g] = Geodesics(a ,b ,c ,d ,alpha ,beta ,z)

N =length(z)/2;
phi = 0;

x(1) = a;
y(1) = b;
x(N+2) = c;
y(N+2) = d;
for i = 1:N
    x(i+1) = z(2*i-1);
    y(i+1) = z(2*i);
end

rho = @(m,n)1 + alpha * exp(-beta * (m^2 + n^2));
rhox = @(m,n)-2 * alpha * beta * m * exp(-beta * (m^2 + n^2));
rhoy = @(m,n)-2 * alpha * beta * n * exp(-beta * (m^2 + n^2));

for i = 1:N+1
    phi = phi + (N+1) * rho(x(i),y(i)) * ((x(i+1) - x(i))^2 + (y(i+1) - y(i))^2);
end

for i = 2:N+1
    g(2*i-3,1) = (N+1)*(rho(x(i-1),y(i-1))*2*(x(i)-x(i-1))+rho(x(i),y(i))*2*(x(i)-x(i+1))+rhox(x(i),y(i))*((x(i+1)-x(i))^2+(y(i+1)-y(i))^2));
    g(2*i-2,1) = (N+1)*(rho(x(i-1),y(i-1))*2*(y(i)-y(i-1))+rho(x(i),y(i))*2*(y(i)-y(i+1))+rhoy(x(i),y(i))*((x(i+1)-x(i))^2+(y(i+1)-y(i))^2));
end
