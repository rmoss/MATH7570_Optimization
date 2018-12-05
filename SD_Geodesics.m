function  [z_approx ,phi ,g ,k] = SD_Geodesics(a ,b ,c ,d ,alpha ,beta ,z)

[phi ,g] = Geodesics(a ,b ,c ,d ,alpha ,beta ,z);

epsilon = 0.00001;
B = 0.1;
tau = 0.5;
k = 0;

while (norm(g) > epsilon)
    p = -g;
    
    A = 1;
    while (Geodesics(a ,b ,c ,d ,alpha ,beta ,z + A * p) > phi + A * B * g' * p)
        A = A * tau;
    end

    z = z + A * p;
    [phi ,g] = Geodesics(a ,b ,c ,d ,alpha ,beta ,z);
    k = k+1;
end

z_approx = z;
N =length(z)/2;
x(1) = a;
y(1) = b;
x(N+2) = c;
y(N+2) = d;
for i = 1:N
    x(i+1) = z(2*i-1);
    y(i+1) = z(2*i);
end
plot(x,y);

