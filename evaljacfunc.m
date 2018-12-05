function [fx,jacx] = evaljacfunc(grad,jac,x)
%func Summary of this function goes here
%   Detailed explanation goes here
% f1 = @(b)b(1)^2-2*b(1)-b(2)+1;
% f2 = @(b)b(1)^2 + b(2)^2-1;
% fx = [f1(x);f2(x)];
% syms x1 x2;
% jac = jacobian([x1^2-2*x1-x2+1,x1^2 + x2^2-1], [x1,x2]);
% x = x.';
% jaCx = double(subs(jac, [x1,x2], x));
syms x1 x2;
x = x.';
fx1 = subs(grad(1), [x1,x2],x);
fx2 = subs(grad(2), [x1,x2],x);
fx = [fx1;fx2];
jacx = double(subs(jac, [x1,x2], x));


% f1 = @(b)(b(1)-3)^2+(b(2)-2)^2;
% fx = [f1(x)];
% syms x1 x2;
% jac = jacobian([(x1-3)^2+(x2-2)^2], [x1,x2]);
% x = x.';
% jaCx = double(subs(jac, [x1,x2], x));

%x1 = x(1)
%x2 = x(2)
%jaCx = jac


end

