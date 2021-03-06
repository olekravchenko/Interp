%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Solver of \nabla^2 = f(x,y) equation 
%                   in rectangle [a,b]x[c,d]
%                   by atomic function \fup_2(x,y) = \fup_2(x)\fup_2(y)
% 
%                    coded by Oleg Kravchenko, Vasily Bondarenko 2016.04.30
%                   UPD1: 2016.04.30
%                   UPD2: 2016.05.01 2AM (BVV)
%                   UPD3: 2016.05.10 3PM (BVV) (Added Eq. Example, Added TODO)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] A. Dabagyan, and E. Fedotova, "Algorithm of interpolation of function
%       of two variables with help of atomic functions,"
%       M.M.A.D.S., No. 1, 38-44, 2003.
%       link: http://fizmathim.com/matematicheskoe-obespechenie-evm-dlya-interpolyatsii-i-approksimatsii-resheniy-kraevyh-zadach-matematicheskoy-fiziki-s-po
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   TODO
% 1. Merge ToDo & boundary former from solve_0bs
% 
% 3. Fix boundary region problem.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%addpath 'blktridiag'

% rectangular area
a = 0; b = 1;
c = 0; d = 1;
% number of gird points
n = 16; m = n;
hx = (b - a) / (n+1);
b=b-hx;
hy = (d - c) / (m+1);
d=d-hy;
% gtid
% x = a:hx:b;
% y = c:hy:d
i = 0:n+1;
j = 0:m+1;
x = a + hx*i;
y = c + hy*j;
% id of the function
id = '15';                             % id numbers '0'..'11'

[xx, yy] = meshgrid(x,y);
[f, f_info] = func(id,xx,yy);


n = n + 1;
m = m + 1;
% M matrix

%As          =  blktridiag(1, 5/26, 5/26, n+3);
%As(1,2)     = -2; 
%As(1,3)     =  1;
%As(n+3,m+2) = -2; 
%As(n+3,m+1) =  1;

a1 = ((1/hx)^2 + (1/hy)^2)*(288/13/16)^2;
a2 = ((1/hx)^2 - 2*(1/hy)^2)*(288/13/16)^2;
a3 = ((1/hx)^2 + (1/hy)^2)*(-4)*(288/13/16)^2;

b1 = 5/26;
b2 = 1;

c1 =  288/13/16;
c2 = -576/13/16;

A          = blktridiag(a3, a2, a2, n+3);
A(1,1)     = b1; 
A(1,2)     = b2; 
A(1,3)     = b1;
A(n+3,m+2) = b1; 
A(n+3,m+1) = b2;
A(n+3,m)   = b1;

A = full(A);

%B = A;
B          = blktridiag(a2, a1, a1, n+3);
B(1,1)     = b1;
B(1,2)     = b2;
B(1,3)     = b1;
B(n+3,m+2) = b1;
B(n+3,m+1) = b2;
B(n+3,m)   = b1;

B = full(B);


As         = blktridiag(b2, b1, b1, n+3);
As(1,1)    = c1;
As(1,2)    = c2;
As(1,3)    = c1;
As(n+3,m+3)= c1;
As(n+3,m+2)= c2;
As(n+3,m)  = c1;

As = full(As);

%Bs  = b1 * As;
Bs         = blktridiag(b2, b1, b1, n+3);
Bs(1,1)    = c1;
Bs(1,2)    = c2;
Bs(1,3)    = c1;

Bs(n+3,m+3)= c1;
Bs(n+3,m+2)= c2;
Bs(n+3,m)  = c1;

%B   = As;
%A   = (288/13) * As;

M2 = blktridiag(B, A, A, n+3);

M2(1:(n+3), 1:m+3)                                      = As;
M2(1:(n+3),(n+3)+1:2*(n+3))                             = Bs;
M2(1:(n+3),2*(n+3)+1:3*(n+3))                           = As;

M2((n+3)*(m+2)+1:(n+3)*(m+3),(n+3)*(m)+1:(n+3)*(m+1))   = As;
M2((n+3)*(m+2)+1:(n+3)*(m+3),(n+3)*(m+1)+1:(n+3)*(m+2)) = Bs;
M2((n+3)*(m+2)+1:(n+3)*(m+3),(n+3)*(m+2)+1:(n+3)*(m+3)) = As;

M = full(M2);

%M = sparse(M2);


% tic
% M = sparse(M);
% toc

% f vector

ff = zeros((n+3),(m+3));
% ff(end,2:end-1) = 2.0402;
% ff(end,end) = 4.0804;
% ff(end,1) = 4.0804;
% ff(2:end-1,end) = 2.0402;
% ff(1,end) = 4.0804;
% ff(end,end) = 4.0804;

ff(2:n+2,2:n+2) = f;
ff = reshape(ff,(n+3)*(m+3),1);
%[L,U] = lu(M);
coeff = M \ ff;
%coeff = bicg(M,ff);
%coeff = U\(L\ff);
coeff = reshape(coeff,(n+3),(m+3));


n1 = n-1;
m1 = n1;
hx = (b - a) / n1;
hy = (d - c) / m1;
i = 0:n1+1;
j = 0:m1+1;
x = a + hx*i;
y = c + hy*j;
s = zeros(n1+2,m1+2);
for i = -1:n+1
 for j = -1:m+1
     s = s + coeff(i+2,j+2) * fup2_small((x-a)/hx - i)' * fup2_small((y-c)/hy - j);
 end
end



% figure(1)
subplot(2,2,1)

surf(xx,yy,f);
shading interp
title(f_info,'Interpreter','tex', 'FontSize',12)

% figure(2)
subplot(2,2,2);
contourf(xx,yy,s,'.')
title('Interpolated function')
box off

subplot(2,2,3)
stem3(xx,yy,s,'Marker','.','LineStyle','none')
surf(xx,yy,s)
shading interp
title('Numeric Solution')

% figure(2)
subplot(224)
surf(xx,yy,s/sin(xx).*sin(yy))
shading interp
title('Error function')


 %Dense Grid
 nxa = n*8; nya = nxa;
 hxa = hx; hya = hy;
 xdence = linspace(a,b,nxa);
 ydence = linspace(c,d,nya);
 sdence = zeros(nxa,nya);
 [xxd,yyd] = meshgrid(xdence,ydence);


 for i = 1:n+3
    for j = 1:m+3
        sdence = sdence + coeff(i,j) * fup2_small((xdence-a)/hxa/(b-a) - i+2)' * fup2_small((ydence-c)/hya/(d-c) - j+2);
    end
 end

figure(2)
%surf(xx,yy,f)
%alpha(0.8)

surf(xxd,yyd,sdence)
title('Numeric Solution')


% er = norm(s(2:n-1,2:m-1) - f(2:n-1,2:m-1),2) / (n*m);
er1 = norm(s - func(id,xx,yy),Inf) / ((n+2)*(m+2));
er2 = norm(s - func(id,xx,yy),2) / ((n+2)*(m+2));
% format long
disp(['Sample #' int2str(n) ', Max-norm: ' num2str(er1)])
disp(['Sample #' int2str(n) ', L2-norm: ' num2str(er2)])

