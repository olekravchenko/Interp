%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Solver of \nabla^2 = f(x,y) equation 
%                   in rectangle [a,b]x[c,d]
%                   by B_3 (Cubic B-Spline)
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
% 1. Replace Fixed coefficients with values of particular function (+/-)
% 1.1 Review coefficients at boundary region.
% 2. Fix boundary region problem.
% 3. Add Possibilyty to use not only NULL value at boundary region
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%addpath 'blktridiag'

% rectangular area
a = 0; b = 1;
c = a; d = b;
% number of gird points
n = 20; m = n;
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

%пока что hx==hy,  поэтому несколько пренебрежём стандартными условиями

a1 = -6./hx^2;
a2 = 0.75/hx^2;

b1 = 1;
b2 = 0.25;
b3 = 0.0625;

%c1 =  288/13/16;
%c2 = -576/13/16;

A          = blktridiag(a1, a2, a2, n+3);
A(1,1)     = b2; 
A(1,2)     = b1; 
A(1,3)     = b2;

A(n+3,m+2) = b2; 
A(n+3,m+1) = b1;
A(n+3,m)   = b2;

A = full(A);

%B = A;
B          = blktridiag(a2, a2, a2, n+3);
B(1,1)     = b3;
B(1,2)     = b2;
B(1,3)     = b3;

B(n+3,m+2) = b3;
B(n+3,m+1) = b2;
B(n+3,m)   = b3;

B = full(B);


As         = blktridiag(b1, b2, b2, n+3);
As(1,1)    = a2;
As(1,2)    = a1;
As(1,3)    = a2;

As(n+3,m+3)= a2;
As(n+3,m+2)= a1;
As(n+3,m)  = a2;

As = full(As);

Bs         = blktridiag(b2, b3, b3, n+3);
Bs(1,1)    = a2;
Bs(1,2)    = a2;
Bs(1,3)    = a2;

Bs(n+3,m+3)= a2;
Bs(n+3,m+2)= a2;
Bs(n+3,m)  = a2;

Bs = full(Bs);

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
     s = s + coeff(i+2,j+2) * bspline((x-a)/hx - i)' * bspline((y-c)/hy - j);
 end
end



% figure(1)
subplot(2,2,1)

surf(xx,yy,f);
% zlim([min(min(f)) max(max(f))])
shading interp
%lighting phong
title(f_info,'Interpreter','tex', 'FontSize',12)
% title('Original function')
% colormap jet

% figure(2)
subplot(2,2,2);
contourf(xx,yy,s,'.')
% colormap cool
% zlim([min(min(s)) max(max(s))])
% shading interp
% lighting phong
title('Interpolated function')
box off

subplot(2,2,3)
%surf(xx,yy,f)
%alpha(0.8)
%hold on
stem3(xx,yy,s,'Marker','.','LineStyle','none')
surf(xx,yy,s)
% zlim([min(min(s)) max(max(s))])
shading interp
%lighting phong
%alpha(.4)
title('Numeric Solution')

% figure(2)
subplot(224)
surf(xx,yy,s/sin(xx).*sin(yy))
shading interp
%lighting phong
%alpha(.4)
title('Error function')
%axis([a b c d min(min(s-func(id,xx,yy))) max(max(s-func(id,xx,yy)))]);

 %Dense Grid
 nxa = n*8; nya = nxa;
 hxa = hx; hya = hy;
 xdence = linspace(a,b,nxa);
 ydence = linspace(c,d,nya);
 sdence = zeros(nxa,nya);
 [xxd,yyd] = meshgrid(xdence,ydence);


 for i = 1:n+3
    for j = 1:m+3
        sdence = sdence + coeff(i,j) * bspline((xdence-a)/hxa/(b-a) - i+2)' * bspline((ydence-c)/hya/(d-c) - j+2);
    end
 end

figure(2)
%surf(xx,yy,f)
%alpha(0.8)

surf(xxd,yyd,sdence)
% zlim([min(min(s)) max(max(s))])
%shading interp
%lighting phong
%alpha(.4)
title('Numeric Solution')


% er = norm(s(2:n-1,2:m-1) - f(2:n-1,2:m-1),2) / (n*m);
er1 = norm(s - func(id,xx,yy),Inf) / ((n+2)*(m+2));
er2 = norm(s - func(id,xx,yy),2) / ((n+2)*(m+2));
% format long
disp(['Sample #' int2str(n) ', Max-norm: ' num2str(er1)])
disp(['Sample #' int2str(n) ', L2-norm: ' num2str(er2)])

