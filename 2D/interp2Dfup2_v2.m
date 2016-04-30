%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Interpolation of f=f(x,y) function 
%                   in rectangle [a,b]x[c,d]
%                   by atomic function \fup_2(x,y) = \fup_2(x)\fup_2(y)
% 
%                    coded by Oleg Kravchenko 2014.04.22
%                   UPD1: 2014.07.30
%                   UPD2: 2016.04.30
%                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] A. Dabagyan, and E. Fedotova, "Algorithm of interpolation of function
%       of two variables with help of atomic functions,"
%       M.M.A.D.S., No. 1, 38-44, 1977.
%       link: http://fizmathim.com/matematicheskoe-obespechenie-evm-dlya-interpolyatsii-i-approksimatsii-resheniy-kraevyh-zadach-matematicheskoy-fiziki-s-po
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% rectangular area
a = 0; b = 10;
c = 0; d = 10;
% number of gird points
n = 2; m = n;
% n = nx + 3; m = ny + 3;
hx = (b - a) / n;
hy = (d - c) / m;
% gtid
% x = a:hx:b;
% y = c:hy:d
i = 0:n+1;
j = 0:m+1;
x = a + hx*i;
y = c + hy*j;
% id of the function
id = '11';                             % id numbers 1..9

[xx, yy] = meshgrid(x,y);
[f, f_info] = func(id,xx,yy);


n = n + 1;
m = m + 1;
% M matrix
M = zeros((n+3)*(m+3),(n+3)*(m+3));
As = diag(ones(1,n+3),0) + diag((5/26)*ones(1,n+2),-1) + diag((5/26)*ones(1,n+2),1);
As(1,2) = -2; As(1,3) = 1;
As(n+3,m+2) = -2; As(n+3,m+1) = 1;
Bs = -2 * As;
B = As;
A = (5/26) * As;
% M matrix
for k = 1:n+1
    M(k*(n+3)+1:(k+1)*(n+3), (k-1)*(m+3)+1:k*(m+3)) = A;
    M(k*(n+3)+1:(k+1)*(n+3), k*(m+3)+1:(k+1)*(m+3)) = B;
    M(k*(n+3)+1:(k+1)*(n+3), (k+1)*(m+3)+1:(k+2)*(m+3)) = A;
end

M(1:(n+3), 1:m+3) = As;
M(1:(n+3),(n+3)+1:2*(n+3)) = Bs;
M(1:(n+3),2*(n+3)+1:3*(n+3)) = As;
M((n+3)*(m+2)+1:(n+3)*(m+3),(n+3)*(m)+1:(n+3)*(m+1)) = As;
M((n+3)*(m+2)+1:(n+3)*(m+3),(n+3)*(m+1)+1:(n+3)*(m+2)) = Bs;
M((n+3)*(m+2)+1:(n+3)*(m+3),(n+3)*(m+2)+1:(n+3)*(m+3)) = As;

tic
M = sparse(M);
toc

% f vector

ff = zeros((n+3), (m+3));
% ff(end,2:end-1) = 2.0402;
% ff(end,end) = 4.0804;
% ff(end,1) = 4.0804;
% ff(2:end-1,end) = 2.0402;
% ff(1,end) = 4.0804;
% ff(end,end) = 4.0804;

ff(2:n+2,2:n+2) = f;
ff = reshape(ff,(n+3)*(m+3),1);

% coeff = M \ ff;
coeff = cgs(M, ff);
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
subplot(131)
surf(xx,yy,f);
% zlim([min(min(f)) max(max(f))])
shading interp
lighting phong
title(f_info,'Interpreter','latex', 'FontSize',12)
% title('Original function')
% colormap jet

% figure(2)
% subplot(222);
% contourf(xx,yy,s,'.')
% % colormap cool
% % zlim([min(min(s)) max(max(s))])
% % shading interp
% % lighting phong
% title('Interpolated function')
% box off

subplot(132)
surf(xx,yy,f)
alpha(0.8)
hold on
stem3(xx,yy,s,'Marker','.','LineStyle','none')
% zlim([min(min(s)) max(max(s))])
shading interp
lighting phong
alpha(.4)
title('Interpolated function')


subplot(133)
surf(xx,yy,s-func(id,xx,yy))
shading interp
lighting phong
alpha(.4)
title('Error function')
axis([0 1 0 1 min(min(s-func(id,xx,yy))) max(max(s-func(id,xx,yy)))]);

% er = norm(s(2:n-1,2:m-1) - f(2:n-1,2:m-1),2) / (n*m);
er = norm(s - func(id,xx,yy),2) / (n*m);
% format long
disp(['Sample #' int2str(n) ', L2-norm: ' num2str(er)])

