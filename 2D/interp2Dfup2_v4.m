%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Interpolation of f=f(x,y) function 
%                   in rectangle [a,b]x[c,d]
%                   by atomic function \fup_2(x,y) = \fup_2(x)\fup_2(y)
% 
%                    coded by Oleg Kravchenko 2014.04.22
%                   UPD1: 2014.07.30 - Different type of norm's computation
%                                      (not bult-in by norm function)
%                   UPD2: 2016.04.30
%                    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] A. Dabagyan, and E. Fedotova, "Algorithm of interpolation of function
%       of two variables with help of atomic functions,"
%       M.M.A.D.S., No. 1, 38-44, 2003.
%       link: http://fizmathim.com/matematicheskoe-obespechenie-evm-dlya-interpolyatsii-i-approksimatsii-resheniy-kraevyh-zadach-matematicheskoy-fiziki-s-po
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

mypath = 'img\';

%% Initial Settings
% rectangular area
a = 0; 
b = 10;
c = 0; 
d = 10;

% function's id
id = '3';                             % id numbers '1'..'12'


% number of gird points
n   = 11; 
m   = n;
hx  = (b - a) / (n - 1);
hy  = (d - c) / (m - 1);
x   = a:hx:b;
y   = c:hy:d;

[xx, yy]    = meshgrid(x, y);
[f, f_info] = func(id,xx,yy);

% surf(xx, yy, f)     % test plot suft

%% matrix assembling
As          =  blktridiag(1, 5/26, 5/26, n+2);
As(1,2)     = -2; 
As(1,3)     =  1;
As(n+2,m+1) = -2; 
As(n+2,m+0) =  1;

Bs          = -2 * As;
B           = As;
A           = (5/26) * As;

M2                                                      = blktridiag(B, A, A, n+2);
M2(1:(n+2), 1:m+2)                                      = As;
M2(1:(n+2),(n+2)+1:2*(n+2))                             = Bs;
M2(1:(n+2),2*(n+2)+1:3*(n+2))                           = As;
M2((n+2)*(m+1)+1:(n+2)*(m+2),(n+2)*(m-1)+1:(n+2)*(m+0)) = As;
M2((n+2)*(m+1)+1:(n+2)*(m+2),(n+2)*(m+0)+1:(n+2)*(m+1)) = Bs;
M2((n+2)*(m+1)+1:(n+2)*(m+2),(n+2)*(m+1)+1:(n+2)*(m+2)) = As;

M2 = sparse(M2);

ff = zeros(n+2, m+2);
% BC
% ff(1,2:end-1)     = - hy^2 / 50;
% ff(end,2:end-1)   = - hy^2 / 50;
% ff(2:end-1,1)     = + hx^2 / 72;
% ff(2:end-1,end)   = + hx^2 / 72;

% ff(1,2:end-1)     = + hx^2 / 72;
% ff(end,2:end-1)   = + hx^2 / 72;
% ff(2:end-1,1)     = - hy^2 / 50;
% ff(2:end-1,end)   = - hy^2 / 50;

% ff(end,2:end-1) = 2.0402;
% ff(end,1) = 4.0804;
% ff(2:end-1,end) = 2.0402;
% ff(1,end) = 4.0804;
% ff(end,end) = 4.0804;

ff(2:n+1,2:n+1) = f;
ff = reshape(ff,(n+2)*(m+2),1);

% coeff = M2 \ ff;

coeff = cgs(M2, ff);
coeff = reshape(coeff,(n+2),(m+2));

%% reconstruction
ntx = 800;
nty = ntx;
tx  = a:(b-a)/(ntx - 1):b;
ty  = c:(d-c)/(nty - 1):d;
[txx, tyy] = meshgrid(tx, ty);

% ntx = length(tx);
% nty = length(ty);
s   = zeros(ntx, nty);
fs  = func(id, txx, tyy) ;

tic
for i = -1:n
    for j = -1:m
        s = s + coeff(i+2, j+2) * fup2_small((tx-a)/hx - i)' * fup2_small((ty-c)/hy - j);
    end
end
toc;

%% Plotting
plt = 1;
if plt == 1
    hfig = figure(1);
    %fig1
%     subplot(131)
    surf(xx,yy,f);
    offset = 0.05;
    zlim([min(min(f))-offset max(max(f))+offset])
    shading interp, lighting phong
    title(f_info,'Interpreter','latex', 'FontSize',12)
    saveas(gcf, [mypath 'png\' id '_1.png'])
    saveas(gcf, [mypath 'pdf\' id '_1.pdf'])
    saveas(gcf, [mypath 'tiff\' id '_1.tiff'])


    %fig2
%     subplot(132)
    figure(2);
    surf(xx,yy,f)   
    zlim([min(min(s))-offset max(max(s))+offset])
    shading interp, lighting phong
    alpha(.25)
    title('Interpolated function')

    hold on
    stem3(xx,yy,f,'Marker','.','LineStyle','none')
    saveas(gcf, [mypath 'png\' id '_2.png'])
    saveas(gcf, [mypath 'pdf\' id '_2.pdf'])
    saveas(gcf, [mypath 'tiff\' id '_2.tiff'])
    
    %fig3
%     subplot(133)
    figure(3);
    surf(txx,tyy,s-fs)
    alpha(0.8);
    shading interp, lighting phong
    title('Error function')
    saveas(gcf, [mypath 'png\' id '_3.png'])
    saveas(gcf, [mypath 'pdf\' id '_3.pdf'])
	saveas(gcf, [mypath 'tiff\' id '_3.tiff'])
end


%% Errors Disp
er_norm1 = max(max(abs(s - fs)));
er_norm2 = norm(s - fs, 2) / sqrt(ntx*nty);



disp(['Sample #' int2str(n) ', Max-norm: ' num2str(er_norm1)])
disp(['Sample #' int2str(n) ', L2-norm: ' num2str(er_norm2)])

sprintf('%10.8f', er_norm1)
sprintf('%10.8f', er_norm2)