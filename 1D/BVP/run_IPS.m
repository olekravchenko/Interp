%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Numerical solution of Poisson equation by 
%                       IDO and CIB-BS^K methods 
%                       f''(x) = g(x),  x \in [a,b]                  
%                       f(a) = f_1, f(b) = f_2.
%                    coded by Oleg Kravchenko 2016.01.28
%                    UPD1: 01/05/2016: Comments updated
%                    UPD2: 02/05/2016: Unified structure implemented for
%                                       IDO, CIP-BS^0,1 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] T. Aoki, "Interpolated Differential Operator (IDO) scheme for solving
%       partial differential equations," Comp. Phys. Comm., Vol. 102, 
%       132-146, 1997.
%       DOI: http://dx.doi.org/10.1016/S0010-4655(97)00020-9
% [2] T. Utsumi, T. Yabe, J. Koga, T. Aoki, M. Sekine, Y. Ogata, 
%     E. Matsunaga, "A note on the basis set approach in the
%       constrained interpolation profile method," Journal of Comp. Phys., Vol. 196, 
%       1-7, 2004.
%       DOI: http://dx.doi.org/10.1016/j.jcp.2003.10.019
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

addpath('core')

%% Exact solution
% f = @(t) sin(2*pi*t);
% f = @(t) -(1/(2*pi)^2)*sin(2*pi*t);
% f = @(x) x.^3 .* sin(2*pi*x);
% f = @(xi) -(1/pi^2)*sin(pi*xi);
f = @(xi) sin(4*xi);

%% Computation
xL = 0; xR = 1; nx = 8;
method = 'CIPBS1';                          % 'IDO', 'CIPBS0', 'CIPBS1'
tic
[x, y, d1y] = InterpPoissonSolver(method,f,xL,xR,nx);
toc

%% Dence grid
nx = 2^9;
xdence = linspace(xL,xR,nx);
fdence = f(xdence);

%% Draw plot
% Set Default Figure Options
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',13)

figure,
plot(xdence,fdence,x,y,'o-')
xlabel('$x$'); ylabel('$f(x)$');
title(['Numerical solution of 1D Poisson equation by ' method])