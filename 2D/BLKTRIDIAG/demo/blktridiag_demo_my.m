%% A 10x10 tridiagonal matrix, with 2 on the diagonal, -1 on the off diagonal.

close all;
clear all;
clc;

addpath '..'

A = blktridiag(2,-1,-1,10);

% The sparsity pattern is correct
spy(A)

% and the elements are as designated
full(A)

%% A lower block bidiagonal matrix with replicated blocks

% with 2x2 blocks of ones on the main diagonal, and
% 2x2 blocks of twos on the sub-diagonal

A = blktridiag(ones(3),2*ones(3),zeros(3)+.1,6);

spy(A)
full(A)

%% A block tridiagonal matrix with replicated blocks

Amd     = reshape(1:9,3,3);
Asub    = reshape(11:19,3,3);
Asup    = reshape(21:29,3,3);
A       = blktridiag(Amd,Asub,Asup,4);

spy(A)
full(A)


