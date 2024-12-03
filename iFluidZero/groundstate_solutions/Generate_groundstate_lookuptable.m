clear all; close all; clc;

%% Calculate table of ground state solutions of LL model
N         = 2^10; % resolution of rapidity grid 
M         = 50000; % number of table entries 

lambda    = logspace(-2, 4, M);
gamma     = zeros(1, M);
e         = zeros(1, M);
K         = zeros(1, M); % note, Fermi point K in units of density

tic
for i = 1:M
    [~, gamma(i), e(i), K(i)] = solve_LL_groundstate(lambda(i), N);
end
toc

mu = 3*e - gamma.*gradient(e,gamma);
        
%% Save tables to disk

filename    = 'LL_GS_tables.mat';

save(filename, 'gamma','lambda', 'e', 'K', 'mu')