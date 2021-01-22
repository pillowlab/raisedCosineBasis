% Demo script 1: 
%
% Generate a raised-cosine basis with log scaling of x axis

% Set up parameters governing basis
Bprs = struct;  % make new struct 'Bprs'
Bprs.nBasis = 8;  % number of basis vectors
Bprs.peakRange = [0, 10]; % location of 1st and last cosine peaks
Bprs.dt = 0.1; % time bin size
Bprs.logScaling = 'log';  % 'log' or 'linear'
Bprs.logOffset = 1.25;  % nonlinear stretch factor (larger => more linear)

% Make basis
[cosBasis,tt] = makeRaisedCosBasis(Bprs);

% ---- make plots of basis and summed basis --------
subplot(211);
plot(tt,cosBasis, 'linewidth', 2);
set(gca,'tickdir', 'out');
xlabel('time lag');  ylabel('basis function value');
title('basis');

subplot(212);  % summed basis (to show tiling behavior: sums to constant)
plot(tt,sum(cosBasis,2), 'linewidth', 2);
set(gca,'tickdir', 'out');
xlabel('time lag');  
ylabel('sum');
title('summed basis vectors');


