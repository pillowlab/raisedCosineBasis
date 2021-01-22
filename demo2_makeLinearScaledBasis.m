% Demo script 2: 
%
% Generate raised-cosine basis with linear scaling of time axis

% Set up parameters governing basis
Bprs = struct;  % make new struct 'Bprs'
Bprs.nBasis = 8;  % number of basis vectors
Bprs.peakRange = [1, 80]; % location of 1st and last cosine peaks
Bprs.timeRange = [1, 100]; % 1st and last time bins for basis
Bprs.dt = 1; % time bin size
Bprs.logScaling = 'linear';  % 'log' or 'linear'

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

