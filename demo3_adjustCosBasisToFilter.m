% Demo 3:  generate raised cosine basis that attempts to match some known filter

% 1. Make test filter
% -------------------
dt = 0.01; % time bin size
tgrid = (0:dt:5)'; % time grid for filter
kfilt = exp(-0.5*tgrid).*sin(3*tgrid.^0.5); % test filter

% 2. Set up params governing basis 
% ------------------------------------------------------------------

% (Adjust 'nBasis', 'peakRange', and 'logOffset' if desired)

Bprs = struct;  % clear Bprs struct
Bprs.nBasis = 5;  % number of basis vectors
Bprs.peakRange = [0.1, 3]; % location of 1st and last cosine peaks
Bprs.logOffset = .1;  % nonlinear stretch factor (larger => more linear)
Bprs.dt = dt; % time bin size
Bprs.logScaling = 'log';  % 'log' or 'linear'
Bprs.timeRange = [0,tgrid(end)]; % set time range

% Make the basis
[cosBasis,tt] = makeRaisedCosBasis(Bprs);
plot(tt,cosBasis);

% 3. Fit filter in the basis using least-squares
% -----------------------------------------------
kfilt_intrp = interp1(tgrid,kfilt,tt,'spline'); % interpolate filter to basis time points
wts  = (cosBasis'*cosBasis)\(cosBasis'*kfilt_intrp);
khat = cosBasis*wts;

% 4. Make plot
% -------------
subplot(211);  % plot basis
plot(tt,cosBasis, 'linewidth', 2); axis tight;
title('raised cosine basis'); box off;

subplot(212); % plot it
plot(tgrid,kfilt,tt,khat, '--', 'linewidth', 2); hold on;
plot(tgrid,tgrid*0, 'k'); hold off; axis tight;
xlabel('time lag (s)'); 
legend('true filter', 'basis fit');
title('true filter and basis reconstruction');
box off;

