% Test script to illustrate generation of raised cosine basis 


%% Example 1: generate and plot an exable basis

% Set up parameters governing basis
Bprs.nh = 8;  % number of basis vectors
Bprs.endpoints = [0, 10]; % location of 1st and last cosines
Bprs.b = 1.5;  % nonlinear stretch factor (larger => more linear)
Bprs.dt = 0.1; % time bin size

% Make basis
[tt, Borth, B] = makeRaisedCosBasis(Bprs);

% Plot basis and summed basis
subplot(211);
plot(tt,B);
xlabel('time lag');  ylabel('coefficient');
subplot(212);  % summed basis (to show tiling behavior: sums to constant)
plot(tt,sum(B,2));
xlabel('time lag');  ylabel('summed basis vectors');


%% Example 2: generate raised cosine basis that attempts to match some known filter

% Make test filter
dt = 0.01;
tgrid = (0:dt:5)';
kfilt = exp(-1*tgrid).*sin(5*tgrid.^0.5);

% Set up / tweak parameters governing basis
Bprs.nh = 8;  % number of basis vectors
Bprs.endpoints = [.05, 2.4]; % location of 1st and last cosines
Bprs.b = .1;  % nonlinear stretch factor (larger => more linear)
Bprs.dt = dt; % time bin size

% Make basis
[tt, Borth, B] = makeRaisedCosBasis(Bprs);
subplot(211);  % plot basis
plot(tt,B); axis tight;
title('raised cosine basis');

% Now fit filter in basis
kfilt_intrp = interp1(tgrid,kfilt,tt,'spline'); % interpolate filter to basis time points
wts  = (B'*B)\(B'*kfilt_intrp);
khat = B*wts;

subplot(212); % plot it
plot(tgrid,kfilt,tt,khat, '--'); axis tight;
xlabel('time lag (s)'); 
legend('true filter', 'basis fit');
title('true filter and basis reconstruction');

