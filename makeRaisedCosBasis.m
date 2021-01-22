function [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(BasisPrs)
% [cosBasis, tgrid, basisCtrs] = makeRaisedCosBasis(BasisPrs)
%
% Make basis of raised cosines with a logarithmic or linear time axis
%
% INPUT:  
% ------
% BasisPrs = struct with fields:
% 
%  .nBasis [1]       - number of basis vectors
%  .peakRange [1 2]  - position of 1st and last cosine peaks
%  .dt [1]           - time bin size of bins representing basis [default=1]
%  .logScaling       - 'log' [default] or 'linear', for scaling of time axis
%  .logOffset        - offset for nonlinear stretching of t axis: 
%                        y = log(t+logOffset) 
%                        (larger -> more nearly linear stretching)
%  .timeRange [1 2]  - time range of basis [OPTIONAL]
%  .initialOnesFlag  - flag that (if set to 1) sets first basis vector to 
%                         all ones prior to the 1st peak. [OPTIONAL]
%                        (useful basis for refractory effects in GLM)
%
% OUTPUT:
% -------
%      Basis [nT x nB] - cosine basis vectors
%     ttgrid [nT x 1]  - time lattice on which basis is defined
%   Bcenters [nB x 1]  - centers of each cosine basis function
%
%
%  Example call:
%  -------------
% Bprs.nBasis = 8;  % number of basis vectors
% Bprs.peakRange = [0, 10]; % location of 1st and last peaks
% Bprs.dt = 0.1; % time bin size
% Bprs.logScaling = 'log';  % specify log scaling of time axis
% Bprs.logOffset = 1.5;  % nonlinear stretch factor (larger => more linear)
%
% [BB, tt] = makeRaisedCosBasis(Bprs); 

% ===================== Process inputs ==============

% ---  Extract mandatory fields ----
nB = BasisPrs.nBasis; % number of basis functions
peakRange = BasisPrs.peakRange; % location of first and last basis funcs

% --- Process optional fields ------
try dt = BasisPrs.dt;
catch, dt = 1;
end
try logOffset = BasisPrs.logOffset; 
catch, logOffset = 1;
end
try logScalingStr = BasisPrs.logScaling;
    if strcmpi(logScalingStr,'log')
        logScaling = 1;
    elseif strcmpi(logScalingStr,'linear')
        logScaling = 0;
    else
        warning('Unrecognized logScaling field (valid choices: ''log'' or ''linear''); Defaulting to ''log''');
        logScaling = 1;
    end
catch, logScaling = 1;
end
try timeRange = BasisPrs.timeRange;
    minT = timeRange(1);
    maxT = timeRange(2);
catch, timeRange = [];
end
try initialOnesFlag = BasisPrs.initialOnesFlag;
catch, initialOnesFlag = 0;
end

% Define function for single raised cosine basis function
raisedCosFun = @(x,ctr,dCtr)((cos(max(-pi,min(pi,(x-ctr)*pi/dCtr/2)))+1)/2);

if logScaling
   % ------------- Use log scaling of x axis --------

   % Check that logOffset is positive
    if logOffset <= 0, error('logOffset must be > 0');
    end
    
    % Define nonlinear time axis stretching function and its inverse
    nlin = @(x)(log(x+1e-20));
    invnl = @(x)(exp(x)-1e-20);
    
    % Compute location for cosine basis centers
    logPeakRange = nlin(peakRange+logOffset);   % 1t and last cosine peaks in stretched coordinates
    dCtr = diff(logPeakRange)/(nB-1);   % spacing between raised cosine peaks
    Bctrs = logPeakRange(1):dCtr:logPeakRange(2);  % peaks for cosine basis vectors
    basisPeaks = invnl(Bctrs);  % peaks for cosine basis vectors in rescaled time
    
    % Compute time grid points
    if isempty(timeRange)
        minT = 0; % minimum time bin (where first basis function starts)
        maxT = invnl(logPeakRange(2)+2*dCtr)-logOffset; % maximum time bin (where last basis function stops)
    end
    tgrid = (minT:dt:maxT)'; % time grid
    nT = length(tgrid);   % number of time points in basis
    
    % Make the basis
    cosBasis = raisedCosFun(repmat(nlin(tgrid+logOffset), 1, nB), repmat(Bctrs, nT, 1), dCtr);
    
else
   % ------------- Use linear scaling of x axis --------

   % Compute location for cosine basis centers
    dCtr = diff(peakRange)/(nB-1);  % spacing between raised cosine peaks
    Bctrs = peakRange(1):dCtr:peakRange(2);  % peaks for cosine basis vectors
    basisPeaks = Bctrs;  % vector of raised cosine centers
    
    if isempty(timeRange)
        minT = peakRange(1)-2*dCtr; % min time bin (where 1st basis vector starts)
        maxT = peakRange(2)+2*dCtr; % max time bin (where last basis vector stops)
    end
    tgrid = (minT:dt:maxT)'; % time grid
    nT = length(tgrid);   % number of time points in basis
   
    % Make the basis
    cosBasis = raisedCosFun(repmat(tgrid, 1, nB), repmat(Bctrs, nT, 1), dCtr);
     
end

% If necessary, set first basis vector to 1 before first peak
if initialOnesFlag == 1
    ii = tgrid<=peakRange(1);  % indices to set to 1
    cosBasis(ii,1) = 1; % set first basis vector to constant before peak
end

% Check condition number
cc = cond(cosBasis);  
if cc > 1e12
    warning(sprintf('Raised cosine basis is poorly conditioned (cond # = %.2f)',cc));
end

