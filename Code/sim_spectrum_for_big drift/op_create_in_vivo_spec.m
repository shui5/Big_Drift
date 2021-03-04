function [ out] = op_create_in_vivo_spec(pathBasis,includeMetabs,spectralwidth,npoints)

%% 1. Load basis set
% Load the specified basis set or the user basis set file
basisSet = load(pathBasis);
basisSet = basisSet.BASIS;
    
% Generate the list of basis functions that are supposed to be included in
% the basis set

% To do: Interface with interactive user input
metabList = fit_createMetabList(includeMetabs);
% No MMs to be included for now
fitMM = 0;
% Create the modified basis set
basisSet = fit_selectMetabs(basisSet, metabList, fitMM);


%% 2. Apply typical line broadening etc.
% ... settings:
RangePPM         = [basisSet.ppm(1), basisSet.ppm(end)];

% ... fit parameters
nMets       = length(includeMetabs);
nMM         = 0;
nBasisFcts  = nMets + nMM; % number of basis functions

load('Lookup_LCM_PRESS_Philips.mat'); % Load Mean values from the Big PRESS Philips dataset

%Remove unwanted entries
delete = ones(length(names),1);
for kk = 1 : length(includeMetabs)
    idx = find(strcmp(names, includeMetabs{kk}));
    if ~isempty(idx)
        delete(idx) = 0;
    end
end


lorentzLB(logical(delete)) = [];
freqShift(logical(delete)) = [];
ampl(logical(delete)) = [];
ph0         = ph0 * pi/180; % zero-order phase correction [convert from deg to rad]
ph1         = ph1 * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]

% Or set the values yourself
% refFWHM     = 12; %FWHM of the spectrum used for the convolution 
ph0         = 0 * pi/180; % zero-order phase correction [convert from deg to rad]
ph1         = 0 * pi/180; % first-order phase correction [convert from deg/ppm to rad/ppm]
% gaussLB     = 11.01; % Gaussian damping [Hz]
% lorentzLB   = fitParams.lorentzLB; % Lorentzian damping [Hz] for each basis function
% freqShift   = zeros(length(metabList),1); % Frequency shift [Hz] for each basis function
% ampl        = fitParams.ampl; % Amplitudes for metabolite/MM/lipid basis functions
% linShape    = lineshape; % Convolution lineshape

% Normalize
if exist('lineShape')
    lineShape = lineShape/sum(lineShape);
end



%%% 2. APPLY THE NON-LINEAR PARAMETERS %%%
% Run the time-domain operations on the metabolite basis functions
% (frequency shift, Lorentzian dampening, Gaussian dampening, zero phase shift)
t = basisSet.t;
for ii=1:nBasisFcts
    basisSet.fids(:,ii) = basisSet.fids(:,ii) .* exp(-1i*freqShift(ii).*t)' .* exp(-lorentzLB(ii).*t)' .* exp(-gaussLB.*t.*t)';    
    basisSet.fids(:,ii) = basisSet.fids(:,ii) * exp(1i*ph0);
end
basisSet.specs = fftshift(fft(basisSet.fids,[],1),1);

% Run the frequency-domain operations on the basis functions
% (first order phase correction)
% Cut out the frequency range of the basis set
basisSet = op_freqrange(basisSet,RangePPM(1),RangePPM(end));
% Create a ppm vector around a pivot point (water)
ppm_ax = basisSet.ppm';
pivotPoint = 4.68;
multiplier = ppm_ax - pivotPoint;
% Apply the linear phase correction
for ii=1:nBasisFcts
    basisSet.specs(:,ii) = basisSet.specs(:,ii) .* exp(1i*ph1*multiplier);
end
basisSet.fids = ifft(fftshift(basisSet.specs,1),[],1);

%%% 3. APPLY THE LINEAR PARAMETERS %%%
% Convolve the lineshape with the metabolite basis functions only
% (NOT the macromolecules or lipids or baseline splines).
if exist('lineShape')
    for kk = 1:basisSet.nMets
        basisSet.specs(:,kk) = conv(basisSet.specs(:,kk), lineShape, 'same');
    end
end

% Calculate the final spectrum
out.specs = basisSet.specs * ampl;
out.fids = ifft(fftshift(out.specs,1),[],1);
out.ppm   = ppm_ax;
out.t = t;
out.sw = spectralwidth;
out.sz = size(out.specs);
out.Bo = basisSet.Bo;
out.te = basisSet.te;
out.centerFreq = basisSet.centerFreq;

%% 3. Resample basis set
% Determine the ppm ranges of both the data and the basis functions.
dwelltime = 1/spectralwidth; % dwelltime [s]
% Calculate t and ppm arrays using the calculated parameters:
f = [(-spectralwidth/2)+(spectralwidth/(2*npoints)):spectralwidth/(npoints):(spectralwidth/2)-(spectralwidth/(2*npoints))];
ppmRangeData = f/(out.Bo*42.577);
% Philips data assumes the center frequency to be 4.68 ppm:
centerFreq = 4.68;
ppmRangeData=ppmRangeData + centerFreq;

ppmRangeData        = ppmRangeData';
ppmRangeBasis       = out.ppm;
ppmIsInDataRange    = (ppmRangeBasis < ppmRangeData(1)) & (ppmRangeBasis > ppmRangeData(end));
if sum(ppmIsInDataRange) == 0
    ppmIsInDataRange    = (ppmRangeBasis > ppmRangeData(1)) & (ppmRangeBasis < ppmRangeData(end));
end

% Now resample the basis functions to match the resolution and frequency
% range (ppm axis) of the data.
fids_interp     = zeros(length(ppmRangeData), 1); % allocate memory
specs_interp    = zeros(length(ppmRangeData), 1); % allocate memory

specs_interp(:,1)      = interp1(ppmRangeBasis(ppmIsInDataRange), out.specs(ppmIsInDataRange,1), ppmRangeData, 'pchip', 'extrap');
%convert back to time domain
%if the length of Fids is odd, then you have to do a circshift of one to
%make sure that you don't introduce a small frequency shift into the fids
%vector.
if mod(length(out.specs),2)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fids_interp(:,1)   = ifft(fftshift(specs_interp(:,1),1),[],1);
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fids_interp(:,1)   = ifft(circshift(fftshift(specs_interp(:,1),1)),[],1);
end


% Create the output resampled basis set container
out.ppm     = ppmRangeData;
out.specs   = specs_interp;
out.fids    = fids_interp;
out.sz      = size(out.fids);
out.n       = length(out.fids);
out.dwelltime = dwelltime;
% Calculate the new spectral width and dwelltime:
dppm                        = abs(out.ppm(2)-out.ppm(1));
ppmrange                    = abs((out.ppm(end)-out.ppm(1)))+dppm;
out.sw   = ppmrange*out.Bo*42.577;
out.dwelltime       = 1/out.sw;
% Calculate the time scale
out.t = (0:out.dwelltime:(out.sz(1)-1)*out.dwelltime);


end