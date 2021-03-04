function MRS_struct = bd_PhilipsReadRaw(fname, spar_fname)
%   Reads Philips DATA/LIST files into Gannet.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-02)
%       goeltzs1@jhmi.edu
%
%   Credits:
%       This code uses the function
%       loadRawKspace.m
%       from the excellent "Matlab raw kspace tools" toolbox
%       (Wouter Potters, Academic Medical Center, Amsterdam, NL)
%       https://bitbucket.org/wpotters/matlab-raw-kspace-tools
%
%   History:
%       2018-03-02: First version.
%       2018-03-26: Fixed bug - phase correction is NOT performed if coil
%                   combination is based on SENSE reconstruction (e.g.
%                   MEGA-PRIAM)

% Extract information from SPAR files that is not in DATA/LIST

ii = 1;
MRS_struct.filename = fname;

% Open spar file and read parameters
sparname   = fopen(spar_fname,'r');
sparheader = textscan(sparname, '%s');
sparheader = sparheader{1};
sparidx    = find(ismember(sparheader, 'samples') == 1); % number of data points
MRS_struct.p.npoints(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'rows') == 1); % number of rows
MRS_struct.p.nrows(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'averages') == 1); % number of averages
MRS_struct.p.Navg(ii) = MRS_struct.p.nrows(ii) * str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'repetition_time') == 1); % TR
MRS_struct.p.TR(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'echo_time') == 1); % TE
MRS_struct.p.TE(ii) = str2double(sparheader{sparidx+2}); %
sparidx = find(ismember(sparheader, 'synthesizer_frequency') == 1); % F0
MRS_struct.p.LarmorFreq(ii) = str2double(sparheader{sparidx+2})/1e6;
sparidx = find(ismember(sparheader, 'sample_frequency') == 1); % readout bandwidth
MRS_struct.p.sw(ii) = str2double(sparheader{sparidx+2});

MRS_struct.p.dt(ii)      = 1/MRS_struct.p.sw(ii);
MRS_struct.p.specRes(ii) = MRS_struct.p.sw(ii) / MRS_struct.p.npoints(ii);
MRS_struct.p.Tacq(ii)    = MRS_struct.p.npoints(ii) * MRS_struct.p.dt(ii);

% Read voxel geometry information.
% THIS IS IN THE ORDER LR-AP-FH!
sparidx = find(ismember(sparheader, 'lr_size') == 1); % voxel size
MRS_struct.p.voxdim(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'ap_size') == 1);
MRS_struct.p.voxdim(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_size') == 1);
MRS_struct.p.voxdim(ii,3) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_off_center') == 1); % voxel center offset
MRS_struct.p.voxoff(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'ap_off_center') == 1);
MRS_struct.p.voxoff(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_off_center') == 1);
MRS_struct.p.voxoff(ii,3) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_angulation') == 1); % voxel angulation (radians)
MRS_struct.p.voxang(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'ap_angulation') == 1);
MRS_struct.p.voxang(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_angulation') == 1);
MRS_struct.p.voxang(ii,3) = str2double(sparheader{sparidx+2});
fclose(sparname);

% Load and process DATA/LIST

data = bd_loadRawKspace(fname);

% Determine number of channels
MRS_struct.p.nreceivers(ii) = length(unique(data.chan));

% Determine number of mixes
% n_mixes = length(unique(data.mix));

% Determine number of averages per mix
n_averages = data.kspace_properties.number_of_signal_averages;
n_av_edit  = n_averages(1);
% if n_mixes == 2
%     n_av_water = n_averages(2); % if second mix exists, it is water
% end
% This will be the NSA as specified on the exam card for the water-suppressed mix (mix = 0)
% and the NSA as specified on the exam card for the water-suppressed mix (mix = 1);

% Determine number of dynamics per NSA. It is not stored in the dynamics
% field, but rather in extra attributes. Could be different for different
% software versions, need to check!
n_dyns      = data.kspace_properties.number_of_extra_attribute_1_values;
n_dyns_edit = n_dyns(1);
%if n_mixes == 2
%    n_dyns_water = n_dyns(2); % if second mix exists, it is water
%end
% Since dynamics are set globally on the exam card, this will be the same
% for both - it will only be the true value for the water-suppressed mix

% Determine number of data points per scan
n_points = data.kspace_properties.F_resolution(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start splitting the list of total scans into its parts:
% Noise scans, water-suppressed scans, and water-unsuppressed scans.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Noise scans have the data type 'NOI'
% isnoise = strcmp(data.typ,'NOI');
% fids_noise = cell2mat(data.complexdata(isnoise));
% n_av_noise = size(fids_noise,1) ./ MRS_struct.p.nreceivers(ii);

% Separate the water-suppressed from the water-unsuppressed scans.
% Water-suppressed scans have the data type 'STD' and mix index 0:
isedit = strcmp(data.typ,'STD') & (data.mix == 0);
fids   = cell2mat(data.complexdata(isedit));
fids   = reshape(fids, [MRS_struct.p.nreceivers(ii) n_av_edit * n_dyns_edit n_points]);
fids   = permute(fids, [1 3 2]);

% Perform coil combination
[~,ind]         = max(abs(mean(fids,3)),[],2);
ind             = mode(ind);
maxPoint        = conj(fids(:,ind,:));
channels_scale  = squeeze(sqrt(sum(maxPoint .* conj(maxPoint),1)));
channels_scale  = repmat(channels_scale, [1 MRS_struct.p.nreceivers(ii) MRS_struct.p.npoints(ii)]);
channels_scale  = permute(channels_scale, [2 3 1]);
maxPoint        = repmat(maxPoint, [1 MRS_struct.p.npoints(ii) 1]) ./ channels_scale;
fids            = fids .* maxPoint;
MRS_struct.fids = conj(squeeze(sum(fids,1)));

end



