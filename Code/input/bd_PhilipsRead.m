function MRS_struct = bd_PhilipsRead(fname, sparname)
% RE/CJE Parse SPAR file for header info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PhilipsRead is designed to handle 'off-first' data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 1;
MRS_struct.filename = fname;

% Work out data header name
sparname = fopen(sparname,'r');
sparheader = textscan(sparname, '%s');
sparheader = sparheader{1};
sparidx = find(ismember(sparheader, 'samples') == 1);
MRS_struct.p.npoints(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'rows') == 1);
MRS_struct.p.nrows(ii) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'averages') == 1);
MRS_struct.p.Navg(ii) = MRS_struct.p.nrows(ii) * str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'repetition_time') == 1);
MRS_struct.p.TR(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'echo_time') == 1);
MRS_struct.p.TE(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'synthesizer_frequency') == 1);
MRS_struct.p.LarmorFreq(ii) = str2double(sparheader{sparidx+2})/1e6;
sparidx = find(ismember(sparheader, 'sample_frequency') == 1);
MRS_struct.p.sw(ii) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'ap_size') == 1);
MRS_struct.p.voxdim(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_size') == 1);
MRS_struct.p.voxdim(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_size') == 1);
MRS_struct.p.voxdim(ii,3) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'ap_off_center') == 1);
MRS_struct.p.voxoff(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_off_center') == 1);
MRS_struct.p.voxoff(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_off_center') == 1);
MRS_struct.p.voxoff(ii,3) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'ap_angulation') == 1);
MRS_struct.p.voxang(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_angulation') == 1);
MRS_struct.p.voxang(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_angulation') == 1);
MRS_struct.p.voxang(ii,3) = str2double(sparheader{sparidx+2});

MRS_struct.fids = bd_SDATreadMEGA(fname, MRS_struct.p.npoints(ii), MRS_struct.p.nrows(ii));

% Undo phase cycling
corrph = repmat([-1 1], [1 size(MRS_struct.fids,2)/2]);
corrph = repmat(corrph, [size(MRS_struct.fids,1) 1]);
MRS_struct.fids = MRS_struct.fids .* corrph;

% Re-introduce initial phase step
MRS_struct.fids = MRS_struct.fids .* ...
    repmat(conj(MRS_struct.fids(1,:)) ./ abs(MRS_struct.fids(1,:)), [MRS_struct.p.npoints(ii) 1]);

MRS_struct.fids = conj(MRS_struct.fids);

end



