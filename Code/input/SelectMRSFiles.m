function [metab, water, struc] = SelectMRSFiles(fileExtension)

c1    = 0;
c2    = 0;
metab = {};
water = {};

switch fileExtension
    
    case {'sdat','SDAT'}
        flist = dir('*.sdat');
        if isempty(flist)
            flist = dir('*.SDAT');
        end
        flist = flist(cellfun(@isempty, strfind({flist.name}, '._'))); %#ok<*STRCLFH>
        for ii = 1:length(flist)
            if contains(cellstr(flist(ii).name),'act')
                c1 = c1 + 1;
                metab(c1) = cellstr(flist(ii).name); %#ok<*AGROW>
            elseif contains(cellstr(flist(ii).name),'ref')
                c2 = c2 + 1;
                water(c2) = cellstr(flist(ii).name);
            end
        end
        
        
    case 'dat'
        flist = dir('*.dat');
        flist = flist(cellfun(@isempty, strfind({flist.name}, '._')));
        for ii = 1:length(flist)
            if contains(cellstr(flist(ii).name),'dat')
                c1 = c1 + 1;
                metab(c1) = cellstr(flist(ii).name);
            elseif contains(cellstr(flist(ii).name),'ref')
                c2 = c2 + 1;
                water(c2) = cellstr(flist(ii).name);
            end
        end
        
        
    case '7'
        flist = dir('*.7');
        flist = flist(cellfun(@isempty, strfind({flist.name}, '._')));
        for ii = 1:length(flist)
            metab(ii) = cellstr(flist(ii).name);
        end
        
        
    case 'data'
        flist = dir('*.data');
        flist = flist(cellfun(@isempty, strfind({flist.name}, '._')));
        for ii = 1:length(flist)
            metab(ii) = cellstr(flist(ii).name);
        end
        
        
    case 'rda'
        flist = dir('*.rda');
        flist = flist(cellfun(@isempty, strfind({flist.name}, '._')));
        for ii = 1:length(flist)
            metab(ii) = cellstr(flist(ii).name);
        end
        
    case {'ima','IMA'}
        flist = dir('*.IMA');
        flist = flist(cellfun(@isempty, strfind({flist.name}, '._')));
        metab = cellstr(flist(1).name);
        
    otherwise
        error('Invalid extension entered.')        
        
end

if nargout > 2
    if strcmp(fileExtension,'7')
        flist = dir;
    else
        flist = dir('*.nii');
    end
    flist = flist(~ismember({flist.name}, {'.','..','.DS_Store'}));
    for ii = 1:length(flist)
        struc(ii) = cellstr(flist(ii).name);
    end
end



