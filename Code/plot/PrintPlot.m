function PrintPlot(outputName)

if nargin < 1
    outputName = '~/Desktop/plot';
end

format      = '-png';
resolution  = '-r300';
padding     = '-p0.01';
crop        = '';
transparent = '';

eval(['export_fig ' outputName ' ' format ' ' resolution ' ' padding ' ' crop ' ' transparent]);