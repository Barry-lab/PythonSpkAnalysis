function generate_fixtures(trialPath, tet, cell, outFile)
% GENERATE_FIXTURES  Dump golden-master fixtures from the CURRENT MATLAB code.
%
% Run this in MATLAB with Universal_Matlab on the path to produce a reference
% .mat that the Python regression suite (barrelab) validates against. Unlike the
% historic example .mat, this captures the *current* behaviour of read_DACQ.m,
% postprocess_pos_data.m and make_smooth_ratemap.m, so it can bit-validate the
% head-position weighting and the smoothing (see FINDINGS.md).
%
% USAGE
%   generate_fixtures('/data/r1838_10-12-14/83812141', 2, 3, 'fixture_t2c3.mat')
%
% ARGS
%   trialPath : path to the trial basename (no extension)
%   tet, cell : tetrode number and cluster id to build a rate map for
%   outFile   : output .mat path
%
% The Python side loads outFile and asserts np.allclose against each field.

if nargin < 4, outFile = sprintf('fixture_t%dc%d.mat', tet, cell); end

vars = default_vars;
[data, cut] = read_DACQ(trialPath, vars, 0);   % 0 = don't load waveforms

% --- position ---
fix.xy        = data.pos.xy;
fix.direction = data.pos.dir(:);
fix.speed     = data.pos.speed(:);
fix.ppm       = key_value('pixels_per_metre', data.pos.header, 'num');
fix.sample_rate = key_value('sample_rate', data.pos.header, 'num');

% --- spikes / cut ---
exactCut       = cut{tet};
fix.exactCut   = exactCut;
fix.pos_sample = data.tetrode(tet).pos_sample;
ind            = exactCut == cell;
ind            = ind(1:length(data.tetrode(tet).pos_sample));
fix.spikePosSamp = data.tetrode(tet).pos_sample(ind);

% --- binned maps (unsmoothed) + smoothed rate map ---
binPix        = (fix.ppm/100) * vars.rm.binSizePosCm;   % cm bins -> pixels
fix.binPix    = binPix;
fix.occupancy = bin_data('dwell',  'position', binPix, data.pos);
fix.spikeMap  = bin_data('spikes', 'position', binPix, data.pos, fix.spikePosSamp);
fix.smthRm    = make_smooth_ratemap(fix.occupancy, fix.spikeMap, ...
                                    vars.rm.smthKernPos, 'boxcar');

% --- Skaggs spatial information (on adaptively-binned map) ---
[adRate, ~, adPos] = skaggsAdaptiveBinning(fix.occupancy, fix.spikeMap, ...
                                           vars.rm.adaptiveSmoothAlpha);
[fix.bits_per_spike, fix.bits_per_sec] = skaggsInfo(adRate, adPos);

fix.trialPath = trialPath; fix.tet = tet; fix.cell = cell;
save(outFile, '-struct', 'fix');
fprintf('Wrote %s\n', outFile);
end
