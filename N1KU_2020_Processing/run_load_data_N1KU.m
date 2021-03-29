% Script run_load_data
%
% Several examples of using load_data
%
% Author: John Paden
%
% See also load_data, pulse_compress

% =======================================================================
% Setup loading parameters for example 1
%  - Plot data in time-space, time-wavenumber, freq-space, and
%    freq-wavenumber domains
% =======================================================================

param = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'),'20200128_01');

% Determine which records you want to load:
frames = frames_load(param);
frm = 13;
param.load_data.recs = frames.frame_idxs(frm) + 0 + [0 15000];

%   param.load_data.imgs = {[-1j 5]};
%   param.load_data.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 2 9; 2 10; 2 11; 2 12; 2 13; 2 14; 2 15; 2 16]};
param.load_data.imgs                  = {[1 1]};
param.load_data.pulse_comp            = true;
param.load_data.raw_data              = false;
%   param.load_data.ft_wind               = @hanning;
%   param.load_data.combine_rx            = false;

% Load data
[hdr,data] = load_data(param);

% Plot data
img = 1;
wf_adc_idx = 1;
wf = param.load_data.imgs{img}(wf_adc_idx,1);
adc = param.load_data.imgs{img}(wf_adc_idx,2);

% Convert to voltage at ADC input
%data{img}(:,1) = data{img}(:,1) * 10^(param.radar.wfs(wf).adc_gains_dB(adc)/20));

