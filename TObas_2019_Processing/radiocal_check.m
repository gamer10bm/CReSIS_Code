% script RadioCal Check
%
% The purpose of this script is to check the effect of radiometric
% calibration given a qlook product
%
% Author: Bailey Miller
%

%% Settings
%Which radar?
radar = 'accum'; season = '2018_Antarctica_TObas';
%Define the segment and frame of interest
seg = '20190207_02'; frame = 3; wf = 1;
%Define the rangeline of interest
rangeline = 550;
%Which product?
dataprod = 'qlook';

%% Automated
%Load the specified qlook product
if ispc
    basepath = 'X:/ct_data/';
else
    basepath = '/cresis/snfs1/dataproducts/';
end
prodpath = fullfile(basepath,radar,season,sprintf('CSARP_%s',dataprod),seg);
filepath = sprintf('%s/Data_img_%02d_%s_%03d.mat',prodpath,wf,seg,frame);


param = struct('radar_name','accum','season_name','2018_Antarctica_TObas');
param.day_seg = '20190207_02';
frm = 3;
img = 1;
out_path = 'qlook';
fn_dir = ct_filename_out(param,out_path)
fn = fullfile(fn_dir,sprintf('Data_img_%02d_%s_%03d.mat',img,param.day_seg,frm))


open(filepath);
%Put data into log scale
logData = 10.*log10(Data);
%Load speed of light 
c = physconst('LightSpeed');
%Convert range bins from time to distance
Range = 2*c*Time;
%Convert range into 1/R^2 relationship in dB
R_2 = 1./Range.^2;
R_2dB = 10.*log10(R_2);

%For the given rangeline, find the maximum log value
[peakval, peakind] = max(logData(:,rangeline));
peakval = peakval
%Find the 1/R^2 for that peak
peakR_2dB = R_2dB(peakind)

%Find the max for each range bin
[peaks, peaksind] = max(logData);
peaksR_2dB = R_2dB(peaksind);
%% Plotting
%Plot the echogram
figure(1)
imagesc(logData)
hold on
plot([rangeline rangeline],[0,size(Data,1)],'r','LineWidth',2)
hold off

%Plot the response for the given rangeline
figure(2)
plot(R_2dB,logData(:,rangeline))
xlabel('1/R^2 (dB)')
ylabel('Response (dB)')
grid on

%Plot the peak for each rangeline
figure(3)
plot(peaks)
hold on
plot(peaksR_2dB)
hold off
ylabel('dB')
grid on