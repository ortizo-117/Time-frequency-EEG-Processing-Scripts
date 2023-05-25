%% SFN extracting data from TF spectra
clear all
close all


%% Parameters
n = 3;
n_intensities = 4;
n_chan = 31;
n_times = 400;
n_lf_frx = 8;
n_lg_frx = 8;
n_hg_frx = 12;

tw_lf = [100 300];
tw_lg = [100 200];
tw_hg = [100 200];

subs_names = {'P01' , 'S01' , 'S02'};
intensities = {'I1' 'I2' 'I3' 'I4'};
tf_bands = {'lf' 'lg' 'hg'}


%% Importing ratings
file_loc = 'Z:\30_Oscar_Ortiz\GluEEG\ratings_long.xlsx';

raitings = readtable(file_loc);

%% Importing data
%S01
path_data = 
filename_EEG1 = 
[EEGdata1, chan_names1 ] = readEEGfile([path filenameEEG1]);
