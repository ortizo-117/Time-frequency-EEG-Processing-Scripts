%% Purpose: Takes Brainviz output for trials at 4 different stim intensities 
%     and bins trials based on percieved intensities.
clear all 
close all

addpath C:\Users\User\Documents\MATLAB\eeglab2021.0
addpath Z:\30_Oscar_Ortiz\GluEEG\5_Matlab\code_gamma
cd Z:\30_Oscar_Ortiz\GluEEG\3_Exports\8_AB
load('chan_info.mat');

eeglab
s_rate = 1000;
time_range = [-400 1000]; % in s
time_vect = linspace(time_range(1),time_range(2),abs(((time_range(1)-time_range(2)))/1000)*s_rate);
n_trials = 20;
n_chans= 32;
chan_name = 'Cz';
pnts = 1400;
n_sets = 4; % number of datasets being loaded in

%% Loading up MRK files
% First
filter = '*.MRK';
my_title = 'Select Marker File 1 ';
[filenameMRK1,path] = uigetfile(filter,my_title);
[extracted_raitings1, n_raitings1] = readMRKfile([path filenameMRK1]);
extracted_raitings1 = cat(2,extracted_raitings1,2);
extracted_raitings1(end) = [];

% Second
my_title = 'Select Marker File 2 ';
[filenameMRK2,path] = uigetfile(filter,my_title);
[extracted_raitings2, n_raitings2] = readMRKfile([path filenameMRK2]);
extracted_raitings2 = cat(2,extracted_raitings2,2);
%extracted_raitings2(end) = [];

% Third
my_title = 'Select Marker File 3 ';
[filenameMRK3,path] = uigetfile(filter,my_title);
[extracted_raitings3, n_raitings3] = readMRKfile([path filenameMRK3]);
extracted_raitings3 = cat(2,extracted_raitings3,2);
%extracted_raitings3(end) = [];


% Fourth
my_title = 'Select Marker File 4 ';
[filenameMRK4,path] = uigetfile(filter,my_title);
[extracted_raitings4, n_raitings4] = readMRKfile([path filenameMRK4]);
extracted_raitings4 = cat(2,extracted_raitings4,2);
extracted_raitings4(end) = [];



%% Loading up EEG data

% First
filter = '*.csv';
my_title = 'Select EEG File 1';
[filenameEEG1,path] = uigetfile(filter,my_title);
[EEGdata1, chan_names1 ] = readEEGfile([path filenameEEG1]);


% Second
filter = '*.csv';
my_title = 'Select EEG File 2';
[filenameEEG2,path] = uigetfile(filter,my_title);
[EEGdata2, chan_names2 ] = readEEGfile([path filenameEEG2]);

% Third
filter = '*.csv';
my_title = 'Select EEG File 3';
[filenameEEG3,path] = uigetfile(filter,my_title);
[EEGdata3, chan_names3 ] = readEEGfile([path filenameEEG3]);

% Fourth
filter = '*.csv';
my_title = 'Select EEG File 4';
[filenameEEG4,path] = uigetfile(filter,my_title);
[EEGdata4, chan_names4 ] = readEEGfile([path filenameEEG4]);


%% Reshape all EEG files and concatenate them all in one
% Double Check n_chans and n_trials - could be that data was cut off or
% extra channels were exported 

% Reshape
EEG1_reshaped = reshape(EEGdata1',n_chans,pnts,n_trials);
EEG2_reshaped = reshape(EEGdata2',n_chans,pnts,n_trials);
EEG3_reshaped = reshape(EEGdata3',n_chans,pnts,n_trials);
EEG4_reshaped = reshape(EEGdata4',n_chans,pnts,n_trials);

%% Concatenate EEG in matrix B.
B = cat(3,EEG1_reshaped,EEG2_reshaped);
B = cat(3,B,EEG3_reshaped);
B = cat(3,B,EEG4_reshaped);


%% Binning

% Concatenate raitings together
R = cat(2,extracted_raitings1,extracted_raitings2);
R = cat(2,R,extracted_raitings3);
R = cat(2,R,extracted_raitings4);

% Bin into four categories based on the data
binned_raitings = bin_raitings_single_chan(R);

% Get locations for trial binning
R1_loc = find(binned_raitings ==1);
R2_loc = find(binned_raitings ==2);
R3_loc = find(binned_raitings ==3);
R4_loc = find(binned_raitings ==4);

%% Indexing EEG based on raitings

EEG_I1 = B(:,:,R1_loc);
EEG_I2 = B(:,:,R2_loc);
EEG_I3 = B(:,:,R3_loc);
EEG_I4 = B(:,:,R4_loc);


%% Storing in structure

% EEG data
struct1.data= EEG_I1;
struct2.data= EEG_I2;
struct3.data= EEG_I3;
struct4.data= EEG_I4;


% sampling rate
struct1.srate = s_rate;
struct2.srate = s_rate;
struct3.srate = s_rate;
struct4.srate = s_rate;


% times
struct1.times = time_vect;
struct2.times = time_vect;
struct3.times = time_vect;
struct4.times = time_vect;


% electrode locations
struct1.chanlocs = chanlocs;
struct2.chanlocs = chanlocs;
struct3.chanlocs = chanlocs;
struct4.chanlocs = chanlocs;

% chaninfo
struct1.chaninfo = chaninfo;
struct2.chaninfo = chaninfo;
struct3.chaninfo = chaninfo;
struct4.chaninfo = chaninfo;

% xmin
struct1.xmin = min(struct1.times)/1000;
struct2.xmin = min(struct2.times)/1000;
struct3.xmin = min(struct3.times)/1000;
struct4.xmin = min(struct4.times)/1000;

% xmax

struct1.xmax = max(struct1.times)/1000;
struct2.xmax = max(struct2.times)/1000;
struct3.xmax = max(struct3.times)/1000;
struct4.xmax = max(struct4.times)/1000;


% pnts
struct1.pnts = length(time_vect);
struct2.pnts = length(time_vect);
struct3.pnts = length(time_vect);
struct4.pnts = length(time_vect);


%% compute time-freq transform for low frequency components.

lf_range = [0.1 30];

figure
[C4_ersp_I1 C4_itc C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I1] = pop_newtimef( struct1, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct1.chanlocs, 'chaninfo', struct1.chaninfo, 'caption', 'C4 - I1', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I1 C3_itc C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I1] = pop_newtimef( struct1, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct1.chanlocs, 'chaninfo', struct1.chaninfo, 'caption', 'C3 - I1', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);

figure
[C4_ersp_I2 C4_itc C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I2] = pop_newtimef( struct2, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct2.chanlocs, 'chaninfo', struct2.chaninfo, 'caption', 'C4 - I2', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I2 C3_itc C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I2] = pop_newtimef( struct2, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct2.chanlocs, 'chaninfo', struct2.chaninfo, 'caption', 'C3 - I2', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);



figure
[C4_ersp_I3 C4_itc C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I3] = pop_newtimef( struct3, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct3.chanlocs, 'chaninfo', struct3.chaninfo, 'caption', 'C4 - I3', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I3 C3_itc C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I3] = pop_newtimef( struct3, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct3.chanlocs, 'chaninfo', struct3.chaninfo, 'caption', 'C3 - I3', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);


figure
[C4_ersp_I4 C4_itc_I4 C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I4] = pop_newtimef( struct4, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct4.chanlocs, 'chaninfo', struct4.chaninfo, 'caption', 'C4 - I4', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I4 C3_itc_I4 C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I4] = pop_newtimef( struct4, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct4.chanlocs, 'chaninfo', struct4.chaninfo, 'caption', 'C3 - I4', 'baseline',[-400 0], 'freqs', lf_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);

%% Save time frequency data on structures for low frequency components

%C3
struct1.C3_lf_tfdat = C3_tfdata_I1;
struct1.C3_lf_ersp = C3_ersp_I1;

struct2.C3_lf_tfdat = C3_tfdata_I2;
struct2.C3_lf_ersp = C3_ersp_I2;

struct3.C3_lf_tfdat = C3_tfdata_I3;
struct3.C3_lf_ersp = C3_ersp_I3;

struct4.C3_lf_tfdat = C3_tfdata_I4;
struct4.C3_lf_ersp = C3_ersp_I4;

%C4

struct1.C4_lf_tfdat = C4_tfdata_I1;
struct1.C4_lf_ersp = C4_ersp_I1;

struct2.C4_lf_tfdat = C4_tfdata_I2;
struct2.C4_lf_ersp = C4_ersp_I2;

struct3.C4_lf_tfdat = C4_tfdata_I3;
struct3.C4_lf_ersp = C4_ersp_I3;

struct4.C4_lf_tfdat = C4_tfdata_I4;
struct4.C4_lf_ersp = C4_ersp_I4;


% store frq on structures
struct1.lf_frex = frequencies;
struct2.lf_frex = frequencies;
struct3.lf_frex = frequencies;
struct4.lf_frex = frequencies;

% store times on structure
struct1.lf_times = times;
struct2.lf_times = times;
struct3.lf_times = times;
struct4.lf_times = times;


%% compute time-freq transform for gamma

gamma_range = [30 100];

figure
[C4_ersp_I1 C4_itc C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I1] = pop_newtimef( struct1, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct1.chanlocs, 'chaninfo', struct1.chaninfo, 'caption', 'C4 - I1', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I1 C3_itc C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I1] = pop_newtimef( struct1, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct1.chanlocs, 'chaninfo', struct1.chaninfo, 'caption', 'C3 - I1', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);

figure
[C4_ersp_I2 C4_itc C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I2] = pop_newtimef( struct2, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct2.chanlocs, 'chaninfo', struct2.chaninfo, 'caption', 'C4 - I2', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I2 C3_itc C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I2] = pop_newtimef( struct2, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct2.chanlocs, 'chaninfo', struct2.chaninfo, 'caption', 'C3 - I2', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);


figure
[C4_ersp_I3 C4_itc C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I3] = pop_newtimef( struct3, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct3.chanlocs, 'chaninfo', struct3.chaninfo, 'caption', 'C4 - I3', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I3 C3_itc C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I3] = pop_newtimef( struct3, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct3.chanlocs, 'chaninfo', struct3.chaninfo, 'caption', 'C3 - I3', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);


figure
[C4_ersp_I4 C4_itc_I4 C4_powbase times frequencies C4_erspboot C4_itcboot C4_tfdata_I4] = pop_newtimef( struct4, 1, 24, [-350  1000], [0] , 'topovec', 24, 'elocs', struct4.chanlocs, 'chaninfo', struct4.chaninfo, 'caption', 'C4 - I4', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);
figure
[C3_ersp_I4 C3_itc_I4 C3_powbase times frequencies C3_erspboot C3_itcboot C3_tfdata_I4] = pop_newtimef( struct4, 1, 7, [-350  1000], [0] , 'topovec', 7, 'elocs', struct4.chanlocs, 'chaninfo', struct4.chaninfo, 'caption', 'C3 - I4', 'baseline',[-400 0], 'freqs', gamma_range, 'plotphase', 'off', 'scale', 'abs','ntimesout',400,'padratio', 1);

%% Save time frequency data on structures for gamma

%C3
struct1.C3_gamma_tfdat = C3_tfdata_I1;
struct1.C3_gamma_ersp = C3_ersp_I1;

struct2.C3_gamma_tfdat = C3_tfdata_I2;
struct2.C3_gamma_ersp = C3_ersp_I2;

struct3.C3_gamma_tfdat = C3_tfdata_I3;
struct3.C3_gamma_ersp = C3_ersp_I3;

struct4.C3_gamma_tfdat = C3_tfdata_I4;
struct4.C3_gamma_ersp = C3_ersp_I4;

%C4

struct1.C4g_gamma_tfdat = C4_tfdata_I1;
struct1.C4_gamma_ersp = C4_ersp_I1;

struct2.C4_gamma_tfdat = C4_tfdata_I2;
struct2.C4_gamma_ersp = C4_ersp_I2;

struct3.C4_gamma_tfdat = C4_tfdata_I3;
struct3.C4_lf_ersp = C4_ersp_I3;

struct4.C4_gamma_tfdat = C4_tfdata_I4;
struct4.C4_gamma_ersp = C4_ersp_I4;


% store frq on structures
struct1.gamma_frex = frequencies;
struct2.gamma_frex = frequencies;
struct3.gamma_frex = frequencies;
struct4.gamma_frex = frequencies;

% store times on structure
struct1.gamma_times = times;
struct2.gamma_times = times;
struct3.gamma_times = times;
struct4.gamma_times = times;

% store number of chans
nb_chans = 32;
struct1.nbchan = nb_chans;
struct2.nbchan = nb_chans;
struct3.nbchan = nb_chans;
struct4.nbchan = nb_chans;

% store number of trials 
temp = size(struct1.data);
struct1.trials = temp(end);

temp = size(struct2.data);
struct2.trials = temp(end);

temp = size(struct3.data);
struct3.trials = temp(end);

temp = size(struct4.data);
struct4.trials = temp(end);


%% Low pass filter for LEPS

my_locut = 0.1;
my_hicut = 20;

%filter the data
struct_1_filtered = pop_eegfiltnew(struct1, 'locutoff',my_locut,'hicutoff',my_hicut,'plotfreqz',1);
struct_2_filtered = pop_eegfiltnew(struct2, 'locutoff',my_locut,'hicutoff',my_hicut,'plotfreqz',1);
struct_3_filtered = pop_eegfiltnew(struct3, 'locutoff',my_locut,'hicutoff',my_hicut,'plotfreqz',1);
struct_4_filtered = pop_eegfiltnew(struct4, 'locutoff',my_locut,'hicutoff',my_hicut,'plotfreqz',1);

struct1.data_filtered = struct_1_filtered.data;
struct2.data_filtered = struct_2_filtered.data;
struct3.data_filtered = struct_3_filtered.data;
struct4.data_filtered = struct_4_filtered.data;




