%% Variables

s_rate = 1000;
time_range = [-400 1000]; % in s
time_vect = linspace(time_range(1),time_range(2),abs(((time_range(1)-time_range(2)))/1000)*s_rate);
n_trials = 20;
n_chans= 32;
chan_name = 'Cz';
pnts = 1400;





%% Read EEG file for 1st set
filter = '*.csv';
title = 'Select EEG File';
[filenameEEG1,path] = uigetfile(filter,title);
[EEGdata1, chan_names1 ] = readEEGfile(filenameEEG1);




%% Read EEG file for 2nd set
filter = '*.csv';
title = 'Select EEG File';
[filenameEEG2,path] = uigetfile(filter,title);
[EEGdata2, chan_names2 ] = readEEGfile(filenameEEG2);







%% Read EEG file for 3rd set
filter = '*.csv';
title = 'Select EEG File';
[filenameEEG3,path] = uigetfile(filter,title);
[EEGdata3, chan_names3 ] = readEEGfile(filenameEEG3);




%% Read EEG file for 4th set
filter = '*.csv';
title = 'Select EEG File';
[filenameEEG4,path] = uigetfile(filter,title);
[EEGdata4, chan_names4 ] = readEEGfile(filenameEEG4);


%% Getting Cz and reshaping data

Cz_i = find(ismember(chan_names1,chan_name));
Cz_dat = EEGdata1(:,Cz_i);

EEG1 = reshape(Cz_dat,length(time_vect),[],1);

Cz_i = find(ismember(chan_names2,chan_name));
Cz_dat = EEGdata2(:,Cz_i);

EEG2 = reshape(Cz_dat,length(time_vect),[],1);

Cz_i = find(ismember(chan_names3,chan_name));
Cz_dat = EEGdata3(:,Cz_i);

EEG3 = reshape(Cz_dat,length(time_vect),[],1);


Cz_i = find(ismember(chan_names4,chan_name));
Cz_dat = EEGdata4(:,Cz_i);

EEG4 = reshape(Cz_dat,length(time_vect),[],1);


%% concatenate EEG data and raitings

C = cat(2,EEG1,EEG2);
C = cat(2,C,EEG3);
C = cat(2,C,EEG4);



%% Reading binning information
file_name_in = 'JL_pilot.xlsx';
binned_data = get_binning_information(file_name_in);

R1_loc = find(binned_data ==1);
R2_loc = find(binned_data ==2);
R3_loc = find(binned_data ==3);
R4_loc = find(binned_data ==4);

EEG_I1 = C(:,R1_loc);
EEG_I2 = C(:,R2_loc);
EEG_I3 = C(:,R3_loc);
EEG_I4 = C(:,R4_loc);


%% Check means and plot
mean_R1 = squeeze(mean(EEG_I1,2));
mean_R2 = squeeze(mean(EEG_I2,2));
mean_R3 = squeeze(mean(EEG_I3,2));
mean_R4 = squeeze(mean(EEG_I4,2));

figure
hold on;
plot(time_vect,mean_R1,'DisplayName','0-25');
plot(time_vect,mean_R2,'DisplayName','<=25-50');
plot(time_vect,mean_R3,'DisplayName','<=50-75');
plot(time_vect,mean_R4,'DisplayName','<=75-100');
xline(0,'--r','DisplayName','Stim Onset');
lgd = legend;
lgd.NumColumns = 2;
xlabel('Time (ms)');
ylabel('Amplitude (μV)');
set(gca, 'YDir','reverse')
set(gcf,'color','w');

%% Create EEGlab structure phantom

% reshape data

EEG1_reshaped = reshape(EEGdata1',n_chans,pnts,n_trials);
EEG2_reshaped = reshape(EEGdata2',n_chans,pnts,n_trials);
EEG3_reshaped = reshape(EEGdata3',n_chans,pnts,n_trials);
EEG4_reshaped = reshape(EEGdata4',n_chans,pnts,n_trials);

% Concatenate
B = cat(3,EEG1_reshaped,EEG2_reshaped);
B = cat(3,B,EEG3_reshaped);
B = cat(3,B,EEG4_reshaped);

% Binning data
file_name_in = 'JL_pilot.xlsx';
binned_data = get_binning_information(file_name_in);

R1_loc = find(binned_data ==1);
R2_loc = find(binned_data ==2);
R3_loc = find(binned_data ==3);
R4_loc = find(binned_data ==4);

EEG_I1 = B(:,:,R1_loc);
EEG_I2 = B(:,:,R2_loc);
EEG_I3 = B(:,:,R3_loc);
EEG_I4 = B(:,:,R4_loc);


% storing data
struct1.data= EEG_I1;
struct2.data= EEG_I2;
struct3.data= EEG_I3;
struct4.data= EEG_I4;

% storing original filename name
struct1.filename = filenameEEG1;
struct2.filename = filenameEEG2;
struct3.filename = filenameEEG3;
struct4.filename = filenameEEG4;


% store sampling rate
struct1.srate = s_rate;
struct2.srate = s_rate;
struct3.srate = s_rate;
struct4.srate = s_rate;


% store times
struct1.times = time_vect;
struct2.times = time_vect;
struct3.times = time_vect;
struct4.times = time_vect;


% store locations
struct1.chanlocs = chanlocs;
struct2.chanlocs = chanlocs;
struct3.chanlocs = chanlocs;
struct4.chanlocs = chanlocs;

% store chaninfo
struct1.chaninfo = chaninfo;
struct2.chaninfo = chaninfo;
struct3.chaninfo = chaninfo;
struct4.chaninfo = chaninfo;

%store xmin
struct1.xmin = min(struct1.times)/1000;
struct2.xmin = min(struct2.times)/1000;
struct3.xmin = min(struct3.times)/1000;
struct4.xmin = min(struct4.times)/1000;

%store xmax

struct1.xmax = max(struct1.times)/1000;
struct2.xmax = max(struct2.times)/1000;
struct3.xmax = max(struct3.times)/1000;
struct4.xmax = max(struct4.times)/1000;


%store pnts
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






