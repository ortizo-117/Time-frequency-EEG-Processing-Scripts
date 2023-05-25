

%% load in data

S1 = load('Z:\30_Oscar_Ortiz\GluEEG\5_Matlab\2_JL\processed_data_avg.mat');
S2 = load('Z:\30_Oscar_Ortiz\GluEEG\5_Matlab\3_LL\processed_data_avg.mat');
S3 = load('Z:\30_Oscar_Ortiz\GluEEG\5_Matlab\4_JA\processed_data_avg.mat');
S1b= load('Z:\30_Oscar_Ortiz\GluEEG\5_Matlab\5_JL2\processed_data_avg.mat');
%S4 = load('Z:\30_Oscar_Ortiz\GluEEG\5_Matlab\6_JM\processed_data_avg.mat')

%% Define parameters

subs ={'S1' 'S2' 'S3' 'S1b'};

intensities = {'struct1' 'struct2' 'struct3' 'struct4'};
titles = {'I1', 'I2', 'I3', 'I4'};
grand_titles = {'I1 std', 'I1', 'I2 std', 'I2', 'I3 std', 'I3','I4 std', 'I4'};

%% Open up Gamma per subject


for i = 1:length(subs)
    figure
    for j = 1:length(intensities)
        subplot(2,2,j)
         data2plot = eval([char(subs(i)),'.',char(intensities(j)),'.C3_gamma_ersp;']);
         x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
         y = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_frex;']);
         contourf(x,y,data2plot,40,'linecolor','none');
         caxis([0 2.5])
         c = colorbar;
         c.Label.String = 'ERSP (% of baseline)';
         title(char(titles(j)))
         xlabel('Time (ms)')
         ylabel('Frequency (Hz)')
    end
    sgtitle([char(subs(i)) ' Avg Ref'])
end


%% Plot time series for gamma

for i = 1:length(subs)
    figure
    for j = 1:length(intensities)
        subplot(2,2,j)
         temp = eval([char(subs(i)),'.',char(intensities(j)),'.C3_gamma_ersp;']);
         max_data = max(temp);
         min_data = min(temp);
         x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
         plot(x,max_data)
         hold on
         plot(x,min_data)
         ylim([0 4])
         xlabel('Time (ms)')
         ylabel('Amplitude')
         legend('Max','Min')
         title(char(titles(j)))
    end
    sgtitle(char(subs(i)))
end


%% Plot Individual time series but on the same plot


for i = 1:length(subs)
    figure
    for j = 1:length(intensities)
         temp = eval([char(subs(i)),'.',char(intensities(j)),'.C3_gamma_ersp;']);
         max_data = max(temp);
         min_data = min(temp);
         x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
         plot(x,max_data)
         hold on
    end
    ylim([0 4])
    xlabel('Time (ms)')
    ylabel('Max Gamma Activity')
    legend(titles)
    sgtitle(char(subs(i)))
end

%% Plot grand averages
all_data = zeros(length(subs),length(intensities),19,400);


% Put all data into one structure
for i = 1:length(subs)
    for j = 1:length(intensities)
         data2plot = eval([char(subs(i)),'.',char(intensities(j)),'.C3_gamma_ersp;']);
         all_data(i,j,:,:) = data2plot;
    end
%     x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
%     y = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_frex;']);
%     contourf(x,y,data2plot,40,'linecolor','none');
%     caxis([0 3])
%     c = colorbar;
%     c.Label.String = 'ERSP (% of baseline)';
%     title(char(titles(j)))
%     xlabel('Time (ms)')
%     ylabel('Frequency (Hz)')
%     sgtitle(char(subs(i)))
end


% taking the subject average

grand_averages = squeeze(mean(all_data,1));

% plot spectra
figure
for j = 1:length(intensities)
    subplot(2,2,j)
    data2plot = squeeze(grand_averages(j,:,:));
    x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
    y = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_frex;']);
    contourf(x,y,data2plot,40,'linecolor','none');
    caxis([0 2])
    c = colorbar;
    c.Label.String = 'ERSP (% of baseline)';
    title(char(titles(j)))
    xlabel('Time (ms)')
end
sgtitle('Gamma spectra grand averages')


% plot max vals
my_alpha = 0.1
cMap = cool(length(intensities));


figure
for j = 1:length(intensities)
    temp = squeeze(all_data(:,j,:,:));
    data2plot = squeeze(max(temp,[],2));
    x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
    stdshade(data2plot,my_alpha,cMap(j,:),x)
    hold on
end
ylim([0 3.5])
xlim([-200 800])
xlabel('Time (ms)')
ylabel('Max Gamma Activity')
legend(grand_titles,'NumColumns',2)
sgtitle('Max Gamma Activity (n=3)')





%% Open up LF per subject


for i = 1:length(subs)
    figure
    for j = 1:length(intensities)
        subplot(2,2,j)
         data2plot = eval([char(subs(i)),'.',char(intensities(j)),'.C3_lf_ersp;']);
         x = eval([char(subs(i)),'.',char(intensities(j)),'.lf_times;']);
         y = eval([char(subs(i)),'.',char(intensities(j)),'.lf_frex;']);
         contourf(x,y,data2plot,40,'linecolor','none');
         caxis([0 6])
         c = colorbar;
         c.Label.String = 'ERSP (% of baseline)';
         title(char(titles(j)))
         xlabel('Time (ms)')
         ylabel('Frequency (Hz)')
    end
    sgtitle(char(subs(i)))
end


%% Plot time series for lf

for i = 1:length(subs)
    figure
    for j = 1:length(intensities)
        subplot(2,2,j)
         temp = eval([char(subs(i)),'.',char(intensities(j)),'.C3_lf_ersp;']);
         max_data = max(temp);
         min_data = min(temp);
         x = eval([char(subs(i)),'.',char(intensities(j)),'.lf_times;']);
         plot(x,max_data)
         hold on
         plot(x,min_data)
         ylim([0 10])
         xlabel('Time (ms)')
         ylabel('Amplitude')
         legend('Max','Min')
         title(char(titles(j)))
    end
    sgtitle(char(subs(i)))
end


%% Plot Individual time series but on the same plot


for i = 1:length(subs)
    figure
    for j = 1:length(intensities)
         temp = eval([char(subs(i)),'.',char(intensities(j)),'.C3_lf_ersp;']);
         max_data = max(temp);
         min_data = min(temp);
         x = eval([char(subs(i)),'.',char(intensities(j)),'.lf_times;']);
         plot(x,max_data)
         hold on
    end
    ylim([0 15])
    xlabel('Time (ms)')
    ylabel('Max LF Activity')
    legend(titles)
    sgtitle(char(subs(i)))
end

%% Plot grand averages
all_data_lf = zeros(length(subs),length(intensities),8,400);


% Put all data into one structure
for i = 2:length(subs)
    for j = 1:length(intensities)
         data2plot = eval([char(subs(i)),'.',char(intensities(j)),'.C3_lf_ersp;']);
         all_data_lf(i,j,:,:) = data2plot;
    end
%     x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
%     y = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_frex;']);
%     contourf(x,y,data2plot,40,'linecolor','none');
%     caxis([0 3])
%     c = colorbar;
%     c.Label.String = 'ERSP (% of baseline)';
%     title(char(titles(j)))
%     xlabel('Time (ms)')
%     ylabel('Frequency (Hz)')
%     sgtitle(char(subs(i)))
end


% taking the subject average

grand_averages = squeeze(mean(all_data_lf,1));


% plot lf spectra
figure
for j = 1:length(intensities)
    subplot(2,2,j)
    data2plot = squeeze(grand_averages(j,:,:));
    x = eval([char(subs(i)),'.',char(intensities(j)),'.lf_times;']);
    y = eval([char(subs(i)),'.',char(intensities(j)),'.lf_frex;']);
    contourf(x,y,data2plot,40,'linecolor','none');
    caxis([0 5])
    c = colorbar;
    c.Label.String = 'ERSP (% of baseline)';
    title(char(titles(j)))
    xlabel('Time (ms)')
end
sgtitle('LF spectra grand averages')


% plot max vals
my_alpha = 0.1
cMap = cool(length(intensities));


figure
for j = 1:length(intensities)
    temp = squeeze(all_data(:,j,:,:));
    data2plot = squeeze(max(temp,[],2));
    x = eval([char(subs(i)),'.',char(intensities(j)),'.gamma_times;']);
    stdshade(data2plot,my_alpha,cMap(j,:),x)
    hold on
end
ylim([0 3.5])
xlim([-200 800])
xlabel('Time (ms)')
ylabel('Max Gamma Activity')
legend(grand_titles,'NumColumns',2)
sgtitle('Max Gamma Activity (n=3)')


