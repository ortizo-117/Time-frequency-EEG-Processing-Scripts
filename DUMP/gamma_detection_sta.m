
function EEG_out  = gamma_detection_sta(EEG,epoch_length,freqs_in,baseline_in,n_times_out,pads,win_size,NF,c_lims)
%Author: Oscar Ortiz
%Date: SEPT 29, 2022
%Description: This function is to plot gamma activity
%The input of this function is the directory where the .set files are,
%that have already been decomposed by ICA.
problem_v = []; 
 
     
     %%  extracting tf for each channel
num_chan=EEG.nbchan; % get number of channels
for i = 1:num_chan
    if i ==1
        try
            fh1 = figure;
             [ersp, ~, ~, times, frequencies, erspboot, ~, ~] = pop_newtimef( EEG, 1, i,epoch_length, [0] , 'topovec', i, 'elocs', EEG.chanlocs, 'chaninfo', ...
             EEG.chaninfo, 'caption', char(EEG.chanlocs(i).labels), ...
             'baseline',baseline_in,'freqs',freqs_in, 'plotphase', 'on', ...
             'ntimesout', n_times_out, 'padratio', pads,'winsize', win_size,'nfreqs', NF,...
             'erspmax', c_lims);
             m_sizes = size(ersp);
             out_mat = zeros(num_chan,m_sizes(1),m_sizes(2));
             out_mat(i,:,:) = ersp;
             pause(0.1)
             close(fh1)
        catch
              disp(['Problem with channel '  char(EEG.chanlocs(i).labels) '. Potentially all zeros. This is the first channel. Pressing on...'])
              problem_v = [problem_v , i];
        end
    else
        try
             fh1 = figure;
             [ersp, ~, ~, times, frequencies, erspboot, ~, ~] = pop_newtimef( EEG, 1, i,epoch_length, [0] , 'topovec', i, 'elocs', EEG.chanlocs, 'chaninfo', ...
             EEG.chaninfo, 'caption', char(EEG.chanlocs(i).labels), ...
             'baseline',baseline_in,'freqs',freqs_in, 'plotphase', 'on', ...
             'ntimesout', n_times_out, 'padratio', pads,'winsize', win_size,'nfreqs', NF,...
             'erspmax', c_lims);
             out_mat(i,:,:) = ersp;
             pause(0.1)
             close(fh1)
         catch
              disp(['Problem with channel '  char(EEG.chanlocs(i).labels) '. Potentially all zeros. Pressing on...'])
              problem_v = [problem_v , i];
        end
    end
end  
     
     %% put it into EEG structure 
     
EEG.tf_data = out_mat;
EEG.tf_frex = frequencies;
EEG.tf_time = times;
EEG.chans_with_problems = problem_v;
EEG_out = EEG;
          
end  
 


