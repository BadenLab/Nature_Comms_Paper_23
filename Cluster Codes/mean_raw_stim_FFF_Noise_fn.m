function mean_raw_stim_arr = mean_raw_stim_FFF_Noise_fn(stimulus_arr,trig_times_vec,p)
% Calculate the mean of the raw stimuli
% trig_times_vec: the vector input is 'trig_times_vec_trunc'

mean_raw_stim_arr_loop = zeros(p.Num_STE_bins,1);

modNum_FinalFrames = mod(p.Num_trigs,p.stim_frames);

if p.Time_Choice == 1 || p.Mean_Stim_Choice  == 1 % stimulus frames
    
    if modNum_FinalFrames==0 % Every noise chunk is repeated a whole number of times
        
        for i = 1:p.Num_Raw_Stim
            mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + stimulus_arr(i:p.Num_STE_bins+i-1);
        end
        
    else % modNum_FinalFrames~==0 % The final noise chunk does not contain a full set of frames
        
        for i = 1:p.Num_Raw_Stim
            if modNum_FinalFrames >= p.Num_STE_bins+i-1 % i
                mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + p.Num_FNoise_rep_ceil*stimulus_arr(i:p.Num_STE_bins+i-1);
            else % modNum_FinalFrames < p.Num_STE_bins+i-1 % i
                mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + (p.Num_FNoise_rep_ceil - 1)*stimulus_arr(i:p.Num_STE_bins+i-1);
            end
        end
        
    end
    
elseif p.Time_Choice == 2 && p.Mean_Stim_Choice  == 2 % define own time grid
    
    loop_stim_count = 0; % To count the number of stimuli when choose own time and final noise chunk does not contain a full set of frames
    
    for i = 1:p.Num_Raw_Stim
        
        loop_sample_times_vec = p.Num_Raw_Stim_t_vec(i) + p.stim_timesample_vec;
        loop_index_vec        = zeros(p.Num_STE_bins,1);
        
        % Find which stimulis frame was active at each of the sample times
        % Indices are recorded (use sum to find index as
        % number of triggers fired by the j'th time).
        for j = 1:p.Num_STE_bins
            loop_index_vec(j) = sum((trig_times_vec-loop_sample_times_vec(j))<=1e-10);
            % was <0, make <= 0 to allow for equality case which arrises here.
            % make 1e-10 since can get values like 4.4409*1e-16 for
            % difference beteween first entries of 'trig_times_vec' and
            % 'loop_sample_times_vec' and hence loop_index_vec(1) = 0.
        end
        
        if modNum_FinalFrames==0 % Every noise chunk is repeated a whole number of times
            
            mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + stimulus_arr(loop_index_vec);
            
        else % modNum_FinalFrames~==0 % The final noise chunk does not contain a full set of frames
            
            if modNum_FinalFrames >= loop_index_vec(end)
                mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + p.Num_FNoise_rep_ceil*stimulus_arr(loop_index_vec);
                loop_stim_count = loop_stim_count + p.Num_FNoise_rep_ceil;
            else % modNum_FinalFrames < loop_index_vec(end)
                mean_raw_stim_arr_loop = mean_raw_stim_arr_loop + (p.Num_FNoise_rep_ceil - 1)*stimulus_arr(loop_index_vec);
                loop_stim_count = loop_stim_count + (p.Num_FNoise_rep_ceil-1);
            end
            
        end
        
    end
    
end

if modNum_FinalFrames==0 % Every noise chunk is repeated a whole number of times
    
    mean_raw_stim_arr = mean_raw_stim_arr_loop/p.Num_Raw_Stim;
    
else % modNum_FinalFrames~==0 % The final noise chunk does not contain a full set of frames
    
    if p.Time_Choice == 1 || p.Mean_Stim_Choice  == 1 % stimulus frames
        num_final_stim = modNum_FinalFrames - p.Num_STE_bins + 1; % Number of raw stimuli in final noise chunk
        
        if num_final_stim > 0
            Mod_Num_Raw_Stim = p.Num_Raw_Stim*(p.Num_FNoise_rep_ceil - 1) + num_final_stim;
        else
            Mod_Num_Raw_Stim = p.Num_Raw_Stim*(p.Num_FNoise_rep_ceil - 1);
        end
        
    elseif p.Time_Choice == 2 && p.Mean_Stim_Choice  == 2 % define own time grid
        
        Mod_Num_Raw_Stim = loop_stim_count; % Was going to do the following for num_final_stim but incorrect, need to subtract fdifferent number but hard to determine: sum(p.Num_Raw_Stim_t_vec<trig_times_vec(modNum_FinalFrames+1))- p.Num_STE_bins + 1;% number of raw stim that can occur in  truncated final chunk
        
    end
    
    mean_raw_stim_arr = mean_raw_stim_arr_loop/Mod_Num_Raw_Stim;
    
end

