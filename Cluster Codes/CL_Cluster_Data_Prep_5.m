%% Cluster Data Preparation 5
% This code is used to gather and preprocess the data to be clustered
% It then calls the clustering function
% 23,09,2021 Onwards - Adding Chirp 3 stimulus

clear all;
clc;

%% Make Choices

%%% Choose whether to use parallel processing and number of cores
p.Parpool   = 1; % 1: Yes, 2: No.
p.Num_Cores = 6; % 6

%%% Choose data sets
% [set1,set2,...,setN]
% [15_01_2021 Phase_01 BdWth05 Zfish,     03_12_2020 Phase_01 BdWth05 ZFish,...
%  02_12_2020 Phase_00 BdWth05 ZFish,     27_11_2020 Phase_00 BdWth06thr3 Zfish,...
%  18_11_2020 Phase_01 BdWth06thr3 Zfish, 16_11_2020 Phase_00 BdWth06thr3 Zfish,...
%  24_08_2020 Phase_00 BdWth05 Zfish]
data_set_vec     = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % [1,1,0,0,1,1,1] [0,0,0,0,1,0,0] 1 [1 1 1 1] [1 0 0 0] [1 1 1 1 1 1] [1 1 1 1 1] [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]   [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1] [1 1 1 1 1 1 1 1 1 1 1 1]
data_set_indices = find(data_set_vec==1);
p.Num_data_sets  = length(data_set_indices);% length(data_set_vec)

%%% Choose which objects to cluster on
% 1 = yes, 2 = no.
% 1. FFF;                                 (all data sets)
% 2. Chirp;                               (all data sets)
% 3. FFF Noise;                           (most recent data sets)
% 4. Gratings;                            (all data sets)
% 5. CNoise: Full RF Size;                (all data sets)
% 6. CNoise: Full RF Ellipticity;         (all data sets)
% 7. CNoise: Full RF Dominant Axis Angle; (all data sets)
% 8. FFF2                                 (6 colours and triggers);
% 9. Chirp2                               (3 triggers);
% 10. Silent Substitution                 (18 triggers);
% 11. Contrast Steps                      (10 triggers);
% 12. Chirp3                              (4 triggers).
p.Obj_clust_vec = [0,0,1,0,0,0,0,0,0,0,0,0];%[1,1,1,0,1,1,1]; [0,0,1,0,0,0,0,1,0,1,1];
%[0,0,1,0,0,0,0,1,0,0,1,1]
p.Num_clust_obj = sum(p.Obj_clust_vec);

%%% Choose which objects to plot
% 1 = yes, 2 = no.
% 1. Cell positions;                      (scatter plot)
% 2. Mean RF contour;                     (line plot)
% 3. FFF;                                 (trace & heatmap)
% 4. Chirp;                               (trace & heatmap)
% 5. FFF Noise;                           (trace & heatmap)
% 6. Gratings;                            (trace & heatmap)
% 7. CNoise: Full RF Size;                (histogram)
% 8. CNoise: Full RF Ellipticity;         (histogram)
% 9. CNoise: Full RF Dominant Axis Angle; (histogram)
% 10. FFF2 (6 colours and triggers);
% 11. Chirp2 (3 triggers);
% 12. Silent Substitution (18 triggers);
% 13. Contrast Steps (10 triggers);
% 14. Chirp3 (4 triggers);
% 15. Data Set.                           (histogram) --> was 9 then 10. 23,09,2021: now 15
p.Obj_plot_vec     = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,1];%[0,0,1,1,0,1,0,0,0]  [0,0,0,0,1,0,0,0,0,1,0,1,1,0];
%[0,0,0,0,1,0,0,0,0,1,0,0,1,1,1]
p.Plot_per_obj_vec = [1,1,2,2,2,2,1,1,1,2,2,2,2,2,1];
p.Num_plot_obj     = sum(p.Obj_plot_vec);
p.Num_plots        = dot(p.Obj_plot_vec,p.Plot_per_obj_vec);

%%% QI Thresholds
FFF_QI_thresh            = 0.20; % 0.15 (could do 0.12-0.14)
Chirp_QI_thresh          = 0.6;  % 0.4  (could do to 0.35)
FFF_Noise_QI_thresh      = 2.5;  % 0.3  (could do 0.2)
Gratings_400px_QI_thresh = 0.2;  % 0.2  (could do 0.15 (though lets all/most data 1 in)
CNoise_QI_thresh         = 0;    % 1    (minumun number of RF pixels, could do 2/3 - strict ineq so 1 gives RF > 1)
FFF2_QI_thresh           = 0.4; % was 0.2 (09,09,2021)
Chirp2_QI_thresh         = 0.4; % 
SSub_QI_thresh           = 0.2; % 
CSteps_QI_thresh         = 0.4; % PAR Mod 22,09,2021 was '0.4' PAR Mod 23,11,2021 was '0.085' put back to 0.4
Chirp3_QI_thresh         = 0.4; % 

%%% Choose whether to use CNoise QI
% Don't use in general as want even 1 pixel RFs so all spiking cells would get through
% 1. Yes;
% 2. No,
Use_CNoise_QI = 2;


% Data set names
All_DataFiles = {};
if p.Obj_clust_vec(1) == 1 || p.Obj_plot_vec(3) == 1 % FFF
    %%% Zebrafish data
    %FFF_DataFiles = {'FFF'};                                      % Data 1
    %FFF_DataFiles = {'03_12_2020_ZF_Ph01_FFF_All_Cells_Corr'};   % Data 2 ('03_12_2020_ZF_Ph01_FFF_100_Cells_Corr')
    %FFF_DataFiles = {'02_12_2020_ZF_Ph00_FFF_1_All_Cells_Marv'}; % Data 3 ('02_12_2020_ZF_Ph00_FFF_2_All_Cells_Marv')
    %FFF_DataFiles = {'27_11_2020_ZF_Ph00_FFF_1_All_Cells_Marv'}; % Data 4 ('27_11_2020_ZF_Ph00_FFF_2_All_Cells_Marv')
    %FFF_DataFiles = {'18_11_2020_ZF_Ph01_FFF_1_All_Cells'};      % Data 5 ('18_11_2020_ZF_Ph01_FFF_2_All_Cells')
    %FFF_DataFiles = {'16_11_2020_ZF_Ph00_FFF_All_Cells_Marv'};   % Data 6
    %FFF_DataFiles = {'24_08_2020_ZF_Ph00_FFF_All_Cells_Marv'};   % Data 7
    %FFF_DataFiles = {'FFF',...
%                      '03_12_2020_ZF_Ph01_FFF_All_Cells_Corr',...
%                      '02_12_2020_ZF_Ph00_FFF_1_All_Cells_Marv',...
%                      '27_11_2020_ZF_Ph00_FFF_1_All_Cells_Marv'};   % Data 1--4
%     FFF_DataFiles = {'FFF',...
%                      '03_12_2020_ZF_Ph01_FFF_All_Cells_Corr',...
%                      '02_12_2020_ZF_Ph00_FFF_1_All_Cells_Marv',...
%                      '27_11_2020_ZF_Ph00_FFF_1_All_Cells_Marv',...
%                      '18_11_2020_ZF_Ph01_FFF_1_All_Cells',...
%                      '16_11_2020_ZF_Ph00_FFF_All_Cells_Marv',...
%                      '24_08_2020_ZF_Ph00_FFF_All_Cells_Marv'};   % Data 1--7
    %FFF_DataFiles = {'FFF','03_12_2020_ZF_Ph01_FFF_100_Cells_Corr','02_12_2020_ZF_Ph00_FFF_1_All_Cells_Marv','27_11_2020_ZF_Ph00_FFF_1_All_Cells_Marv'};
    %FFF_DataFiles = {'FFF','03_12_2020_ZF_Ph01_FFF_100_Cells_Corr','#','#','18_11_2020_ZF_Ph01_FFF_1_All_Cells','16_11_2020_ZF_Ph00_FFF_All_Cells_Marv','24_08_2020_ZF_Ph00_FFF_All_Cells_Marv'};
    %%% Chicken data
%     FFF_DataFiles = {'06_12_2019_Chick_UV_FFF_Marv','11_12_2019_Chick_UV_FFF_Marv',...
%                      '18_12_2019_Chick_UV_FFF_Marv','20_12_2019_Chick_UV_FFF_Marv'}; % UV Data
%     FFF_DataFiles = {'09_09_2020_Chick_NoUV_FFF_Marv','10_09_2020_Chick_NoUV_FFF_Marv',...
%                      '16_10_2020_Chick_NoUV_FFF_Marv','28_01_2021_Chick_NoUV_FFF_Marv',...
%                      '29_01_2021_Chick_NoUV_FFF_1st_Marv','29_01_2021_Chick_NoUV_FFF_2nd_Marv'}; % No UV Data
    All_DataFiles = [All_DataFiles;FFF_DataFiles];
end
if p.Obj_clust_vec(2) == 1 || p.Obj_plot_vec(4) == 1 % Chirp
    %Chirp_DataFiles = {'Chirp'};                                      % Data 1
    %Chirp_DataFiles = {'03_12_2020_ZF_Ph01_Chirp_1_All_Cells'};      % Data 2 ('03_12_2020_ZF_Ph01_Chirp_1_100_Cells_Corr','03_12_2020_ZF_Ph01_Chirp_2_All_Cells','03_12_2020_ZF_Ph01_Chirp_2_100_Cells_Corr')
    %Chirp_DataFiles = {'02_12_2020_ZF_Ph00_Chirp_1_All_Cells_Marv'}; % Data 3 ('02_12_2020_ZF_Ph00_Chirp_2_All_Cells_Marv','02_12_2020_ZF_Ph00_Chirp_3_All_Cells_Marv')
    %Chirp_DataFiles = {'27_11_2020_ZF_Ph00_Chirp_1_All_Cells_Marv'}; % Data 4 ('27_11_2020_ZF_Ph00_Chirp_2_All_Cells_Marv')
    %Chirp_DataFiles = {'18_11_2020_ZF_Ph01_Chirp_1_All_Cells'};      % Data 5 ('18_11_2020_ZF_Ph01_Chirp_2_All_Cells')
    %Chirp_DataFiles = {'16_11_2020_ZF_Ph00_Chirp_All_Cells_Marv'};   % Data 6
    %                                                                 % Data 7: no Chirp
    %Chirp_DataFiles = {'Chirp',...
%                      '03_12_2020_ZF_Ph01_Chirp_1_All_Cells',...
%                      '02_12_2020_ZF_Ph00_Chirp_1_All_Cells_Marv',...
%                      '27_11_2020_ZF_Ph00_Chirp_1_All_Cells_Marv'};   % Data 1--4
    Chirp_DataFiles = {'Chirp',...
                     '03_12_2020_ZF_Ph01_Chirp_1_All_Cells',...
                     '02_12_2020_ZF_Ph00_Chirp_1_All_Cells_Marv',...
                     '27_11_2020_ZF_Ph00_Chirp_1_All_Cells_Marv',...
                     '18_11_2020_ZF_Ph01_Chirp_1_All_Cells',...
                     '16_11_2020_ZF_Ph00_Chirp_All_Cells_Marv'};   % Data 1--6
    %Chirp_DataFiles = {'Chirp','03_12_2020_ZF_Ph01_Chirp_1_100_Cells_Corr','02_12_2020_ZF_Ph00_Chirp_1_All_Cells_Marv','27_11_2020_ZF_Ph00_Chirp_1_All_Cells_Marv'};
    %Chirp_DataFiles = {'Chirp','03_12_2020_ZF_Ph01_Chirp_1_100_Cells_Corr','#','#','18_11_2020_ZF_Ph01_Chirp_1_All_Cells','16_11_2020_ZF_Ph00_Chirp_All_Cells_Marv',[]};
    All_DataFiles   = [All_DataFiles;Chirp_DataFiles];
end
if p.Obj_clust_vec(3) == 1 || p.Obj_plot_vec(5) == 1 % FFF Noise
    %FFF_Noise_DataFiles = {'FFF_Noise'};                                    % Data 1
    %FFF_Noise_DataFiles = {'03_12_2020_ZF_Ph01_FFF_Noise_All_Cells'};      % Data 2 ('03_12_2020_ZF_Ph01_FFF_Noise_100_Cells_Corr')
    %FFF_Noise_DataFiles = {'02_12_2020_ZF_Ph00_FFFNoise_All_Cells_Marv'};  % Data 3
    %FFF_Noise_DataFiles = {'27_11_2020_ZF_Ph00_FFFNoise_All_Cells_Marv'};  % Data 4
    %                                                                       % Data 5: no FFF Noise
    %                                                                       % Data 6: no FFF Noise
    %                                                                       % Data 7: no FFF Noise
%     FFF_Noise_DataFiles = {'FFF_Noise',...
%         '03_12_2020_ZF_Ph01_FFF_Noise_All_Cells',...
%         '02_12_2020_ZF_Ph00_FFFNoise_All_Cells_Marv',...
%         '27_11_2020_ZF_Ph00_FFFNoise_All_Cells_Marv'};   % Data 1--4
    %FFF_Noise_DataFiles = {'FFF_Noise','03_12_2020_ZF_Ph01_FFF_Noise_100_Cells_Corr','02_12_2020_ZF_Ph00_FFFNoise_All_Cells_Marv','27_11_2020_ZF_Ph00_FFFNoise_All_Cells_Marv'};
    %FFF_Noise_DataFiles = {'FFF_Noise','03_12_2020_ZF_Ph01_FFF_Noise_100_Cells_Corr','#','#',[],[],[]};
    
    %%% Chicken data
%     FFF_Noise_DataFiles = {'10_06_2021_Chick_FFFNoise_Marv_2'};     % 26,05,2021 Onwards Data 2
%     FFF_Noise_DataFiles = {'11_06_2021_2nd_Chick_FFFNoise_Marv_3'}; % 26,05,2021 Onwards Data 3
    %FFF_Noise_DataFiles = {'10_06_2021_Chick_FFFNoise_Marv_2','11_06_2021_2nd_Chick_FFFNoise_Marv_3'}; % 26,05,2021 Onwards Data 2 & 3
    
    %FFF_Noise_DataFiles = {'03_08_2021_Ph01_Chick_FFFNoise_Marv'}; % 03_08_21 (Phase_01) - August Data 1
    %FFF_Noise_DataFiles = {'04_08_2021_Ph01_Chick_FFFNoise_Marv'}; % 04_08_21 (Phase_01) - August Data 2
    %FFF_Noise_DataFiles = {'05_08_2021_Ph01_Chick_FFFNoise_Marv'}; % 05_08_21 (Phase_01) - August Data 3
    %FFF_Noise_DataFiles = {'06_08_2021_Ph00_Chick_FFFNoise_Marv'}; % 06_08_21 (Phase_00) - August Data 4
    %FFF_Noise_DataFiles = {'12_08_2021_Ph02_Chick_FFFNoise_Marv'}; % 12_08_21 (Phase_02) - August Data 5
    %FFF_Noise_DataFiles = {'13_08_2021_Ph00_Chick_FFFNoise_Marv'}; % 13_08_21 (Phase_00) - August Data 6
    
%     FFF_Noise_DataFiles = {'03_08_2021_Ph01_Chick_FFFNoise_Marv',...  % August Data 1-6
%                            '04_08_2021_Ph01_Chick_FFFNoise_Marv',...
%                            '05_08_2021_Ph01_Chick_FFFNoise_Marv',...
%                            '06_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '12_08_2021_Ph02_Chick_FFFNoise_Marv',...
%                            '13_08_2021_Ph00_Chick_FFFNoise_Marv'};
                       
%     FFF_Noise_DataFiles = {'04_08_2021_Ph02_Chick_FFFNoise_Marv',...  % New August Data for 18th Nov Clustering (Clustering Done 22nd Nov)
%                            '05_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '06_08_2021_2nd_Ph00_Chick_FFFNoise_Marv',...
%                            '11_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '12_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '13_08_2021_Ph01_Chick_FFFNoise_Marv',...
%                            '14_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '17_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '19_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '19_08_2021_Ph01_Chick_FFFNoise_Marv',...
%                            '20_08_2021_Ph00_Chick_FFFNoise_Marv',...
%                            '21_08_2021_Ph00_Chick_FFFNoise_Marv'};

FFF_Noise_DataFiles = {'04_08_2021_Ph01_Chick_FFFNoise_Marv',...  % 17 August Datasets (FFFNoise, FFF2, CStep, Chirp3 clustering) - 22,11,2021
                      '04_08_2021_Ph02_Chick_FFFNoise_Marv',...
                      '05_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '05_08_2021_Ph01_Chick_FFFNoise_Marv',...
                      '06_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '06_08_2021_2nd_Ph00_Chick_FFFNoise_Marv',...
                      '11_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '12_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '12_08_2021_Ph02_Chick_FFFNoise_Marv',...
                      '13_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '13_08_2021_Ph01_Chick_FFFNoise_Marv',...
                      '14_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '17_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '19_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '19_08_2021_Ph01_Chick_FFFNoise_Marv',...
                      '20_08_2021_Ph00_Chick_FFFNoise_Marv',...
                      '21_08_2021_Ph00_Chick_FFFNoise_Marv'};
                       
    All_DataFiles       = [All_DataFiles;FFF_Noise_DataFiles];
end
if p.Obj_clust_vec(4) == 1 || p.Obj_plot_vec(6) == 1 % Gratings
    Gratings_400px_DataFiles = {'Gratings_400px'};                                    % Data 1
    %Gratings_400px_DataFiles = {'03_12_2020_ZF_Ph01_Gratings_All_Cells'};            % Data 2 ('03_12_2020_ZF_Ph01_Gratings_100_Cells_Corr')
    %Gratings_400px_DataFiles = {'02_12_2020_ZF_Ph00_Gratings_400px_All_Cells_Marv'}; % Data 3
    %Gratings_400px_DataFiles = {'27_11_2020_ZF_Ph00_Gratings_400px_All_Cells_Marv'}; % Data 4 ('27_11_2020_ZF_Ph00_Gratings_150px_All_Cells_Marv')
    %Gratings_400px_DataFiles = {'18_11_2020_ZF_Ph01_Gratings_150px_All_Cells'};      % Data 5  (NB: this is 150px not 400px)
    %Gratings_400px_DataFiles = {'16_11_2020_ZF_Ph00_Gratings_200px_All_Cells_Marv'}; % Data 6: (NB: this is 200px not 400px)
    %                                                                                 % Data 7: no Gratings (150px or otherwise)
    %Gratings_400px_DataFiles = {'Gratings_400px','03_12_2020_ZF_Ph01_Gratings_100_Cells_Corr','02_12_2020_ZF_Ph00_Gratings_400px_All_Cells_Marv','27_11_2020_ZF_Ph00_Gratings_400px_All_Cells_Marv'};
    %Gratings_400px_DataFiles = {'Gratings_400px','03_12_2020_ZF_Ph01_Gratings_100_Cells_Corr','#','#','18_11_2020_ZF_Ph01_Gratings_400px_All_Cells',[],[]};
    All_DataFiles            = [All_DataFiles;Gratings_400px_DataFiles];
end
if sum(p.Obj_clust_vec(5:7)) > 0 || sum(p.Obj_plot_vec([1,2,7,8,9])) > 0 % CNoise
    %CNoise_20px_DataFiles = {'data_CNoise_20px_Paul_zf_RF_test_20_AllStat_2'};              % Data 1 ('data_CNoise_20px_Paul_zf_RF_test_20_AllStat' old)
    %CNoise_20px_DataFiles = {'data_03_12_2020_ZF_Ph01_CNoise_20px_AC_RF_1_AllStat_2'};      % Data 2 % (data_CNoise_20px_Paul_zf_RF_test_23_AllStat) ('data_03_12_2020_ZF_Ph01_CNoise_20px_AC_RF_1_AllSta' old)
    %CNoise_20px_DataFiles = {'data_02_12_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat_2'}; % Data 3 ('data_02_12_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat' old)
    %CNoise_20px_DataFiles = {'data_27_11_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat_2'}; % Data 4 ('data_27_11_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat' old)
    %CNoise_20px_DataFiles = {'data_18_11_2020_ZF_Ph01_CNoise_20px_AC_Marv_RF_1_AllStat'}; % Data 5
    %CNoise_20px_DataFiles = {'data_16_11_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat_2'}; % Data 6 ('data_16_11_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat' old)
    %CNoise_20px_DataFiles = {'data_24_08_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat'}; % Data 7
    %CNoise_20px_DataFiles = {'data_CNoise_20px_Paul_zf_RF_test_20_AllStat_2',...
%         'data_03_12_2020_ZF_Ph01_CNoise_20px_AC_RF_1_AllStat_2',...
%         'data_02_12_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat_2',...
%         'data_27_11_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat_2'};   % Data 1--4
    %CNoise_20px_DataFiles = {'data_CNoise_20px_Paul_zf_RF_test_20_AllStat','data_03_12_2020_ZF_Ph01_CNoise_20px_AC_RF_1_AllStat',...
    %                         'data_02_12_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat','data_27_11_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat',...
    %                         'data_18_11_2020_ZF_Ph01_CNoise_20px_AC_Marv_RF_1_AllStat','data_16_11_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat',...
    %                         'data_24_08_2020_ZF_Ph00_CNoise_20px_AC_Marv_RF_1_AllStat'};
    
    %%% Chicken data
%     CNoise_20px_DataFiles = {'22_04_2021_Chick_CNoise_Marv_2_RF_1_AllStat'};
%     CNoise_20px_DataFiles = {'10_06_2021_Chick_CNoise_Marv_2_RF_1_AllStat'};     % 26,05,2021 Onwards Data 2
%     CNoise_20px_DataFiles = {'11_06_2021_2nd_Chick_CNoise_Marv_3_RF_1_AllStat'}; % 26,05,2021 Onwards Data 3
    CNoise_20px_DataFiles = {'10_06_2021_Chick_CNoise_Marv_2_RF_1_AllStat','11_06_2021_2nd_Chick_CNoise_Marv_3_RF_1_AllStat'}; % 26,05,2021 Onwards Data 2 & 3
    
    All_DataFiles         = [All_DataFiles;CNoise_20px_DataFiles];
end
if p.Obj_clust_vec(8) == 1 || p.Obj_plot_vec(10) == 1 % FFF2
    %%% Chicken data
    %FFF2_DataFiles = {'22_04_2021_Chick_FFF2_Marv_2'};      % Data 1 % Old: 22_04_2021_Chick_FFF2_Marv
    %     FFF2_DataFiles = {'10_06_2021_Chick_FFF2_Marv_2'};     % Data 2
    %     FFF2_DataFiles = {'11_06_2021_2nd_Chick_FFF2_Marv_3'}; % Data 3
    %FFF2_DataFiles = {'10_06_2021_Chick_FFF2_Marv_2','11_06_2021_2nd_Chick_FFF2_Marv_3'}; % Data 2 & 3
    
    %FFF2_DataFiles = {'03_08_2021_Ph01_Chick_FFF2_Marv'}; % 03_08_21 (Phase_01) - August Data 1
    %FFF2_DataFiles = {'04_08_2021_Ph01_Chick_FFF2_Marv'}; % 04_08_21 (Phase_01) - August Data 2
    %FFF2_DataFiles = {'05_08_2021_Ph01_Chick_FFF2_Marv'}; % 05_08_21 (Phase_01) - August Data 3
    %FFF2_DataFiles = {'06_08_2021_Ph00_Chick_FFF2_Marv'}; % 06_08_21 (Phase_00) - August Data 4
    %FFF2_DataFiles = {'12_08_2021_Ph02_Chick_FFF2_Marv'}; % 12_08_21 (Phase_02) - August Data 5
    %FFF2_DataFiles = {'13_08_2021_Ph00_Chick_FFF2_Marv'}; % 13_08_21 (Phase_00) - August Data 6
    
    %     FFF2_DataFiles = {'03_08_2021_Ph01_Chick_FFF2_Marv',...  % August Data 1-6
    %         '04_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '05_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '06_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '12_08_2021_Ph02_Chick_FFF2_Marv',...
    %         '13_08_2021_Ph00_Chick_FFF2_Marv'};
    
    %         FFF2_DataFiles = {'03_08_2021_Ph01_Chick_FFF2_Marv',...  % 16 August Datasets - 04,11,2021
    %         '04_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '05_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '05_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '06_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '06_08_2021_2nd_Ph00_Chick_FFFNoise_Marv',...
    %         '11_08_2021_Ph00_Chick_FFFNoise_Marv',...
    %         '12_08_2021_Ph02_Chick_FFF2_Marv',...
    %         '13_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '13_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '14_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '17_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '18_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '19_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '20_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '21_08_2021_Ph00_Chick_FFF2_Marv'};
    
    %             FFF2_DataFiles = {'03_08_2021_Ph01_Chick_FFF2_Marv',...  % 20 August Datasets - 15,11,2021
    %         '04_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '04_08_2021_Ph02_Chick_FFF2_Marv',...% New
    %         '05_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '05_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '06_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '06_08_2021_2nd_Ph00_Chick_FFFNoise_Marv',... % Missaved name as FFFNoise, now correced
    %         '11_08_2021_Ph00_Chick_FFFNoise_Marv',...     % Missaved name as FFFNoise, now correced
    %         '12_08_2021_Ph00_Chick_FFF2_Marv',...% New
    %         '12_08_2021_Ph02_Chick_FFF2_Marv',...
    %         '13_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '13_08_2021_Ph01_Chick_FFF2_Marv',...
    %         '14_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '17_08_2021_Ph00_Chick_FFF2_Marv_2',...% Modified
    %         '18_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '18_08_2021_Ph01_Chick_FFF2_Marv',...% New
    %         '19_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '19_08_2021_Ph01_Chick_FFF2_Marv',...% New
    %         '20_08_2021_Ph00_Chick_FFF2_Marv',...
    %         '21_08_2021_Ph00_Chick_FFF2_Marv'};
    
    FFF2_DataFiles = {'04_08_2021_Ph01_Chick_FFF2_Marv',...  % 17 August Datasets (FFFNoise, FFF2, CStep, Chirp3 clustering) - 22,11,2021
                      '04_08_2021_Ph02_Chick_FFF2_Marv',...
                      '05_08_2021_Ph00_Chick_FFF2_Marv',...
                      '05_08_2021_Ph01_Chick_FFF2_Marv',...
                      '06_08_2021_Ph00_Chick_FFF2_Marv',...
                      '06_08_2021_2nd_Ph00_Chick_FFF2_Marv',...
                      '11_08_2021_Ph00_Chick_FFF2_Marv',...
                      '12_08_2021_Ph00_Chick_FFF2_Marv',...
                      '12_08_2021_Ph02_Chick_FFF2_Marv',...
                      '13_08_2021_Ph00_Chick_FFF2_Marv',...
                      '13_08_2021_Ph01_Chick_FFF2_Marv',...
                      '14_08_2021_Ph00_Chick_FFF2_Marv',...
                      '17_08_2021_Ph00_Chick_FFF2_Marv_2',...
                      '19_08_2021_Ph00_Chick_FFF2_Marv',...
                      '19_08_2021_Ph01_Chick_FFF2_Marv',...
                      '20_08_2021_Ph00_Chick_FFF2_Marv',...
                      '21_08_2021_Ph00_Chick_FFF2_Marv'};
    
    All_DataFiles = [All_DataFiles;FFF2_DataFiles];
end
if p.Obj_clust_vec(9) == 1 || p.Obj_plot_vec(11) == 1 % Chirp2
    %%% Chicken data
%     Chirp2_DataFiles = {'22_04_2021_Chick_Chirp2_Marv_2'};      % Data 1
%     Chirp2_DataFiles = {'10_06_2021_Chick_Chirp2_Marv_2'};      % Data 2
%     Chirp2_DataFiles = {'11_06_2021_2nd_Chick_Chirp2_Marv_3'};  % Data 3
    Chirp2_DataFiles = {'10_06_2021_Chick_Chirp2_Marv_2','11_06_2021_2nd_Chick_Chirp2_Marv_3'}; % Data 2 % 3
    
    All_DataFiles = [All_DataFiles;Chirp2_DataFiles];
end
if p.Obj_clust_vec(10) == 1 || p.Obj_plot_vec(12) == 1 % SSub
    %%% Chicken data
    %SSub_DataFiles = {'22_04_2021_Chick_SSub_Marv_2'};            % Data 1
    
    %SSub_DataFiles = {'03_08_2021_Ph01_Chick_SSub_Marv'}; % 03_08_21 (Phase_01) - August Data 1
    %SSub_DataFiles = {'04_08_2021_Ph01_Chick_SSub_Marv'}; % 04_08_21 (Phase_01) - August Data 2
    %SSub_DataFiles = {'05_08_2021_Ph01_Chick_SSub_Marv'}; % 05_08_21 (Phase_01) - August Data 3
    %SSub_DataFiles = {'06_08_2021_Ph00_Chick_SSub_Marv'}; % 06_08_21 (Phase_00) - August Data 4
    %SSub_DataFiles = {'12_08_2021_Ph02_Chick_SSub_Marv'}; % 12_08_21 (Phase_02) - August Data 5
    %SSub_DataFiles = {'13_08_2021_Ph00_Chick_SSub_Marv'}; % 13_08_21 (Phase_00) - August Data 6
    
    SSub_DataFiles = {'03_08_2021_Ph01_Chick_SSub_Marv',...   % August Data 1-6
                      '04_08_2021_Ph01_Chick_SSub_Marv',...
                      '05_08_2021_Ph01_Chick_SSub_Marv',...
                      '06_08_2021_Ph00_Chick_SSub_Marv',...
                      '12_08_2021_Ph02_Chick_SSub_Marv',...
                      '13_08_2021_Ph00_Chick_SSub_Marv'};
    
    All_DataFiles = [All_DataFiles;SSub_DataFiles];
end
if p.Obj_clust_vec(11) == 1 || p.Obj_plot_vec(13) == 1 % CSteps
    %%% Chicken data
    %CSteps_DataFiles = {'22_04_2021_Chick_CSteps_Marv_2'};                                    % Data 1
    
    
    % No data - stim not used % 03_08_21 (Phase_01) - August Data 1
    %CSteps_DataFiles = {'04_08_2021_Ph01_Chick_CStep_Marv'}; % 04_08_21 (Phase_01) - August Data 2
    %CSteps_DataFiles = {'05_08_2021_Ph01_Chick_CStep_Marv'}; % 05_08_21 (Phase_01) - August Data 3
    %CSteps_DataFiles = {'06_08_2021_Ph00_Chick_CStep_Marv'}; % 06_08_21 (Phase_00) - August Data 4
    %CSteps_DataFiles = {'12_08_2021_Ph02_Chick_CStep_Marv'}; % 12_08_21 (Phase_02) - August Data 5
    %CSteps_DataFiles = {'13_08_2021_Ph00_Chick_CStep_Marv'}; % 13_08_21 (Phase_00) - August Data 6
    
%     CSteps_DataFiles = {'04_08_2021_Ph01_Chick_SSub_Marv',...   % August Data 1-6 -- Seem to have been using wrong data here if used it???
%                         '05_08_2021_Ph01_Chick_SSub_Marv',...
%                         '06_08_2021_Ph00_Chick_SSub_Marv',...
%                         '12_08_2021_Ph02_Chick_SSub_Marv',...
%                         '13_08_2021_Ph00_Chick_SSub_Marv'};
                    
    CSteps_DataFiles = {'04_08_2021_Ph01_Chick_CStep_Marv',...  % 17 August Datasets (FFFNoise, FFF2, CStep, Chirp3 clustering) - 22,11,2021
                      '04_08_2021_Ph02_Chick_CStep_Marv',...
                      '05_08_2021_Ph00_Chick_CStep_Marv',...
                      '05_08_2021_Ph01_Chick_CStep_Marv',...
                      '06_08_2021_Ph00_Chick_CStep_Marv',...
                      '06_08_2021_2nd_Ph00_Chick_CStep_Marv',...
                      '11_08_2021_Ph00_Chick_CStep_Marv',...
                      '12_08_2021_Ph00_Chick_CStep_Marv',...
                      '12_08_2021_Ph02_Chick_CStep_Marv',...
                      '13_08_2021_Ph00_Chick_CStep_Marv',...
                      '13_08_2021_Ph01_Chick_CStep_Marv',...
                      '14_08_2021_Ph00_Chick_CStep_Marv',...
                      '17_08_2021_Ph00_Chick_CStep_Marv',...
                      '19_08_2021_Ph00_Chick_CStep_Marv',...
                      '19_08_2021_Ph01_Chick_CStep_Marv',...
                      '20_08_2021_Ph00_Chick_CStep_Marv',...
                      '21_08_2021_Ph00_Chick_CStep_Marv'};
    
    All_DataFiles = [All_DataFiles;CSteps_DataFiles];
end
if p.Obj_clust_vec(12) == 1 || p.Obj_plot_vec(14) == 1 % Chirp3
    %%% Chicken data
%     Chirp3_DataFiles = {'03_08_2021_Ph01_Chick_Chirp3_Marv'};      % 03_08_21 (Phase_01) - August Data 1
%     Chirp3_DataFiles = {'04_08_2021_Ph01_Chick_Chirp3_Marv'};      % 04_08_21 (Phase_01) - August Data 2
%     Chirp3_DataFiles = {'05_08_2021_Ph01_Chick_Chirp3_Marv'};      % 05_08_21 (Phase_01) - August Data 3
%     Chirp3_DataFiles = {'06_08_2021_Ph00_Chick_Chirp3_Marv'};      % 06_08_21 (Phase_00) - August Data 4
%     Chirp3_DataFiles = {'12_08_2021_Ph02_Chick_Chirp3_Marv'};      % 12_08_21 (Phase_02) - August Data 5
%     Chirp3_DataFiles = {'13_08_2021_Ph00_Chick_Chirp3_Marv'};      % 13_08_21 (Phase_00) - August Data 6

%     Chirp3_DataFiles = {'03_08_2021_Ph01_Chick_Chirp3_Marv',...   % August Data 1-6 
%                         '04_08_2021_Ph01_Chick_Chirp3_Marv',...
%                         '05_08_2021_Ph01_Chick_Chirp3_Marv',...
%                         '06_08_2021_Ph00_Chick_Chirp3_Marv',...
%                         '12_08_2021_Ph02_Chick_Chirp3_Marv',...
%                         '13_08_2021_Ph00_Chick_Chirp3_Marv'};

    Chirp3_DataFiles = {'04_08_2021_Ph01_Chick_Chirp3_Marv',...  % 17 August Datasets (FFFNoise, FFF2, CStep, Chirp3 clustering) - 22,11,2021
                      '04_08_2021_Ph02_Chick_Chirp3_Marv',...
                      '05_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '05_08_2021_Ph01_Chick_Chirp3_Marv',...
                      '06_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '06_08_2021_2nd_Ph00_Chick_Chirp3_Marv',...
                      '11_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '12_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '12_08_2021_Ph02_Chick_Chirp3_Marv',...
                      '13_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '13_08_2021_Ph01_Chick_Chirp3_Marv',...
                      '14_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '17_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '19_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '19_08_2021_Ph01_Chick_Chirp3_Marv',...
                      '20_08_2021_Ph00_Chick_Chirp3_Marv',...
                      '21_08_2021_Ph00_Chick_Chirp3_Marv'};
    
    All_DataFiles = [All_DataFiles;Chirp3_DataFiles];
end

% Ref as: FFF_DataFiles{1}.
% Ref as: All_DataFiles{1,2}.
%Gen_Data_Names_vec      = {'data 1','data 2','data 3','data 4','data 5','data 6','data 7'};
%Gen_Data_Names_bar_vec  = categorical({'data 1','data 2','data 3','data 4','data 5','data 6','data 7'});
% Gen_Data_Names_vec      = {'data 1','data 2','data 3','data 4','data 5','data 6'};
% Gen_Data_Names_bar_vec  = categorical({'data 1','data 2','data 3','data 4','data 5','data 6'});
% Gen_Data_Names_vec      = {'data 1','data 2','data 3','data 4'};
% Gen_Data_Names_bar_vec  = categorical({'data 1','data 2','data 3','data 4'});
% Gen_Data_Names_vec      = {'data 1','data 2','data 3','data 4','data 5a','data 5b'};
% Gen_Data_Names_bar_vec  = categorical({'data 1','data 2','data 3','data 4','data 5a','data 5b'});
%Gen_Data_Names_vec      = {'data 1','data 2'};
%Gen_Data_Names_bar_vec  = categorical({'data 1','data 2'});
% Gen_Data_Names_vec      = {'data 1'};
% Gen_Data_Names_bar_vec  = categorical({'data 1'});
%%% August Data Sets 1-6
% Gen_Data_Names_vec      = {'3/8/21 Ph01','4/8/21 Ph01','5/8/21 Ph01','6/8/21 Ph00','12/08/21 Ph02','13/08/21 Ph00'};
% Gen_Data_Names_bar_vec  = categorical({'3/8/21 Ph01','4/8/21 Ph01','5/8/21 Ph01','6/8/21 Ph00','12/08/21 Ph02','13/08/21 Ph00'});
% Gen_Data_Names_bar_vec = reordercats(Gen_Data_Names_bar_vec,{'3/8/21 Ph01','4/8/21 Ph01','5/8/21 Ph01','6/8/21 Ph00','12/08/21 Ph02','13/08/21 Ph00'});
%%% August Data Sets 2-6 (e.g. for CSteps)
% Gen_Data_Names_vec      = {'4/8/21 Ph01','5/8/21 Ph01','6/8/21 Ph00','12/08/21 Ph02','13/08/21 Ph00'};
% Gen_Data_Names_bar_vec  = categorical({'4/8/21 Ph01','5/8/21 Ph01','6/8/21 Ph00','12/08/21 Ph02','13/08/21 Ph00'});
% Gen_Data_Names_bar_vec = reordercats(Gen_Data_Names_bar_vec,{'4/8/21 Ph01','5/8/21 Ph01','6/8/21 Ph00','12/08/21 Ph02','13/08/21 Ph00'});
%%% 16 August Data Sets for FFF2, 04,11,2021
% Gen_Data_Names_vec      = {'3/8/21 Ph01','4/8/21 Ph01','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/08/21 Ph00','12/08/21 Ph02','13/08/21 Ph00','13/08/21 Ph01','14/08/21 Ph00','17/08/21 Ph00','18/08/21 Ph00','19/08/21 Ph00','20/08/21 Ph00','21/08/21 Ph00'};
% Gen_Data_Names_bar_vec  = categorical({'3/8/21 Ph01','4/8/21 Ph01','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/08/21 Ph00','12/08/21 Ph02','13/08/21 Ph00','13/08/21 Ph01','14/08/21 Ph00','17/08/21 Ph00','18/08/21 Ph00','19/08/21 Ph00','20/08/21 Ph00','21/08/21 Ph00'});
% Gen_Data_Names_bar_vec = reordercats(Gen_Data_Names_bar_vec,{'3/8/21 Ph01','4/8/21 Ph01','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/08/21 Ph00','12/08/21 Ph02','13/08/21 Ph00','13/08/21 Ph01','14/08/21 Ph00','17/08/21 Ph00','18/08/21 Ph00','19/08/21 Ph00','20/08/21 Ph00','21/08/21 Ph00'});
%%% 20 August Data Sets for FFF2, 15,11,2021
% Gen_Data_Names_vec      = {'3/8/21 Ph01','4/8/21 Ph01','4/8/21 Ph02','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','12/8/21 Ph02','13/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','18/8/21 Ph00','18/8/21 Ph01','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'};
% Gen_Data_Names_bar_vec  = categorical({'3/8/21 Ph01','4/8/21 Ph01','4/8/21 Ph02','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','12/8/21 Ph02','13/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','18/8/21 Ph00','18/8/21 Ph01','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'});
% Gen_Data_Names_bar_vec = reordercats(Gen_Data_Names_bar_vec,{'3/8/21 Ph01','4/8/21 Ph01','4/8/21 Ph02','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','12/8/21 Ph02','13/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','18/8/21 Ph00','18/8/21 Ph01','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'});
%%% 12 August Data Sets for FFFNoise, 18,11,2021
% Gen_Data_Names_vec      = {'4/8/21 Ph02','5/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'};
% Gen_Data_Names_bar_vec  = categorical({'4/8/21 Ph02','5/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'});
% Gen_Data_Names_bar_vec = reordercats(Gen_Data_Names_bar_vec,{'4/8/21 Ph02','5/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'});
%%% 17 August Data Sets for FFFNoise, FFF2, CStep, Chirp3, 22,11,2021
Gen_Data_Names_vec      = {'4/8/21 Ph01','4/8/21 Ph02','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','12/8/21 Ph02','13/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'};
Gen_Data_Names_bar_vec  = categorical({'4/8/21 Ph01','4/8/21 Ph02','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','12/8/21 Ph02','13/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'});
Gen_Data_Names_bar_vec = reordercats(Gen_Data_Names_bar_vec,{'4/8/21 Ph01','4/8/21 Ph02','5/8/21 Ph00','5/8/21 Ph01','6/8/21 Ph00','6/8/21 Ph00 2nd','11/8/21 Ph00','12/8/21 Ph00','12/8/21 Ph02','13/8/21 Ph00','13/8/21 Ph01','14/8/21 Ph00','17/8/21 Ph00','19/8/21 Ph00','19/8/21 Ph01','20/8/21 Ph00','21/8/21 Ph00'});


Num_Stim_Loaded = size(All_DataFiles,1); % Num stim loaded over union of clustered and plotted stimuli.

%%% Plot Object Names
p.Plot_obj_name_vec = {'full RF positions','mean RF shape','FFF','chirp',...
                       'FFF noise','gratings','full RF size','full RF ellipticity',...
                       'full RF dom. axis angle','FFF2','chirp2','silent sub.','c. steps','chirp3','data set'};

%%% Create Experiment Identitly Vec & Record Number of Cells from Each Experiment
%First_stim = find(p.Obj_clust_vec==1,1,'first');

%data_set_indices

% Make dataset names into array to call using data_set_indices
% (expereiments) and First_data_set (stimulus)

Num_Cell_Per_Exp_vec = NaN(p.Num_data_sets,1);
%if sum([p.Obj_clust_vec(1:4),p.Obj_clust_vec(8:11),p.Obj_plot_vec(3:6),p.Obj_plot_vec(10:13)]) > 0 % sum([p.Obj_clust_vec(1:4),p.Obj_plot_vec(3:6)]) > 0 % PAR Mod 14,06,2021
%if sum([p.Obj_clust_vec(1:4),p.Obj_plot_vec(3:6)]) > 0 || (sum([p.Obj_clust_vec(5:7),p.Obj_plot_vec(7:9)])==0 && sum([p.Obj_clust_vec(8:11),p.Obj_plot_vec(10:13)])>0) % PAR Mod 25,06,2021
if sum([p.Obj_clust_vec(1:4),p.Obj_plot_vec(3:6)]) > 0 || (sum([p.Obj_clust_vec(5:7),p.Obj_plot_vec(7:9)])==0 && sum([p.Obj_clust_vec(8:12),p.Obj_plot_vec(10:14)])>0) % PAR Mod 23,09,2021
    for i = 1:p.Num_data_sets
        load(All_DataFiles{1,i},'spiketimestamps'); %All_DataFiles{First_stim,i}
        Num_Cell_Per_Exp_vec(i) = size(spiketimestamps,2);
        clear spiketimestamps;
    end
else % sum([p.Obj_clust_vec(1:4),p.Obj_clust_vec(8:11),p.Obj_plot_vec(3:6),p.Obj_plot_vec(10:13)]) = 0 was sum([p.Obj_clust_vec(1:4),p.Obj_plot_vec(3:6)]) = 0 % PAR Mod 14,06,2021
    for i = 1:p.Num_data_sets
        load(All_DataFiles{1,i},'spike_times_mat'); %All_DataFiles{First_stim,i}
        Num_Cell_Per_Exp_vec(i) = size(spike_times_mat,2);
        clear spike_times_mat;
    end
end
Total_Num_Cells = sum(Num_Cell_Per_Exp_vec);

% Final number of cells per stimulus vector (after poor quality data removed)
cl_var.Num_Cell_Per_Stim_vec = [];

% Cell identity vector
cl_var.Cell_Ident_vec = (1:1:Total_Num_Cells)';

% Overall Quality Vector (Across all Stimuli)
Overall_Quality_vec = logical(zeros(Total_Num_Cells,1));

% Experiment identity vector (which cell belongs to which experiment)
cl_var.Exp_Ident_vec = NaN(Total_Num_Cells,1);
loop_var = 0;
for i  = 1:p.Num_data_sets
    cl_var.Exp_Ident_vec(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = i*ones(Num_Cell_Per_Exp_vec(i),1);
    loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
end
clear loop_var;

%%% Dimensionality reduction options
if sum(p.Obj_clust_vec([1:4,8:12]))>0 % was 8:11 24,09,2021
    
    %%% Choose dimensionality reduction method
    % 1. PCA global           (across whole object);
    % 2. PCA segmented        (performed separately on each object segment);
    % 3. sparse PCA global    (not currently implemented);
    % 4. sparse PCA segmented (not currently implemented).
    p.Dim_red_meth = 2;
    
    %%% Choose PCA explained variance threshold type
    % 1. percentage variance explained;
    % 2. number of components.
    p.PCA_thresh_type = 1; % 1
    
    %%% Chose PCA explained variance threshold / number of components
    if p.PCA_thresh_type == 1
        p.PCA_ExpVar_thresh = 50; % 90, 50
    else % p.PCA_thresh_type == 2
        p.PCA_comp_num      = 2;  % 1, 5
    end
    
end

%%% NB: I could make the PCA options more complicated, but will leave for now

%%% Choose whether to scale each All_Scores column
% 1 = yes;
% 2 = no.
p.Scale_Scores = 2; % 1

%%% Choose clustering method
% 1. Gaussian Mixture Models;
% 2. Hierarchical.
p.Clust_Meth = 1;

%%% Choose clustering method settings
if p.Clust_Meth == 1     % Gaussian Mixture Models;
    
    %%% GMM Settings
    p.k_vec               = 1:50;                      % numbers of clusters 1:100, 1:30 usually
    p.Sigma               = {'diagonal','full'};       % diagonal and full covariance matrices {'diagonal','full'}
    p.SharedCovariance    = {true,false};              % covariance matrices are either all identical or allowed to differ {true,false}
    p.SCtext              = {'true','false'};          % text for labelling shared covariance option in figures {'true','false'}
    p.RegularizationValue = 1e-5;                      % to avoid badly conditioned covariance matrices. Matlab help: 0.01, Baden et al. (2016): 1e-5.
    p.Max_Iter            = 1e4;                       % max num iter. Default = 100, Matlab help: 10000.
    p.Replicate_Num       = 20;                        % 20 replicates in Baden et al. (2016).
    p.nK                  = numel(p.k_vec);            % number of clusters examined
    p.nSigma              = numel(p.Sigma);            % 2 covariance matrix types examined
    p.nSC                 = numel(p.SharedCovariance); % 2 covarience options
    
    % Use BIC: 1, AIC: 2, or prescribe number of clusters: 3.
    p.Info_crit     = 3; % 1
    p.Num_Clus      = 25;  % For when prescribe cluster number. Number of clusters.
    p.Sig_opt       = 1;   % For when prescribe cluster number. Diagonal: 1 or Full: 2.
    p.Sig_share_opt = 2;   % For when prescribe cluster number. Shared: 1 or Unshared: 2.
    
    % Set Minimum number of ROIs in a Cluster
    p.MinClusROI_thresh = 1; % 1, 10
    
else % p.Clust_Meth == 2 % Hierarchical
    
    %%% Use inconsistency cutoff: 1, distance cutoff: 2, or prescribe number of clusters: 3.
    p.Info_crit     = 2;
    if     p.Info_crit == 1 % inconsistency cutoff
        p.inconsistency_cutoff = 1.165;
    elseif p.Info_crit == 2 % distance cutoff
        p.distance_cutoff = 650; % 20, 19, 12.5 % 3: first data set all stim, 0.6: 2 data sets FFF, 1.2: 2 data sets FFF and Gratings
    else % p.Info_crit == 3 % prescribe number of clusters
        p.Num_Clus = 10; % Number of clusters
    end
    
    %%% pdist option (cityblock is good)
    % 1. euclidean (default);
    % 2. mahalanobis;
    % 3. cityblock;
    % 4. minkowski (default exponent = 2 (same as euclidean when = 2, same as city block when = 1, same as Chebychev when = inf));
    % 5. chebychev;
    % 6. correlation;
    % 7. spearman.
    pdist_opt_choice = 3;
    pdist_opt_vec    = {'euclidean','mahalanobis','cityblock','minkowski','chebychev','correlation','spearman'};
    p.pdist_opt      = pdist_opt_vec{pdist_opt_choice}; % 'cityblock'
    
    %%% Linkage option (average is good)
    % 1. single   (default);
    % 2. average;
    % 3. centroid (with euclidean dist only);
    % 4. complete;
    % 5. median   (with euclidean dist only);
    % 6. ward     (with euclidean dist only;
    % 7. weighted.
    linkage_opt_choice = 2;
    linkage_opt_vec    = {'single','average','centroid','complete','median','ward','weighted'};
    p.linkage_opt      = linkage_opt_vec{linkage_opt_choice}; % 'average'
    
    %%% Set Minimum number of ROIs in a Cluster
    p.MinClusROI_thresh = 1; % 1, 10
    
end

%%% Choose all cells or a subset
% 1. All cells;
% 2. Subset of cells.
Cell_Choice = 1; % 1
if Cell_Choice == 2
    First_Cell  = 1;
    Last_Cell   = 10;
    Cell_Choice_vec = First_Cell:Last_Cell;
end

%%% When calculating cell positions
if p.Obj_plot_vec(1) == 1
    p.Num_rows = 20;
    p.Num_cols = 20;
end


%%% When calculating STAs (Kernels)
if p.Obj_clust_vec(3) == 1 || p.Obj_plot_vec(5) == 1
    
    %%% Choose whether to work with time bins or define own time grid
    % 1. work with stimulus frames;
    % 2. define own time grid.
    p.Time_Choice = 2;
    
    %%% Choose whether to subtract the mean raw stimulus in the STA calculation
    % 1. don't subtract mean raw stim;
    % 2. subtract mean raw stim.
    p.STA_Choice  = 2; % 2 normally
    
    %%% Choose whether to calc mean raw stim from stimulus frames (much quicker)
    %%% or own time grid (much slower but a bit more accurate)
    %%% Always use same as stim frames/own time grid choice
    %%% (This is tech superfluous, but leave for now)
    % 1. Calc from stim frames;
    % 2. Calc from own time grid.
    if p.Time_Choice == 1
        p.Mean_Stim_Choice = 1;
    else % p.Time_Choice == 2
        p.Mean_Stim_Choice = 2;
    end
    
    %%% Choose length of time window in which to find the STEs
    if p.Time_Choice == 2 % define own time grid
        STE_int      = 1; % 0.25 sec = 250 ms % 07,09,2021 tried 0.35 for 05_08_21 (Phase_01) chick data - keep p.Num_STE_bins = 26; for now % PAR Mod 24,11,2021: Inc from 0.35 to 1 as Tom suggested
    end
    
    %%% Choose stimulus resolution (number of points sampled from stimulus time
    % window)
    p.Num_STE_bins = 51; % 5 % If p.Time_Choice = 1 then 5 or 10, if p.Time_Choice = 2 then STE_int*100 + 1. % PAR Mod 24,11,2021: Inc from 26 to 51
    
    %%% For STA normalisation (10,06,2021 Onwards)
    if p.Time_Choice == 1 % work with stimulus frames
        STE_norm_bin_start   = -40; % -40
        STE_norm_bin_end     = -21; % -21
        p.Num_STE_bins_norm  = STE_norm_bin_end - STE_norm_bin_start + 1; % 20
    else % if p.Time_Choice == 2 % define own time grid
        STE_norm_t_start = -2;     % -2 sec
        STE_norm_t_end   = -1;     % -1 sec
        p.Num_STE_bins_norm = 51; % 101 (10 ms intervals as with STE) % 18,11,2021 PAR Mod was 101, now 26 % PAR Mod 24,11,2021: Inc from 26 to 51
        STE_int_norm = STE_norm_t_end - STE_norm_t_start;
    end
    
end

%%% When calculating CNoise: Full RF Size
if p.Obj_clust_vec(5) == 1 || p.Obj_plot_vec(7) == 1
    
    %%% Choose between STA-SD and SC
    % 1. STA-SD;
    % 2. SC.
    RF_Ident_Method_RF_Size = 1;
    
end

%%% When calculating CNoise: Full RF Ellipticity
if p.Obj_clust_vec(6) == 1 || p.Obj_plot_vec(8) == 1
    
    %%% Choose between STA-SD and SC
    % 1. STA-SD;
    % 2. SC.
    RF_Ident_Method_RF_Ellipticity = 1;
    
    %%% Choose whether to threshold and if so at what value:
    % 1. Yes;
    % 2. No.
    p.RF_Ellip_thresh_choice = 1;
    if p.RF_Ellip_thresh_choice == 1
        p.RF_Ellip_thresh    = 20; % 20
    end
    
end

%%% When calculating CNoise: Full RF Dominant Axis Angle
if p.Obj_clust_vec(7) == 1 || p.Obj_plot_vec(9) == 1
    
    %%% Choose between STA-SD and SC
    % 1. STA-SD;
    % 2. SC.
    RF_Ident_Method_RF_Dom_Ax_Ang = 1;
    
end

if p.Obj_plot_vec(7) == 1 || p.Obj_plot_vec(8) == 1 || p.Obj_plot_vec(9) == 1 || p.Obj_plot_vec(15) == 1 % PAR Mod 23,09,2021  14-->15 % PAR Mod 22,09,2021 last was 'p.Obj_plot_vec(10) == 1', changed as have added new stim
    
    %%% Choose y lims on histograms
    % 1. same limits across all panels (for consistency);
    % 2. different limits on each panel (for visibility).
    p.hist_y_lim = 2;
    
end


%% Load Data

%%% Cell positions - needs to be rewritten to work with multiple data sets
%%% (not needed for now as retina unknown and diff orientation in each set)
if p.Obj_plot_vec(1) == 1
    
    cl_var.Full_RF_Pos_mat = NaN(max(Num_Cell_Per_Exp_vec),2);
    
    % Load data
    load(CNoise_20px_DataFiles{i},'STASD_Gaus_FullRF_Gaussian_mean_mat','True_cell_index_vec');%load('data_CNoise_20px_Paul_zf_RF_test_20_AllStat.mat','STASD_Gaus_FullRF_Gaussian_mean_mat');
    %cl_var.Full_RF_Pos_mat = STASD_Gaus_FullRF_Gaussian_mean_mat;
    %cl_var.Full_RF_Pos_mat(True_cell_index_vec,:) = STASD_Gaus_FullRF_Gaussian_mean_mat; % PAR temporaty mod 25,06,2021 for use with a single data set
    cl_var.Full_RF_Pos_mat = STASD_Gaus_FullRF_Gaussian_mean_mat;% PAR temporaty mod 25,06,2021 for use with a single data set
    clear STASD_Gaus_FullRF_Gaussian_mean_mat True_cell_index_vec;
    
    % Plot data
    figure;
    scatter(cl_var.Full_RF_Pos_mat(:,1),cl_var.Full_RF_Pos_mat(:,2),'x','LineWidth',1.5);
    axis equal;
    set(gca,'Ydir','reverse');
    xlim([0.5 (p.Num_rows + 0.5)]);
    ylim([0.5 (p.Num_cols + 0.5)]);
    xlabel('x');
    ylabel('y');
    title('cell positions');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
end



%%% Mean RF contour
if p.Obj_plot_vec(2) == 1
    
%     % Load data
%     load('data_CNoise_20px_Paul_zf_RF_test_20_AllStat.mat','');
%     cl_var.Full_RF_Pos_mat = ;
%     clear ;
%     
%     % Plot data
%     figure;
    
    
    
end



%%% FFF
if p.Obj_clust_vec(1) == 1 || p.Obj_plot_vec(3) == 1
    
    FFF_NumTrigPerRep     = 4; % User defined
    FFF_spike_times_arr   = cell(p.Num_data_sets,1);
    FFF_trig_times_mat    = cell(p.Num_data_sets,1);
    FFF_stim_end_time_vec = NaN(p.Num_data_sets,1);
    FFF_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    
    FFF_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    FFF_Plot_Row = ceil(p.Num_data_sets/FFF_Plot_Dim);
    FFF_Plot_Col = FFF_Plot_Dim;
    
    FFF_QI_vec                     = NaN(Total_Num_Cells,1);
    FFF_Qual_length_ksdensity_grid = 1000; % 100 coule use 1e3 (size(FFF_spike_times_arr{i},1))
    FFF_Qual_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best
    
    FFF_Trig_Fig      = figure;
    FFF_Trig_Dist_Fig = figure;
    FFF_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(FFF_DataFiles{i});
        %load FFF.mat;
        %load 03_12_2020_ZF_Ph01_FFF_100_Cells_Corr.mat;
        %load 03_12_2020_ZF_Ph01_FFF_All_Cells.mat;
        
        % Spike times
        if Cell_Choice == 1 % All cells
            FFF_spike_times_arr{i} = spiketimestamps;
            %FFF_Cell_Choice_vec    = 1:1:size(spiketimestamps,2);
            %FFF_True_cell_index_vec = FFF_Cell_Choice_vec(any(~isnan(FFF_spike_times_arr{i}))); % Record indices of responsive cells
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define FFF version)
            FFF_spike_times_arr{i} = spiketimestamps(:,Cell_Choice_vec);
            %FFF_True_cell_index_vec = Cell_Choice_vec(any(~isnan(FFF_spike_times_arr{i}))); % Record indices of responsive cells
        end
        
        % Remove NaN columns (non-responsive cells)?
        %FFF_spike_times_arr{i}(:,~any(~isnan(FFF_spike_times_arr{i}))) = [];
        %cl_var.FFF_True_Num_Cells = size(FFF_spike_times_arr{i},2);
        
        % Triggers
        FFF_trigCh_vec            = Ch_new.trigger_ch;
        if isa(FFF_trigCh_vec,'cell')
            FFF_trigCh_vec = cell2mat(FFF_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        FFF_min_trigCh_vec        = min(FFF_trigCh_vec);
        FFF_max_trigCh_vec        = max(FFF_trigCh_vec);
        FFF_trigThreshFac         = 0.05;
        FFF_trigHigh_vec          = double(FFF_trigCh_vec > FFF_min_trigCh_vec + FFF_trigThreshFac*(FFF_max_trigCh_vec-FFF_min_trigCh_vec));
        [~,FFF_trig_index_vec]    = findpeaks(FFF_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        FFF_trig_index_vec_length = length(FFF_trig_index_vec);
        
        FFF_NumTrigRep_mod = mod(FFF_trig_index_vec_length,FFF_NumTrigPerRep); % New: 28,05,2021
        
        if FFF_NumTrigRep_mod ~= 0 % New: 28,05,2021 (whole if statment (last line was already there w/o if statement))
            FFF_trig_index_vec = FFF_trig_index_vec(1:end-FFF_NumTrigRep_mod);
            FFF_trig_index_vec_length = length(FFF_trig_index_vec);
            
        else
            FFF_NumTrigRep = FFF_trig_index_vec_length/FFF_NumTrigPerRep;
        end
        
        FFF_sampling_freq = Ch_new.SamplingFrequency;
        FFF_sampling_int  = 1/FFF_sampling_freq;
        FFF_trig_t_vec    = (0:1:length(FFF_trigCh_vec)-1)*FFF_sampling_int;
        
        FFF_trig_times_vec_temp = FFF_trig_t_vec(FFF_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        FFF_trig_times_mat{i}   = FFF_trig_times_vec_temp - FFF_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(FFF_Trig_Fig);
        subplot(FFF_Plot_Row,FFF_Plot_Col,i);
        plot(FFF_trig_t_vec,FFF_trigCh_vec); hold on;
        plot(FFF_trig_times_vec_temp,4085*ones(1,FFF_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % FFF trigger channel
        if (i/FFF_Plot_Col > FFF_Plot_Row-1) || ((i/FFF_Plot_Col > FFF_Plot_Row-2)&&(i/FFF_Plot_Col < FFF_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Plot_Col)~=0)&&(mod(i,FFF_Plot_Col)>mod(p.Num_data_sets,FFF_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,FFF_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        % Find and plot the distribution of trigger intervals
        FFF_actual_stim_int      = diff(FFF_trig_times_mat{i});
        FFF_mean_stim_int_vec(i) = mean(FFF_actual_stim_int);
        FFF_stddev_stim_int      = std(FFF_actual_stim_int);
        figure(FFF_Trig_Dist_Fig);
        subplot(FFF_Plot_Row,FFF_Plot_Col,i);
        histogram(FFF_actual_stim_int); hold on;
        FFF_v1 = vline(FFF_mean_stim_int_vec(i),'r');
        FFF_v2 = vline(FFF_mean_stim_int_vec(i)-FFF_stddev_stim_int,'g');
        FFF_v3 = vline(FFF_mean_stim_int_vec(i)+FFF_stddev_stim_int,'g');
        if (i/FFF_Plot_Col > FFF_Plot_Row-1) || ((i/FFF_Plot_Col > FFF_Plot_Row-2)&&(i/FFF_Plot_Col < FFF_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Plot_Col)~=0)&&(mod(i,FFF_Plot_Col)>mod(p.Num_data_sets,FFF_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,FFF_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % FFF trigger int. dist.
        set(FFF_v1,'LineWidth',1.5);
        set(FFF_v2,'LineWidth',1.5);
        set(FFF_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        
        % I find 40 triggers - so 10 repeats.
        % Trig 2s red 2s dark, Trig 2s green 2s dark, Trig 2s blue 2s dark, Trig 2s UV 2s dark
        % Hence 4 sec between triggers, 16 sec for 1 repeat and 160 sec for 10 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cell's spikes into increasing time order.
        % Check:
        % Paul = [1 2;3 4];
        % Paul(Paul<3) = NaN
        % Paul = [Paul; 1 2]
        % Paul = sort(Paul,1)
        
        % 1. remove spikes before and after the stimulus;
        % Test Before
        % min(FFF_spike_times_arr{i},[],'all')
        % max(FFF_spike_times_arr{i},[],'all')
        FFF_spike_times_arr{i}(FFF_spike_times_arr{i}<=FFF_trig_times_mat{i}(1)) = NaN;
        FFF_spike_times_arr{i}(FFF_spike_times_arr{i}>=(FFF_trig_times_mat{i}(end) + FFF_mean_stim_int_vec(i))) = NaN;
        % Test After
        % min(FFF_spike_times_arr{i},[],'all')
        % max(FFF_spike_times_arr{i},[],'all')
        
        
        if FFF_trig_index_vec_length>FFF_NumTrigPerRep
            FFF_stim_end_time_vec(i) = FFF_trig_times_mat{i}(FFF_NumTrigPerRep+1); % period of a single repeat (Time of start of 2nd repeat (hence length of first repeat))
        else % FFF_trig_index_vec_length<=FFF_NumTrigPerRep
            FFF_stim_end_time_vec(i) = 16; % sec
        end
        
        % 1.5 Calculate Quality Indices 
        FFF_Qual_ksdensity_grid  = linspace(0,FFF_stim_end_time_vec(i),FFF_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(FFF_Qual_length_ksdensity_grid,FFF_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = FFF_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/FFF_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:FFF_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),FFF_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*FFF_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,FFF_Qual_ksdensity_grid,'Bandwidth',FFF_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                FFF_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                FFF_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % FFF_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(FFF_QI_vec);
        %figure;
        %hist(FFF_QI_vec);
        
        % Plot the distribution of QIs
        figure(FFF_QI_Dist_Fig);
        subplot(FFF_Plot_Row,FFF_Plot_Col,i);
        if i == 1
            Plot_vec_loop = FFF_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = FFF_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        FFF_QI_v1 = vline(FFF_QI_thresh,'r','QI thresh.');
        FFF_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/FFF_Plot_Col > FFF_Plot_Row-1) || ((i/FFF_Plot_Col > FFF_Plot_Row-2)&&(i/FFF_Plot_Col < FFF_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Plot_Col)~=0)&&(mod(i,FFF_Plot_Col)>mod(p.Num_data_sets,FFF_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,FFF_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(FFF_QI_v1,'LineWidth',1.5);
        set(FFF_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if FFF_trig_index_vec_length>FFF_NumTrigPerRep
           % FFF_stim_end_time_vec(i) =
           % FFF_trig_times_mat{i}(FFF_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            FFF_spike_times_arr{i}   = mod(FFF_spike_times_arr{i},FFF_stim_end_time_vec(i));
       % else % FFF_trig_index_vec_length<=FFF_NumTrigPerRep % --> commented as now use above
            %FFF_stim_end_time_vec(i) = 16; % sec % --> commented as now use above
        end
        % FFF_NumTrigPerRep
        % FFF_NumTrigRep
        
        % 3. sort each cell's spikes into increasing time order.
        FFF_spike_times_arr{i} = sort(FFF_spike_times_arr{i},1);
        
    end
    
    figure(FFF_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(FFF_Trig_Fig,'color','w');
    set(FFF_Trig_Dist_Fig,'color','w');
    set(FFF_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(FFF_QI_vec); hold on;
    FFF_QI_v1 = vline(FFF_QI_thresh,'r','QI thresh.');
    FFF_QI_v2 = vline(max(FFF_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('FFF QI - All Cells');
    set(FFF_QI_v1,'LineWidth',1.5);
    set(FFF_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    
    % Calculate FFF QI logical vec
    FFF_QI_logical_vec = FFF_QI_vec > FFF_QI_thresh;
    Num_FFF_Qual_Cells = sum(FFF_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_FFF_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,FFF_QI_logical_vec],[],2);
    
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1
    for j = 2:p.Num_data_sets
        for i = 1:FFF_NumTrigPerRep
            if i < FFF_NumTrigPerRep
                indices_loop = find((FFF_spike_times_arr{j}>=FFF_trig_times_mat{j}(i))&(FFF_spike_times_arr{j}<FFF_trig_times_mat{j}(i+1)));
                b_loop       = FFF_trig_times_mat{1}(i+1);
                b_dash_loop  = FFF_trig_times_mat{j}(i+1);
            else %i == FFF_NumTrigPerRep
                if length(FFF_trig_times_mat{j}) > FFF_NumTrigPerRep
                    indices_loop = find((FFF_spike_times_arr{j}>=FFF_trig_times_mat{j}(i))&(FFF_spike_times_arr{j}<=FFF_trig_times_mat{j}(i+1)));
                    b_dash_loop  = FFF_trig_times_mat{j}(i+1);
                else
                    indices_loop = find((FFF_spike_times_arr{j}>=FFF_trig_times_mat{j}(i))&(FFF_spike_times_arr{j}<=FFF_trig_times_mat{j}(i) + FFF_mean_stim_int_vec(j)));
                    b_dash_loop  = FFF_trig_times_mat{j}(i) + FFF_mean_stim_int_vec(j);
                end
                if length(FFF_trig_times_mat{1}) > FFF_NumTrigPerRep % was FFF_trig_times_mat{j} --> changed 08,04,2021
                    b_loop       = FFF_trig_times_mat{1}(i+1);
                else
                    b_loop       = FFF_trig_times_mat{1}(i) + FFF_mean_stim_int_vec(1);
                end
            end
            a_loop       = FFF_trig_times_mat{1}(i);
            a_dash_loop  = FFF_trig_times_mat{j}(i);
            FFF_spike_times_arr{j}(indices_loop) = a_loop + (FFF_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    FFF_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        FFF_num_spikes_vec(i) = size(FFF_spike_times_arr{i},1);
    end
    FFF_max_num_spikes = max(FFF_num_spikes_vec);
    
    cl_var.FFF_spike_times_mat = NaN(FFF_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.FFF_spike_times_mat(1:FFF_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = FFF_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
%    cl_var.FFF_spike_times_mat = ;                         % Just above.
%    cl_var.FFF_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
     cl_var.FFF_trig_times_vec  = FFF_trig_times_mat{1};    % As map subsequent data set times back to first data set times
     cl_var.FFF_stim_end_time   = FFF_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
    
end



%%% Chirp
if p.Obj_clust_vec(2) == 1 || p.Obj_plot_vec(4) == 1
    
    Chirp_NumTrigPerRep     = 2; % User defined
    Chirp_spike_times_arr   = cell(p.Num_data_sets,1);
    Chirp_trig_times_mat    = cell(p.Num_data_sets,1);
    Chirp_stim_end_time_vec = NaN(p.Num_data_sets,1);
    %Chirp_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    Chirp_last_stim_int_vec = NaN(p.Num_data_sets,1);
    
    Chirp_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    Chirp_Plot_Row = ceil(p.Num_data_sets/Chirp_Plot_Dim);
    Chirp_Plot_Col = Chirp_Plot_Dim;
    
    Chirp_QI_vec                     = NaN(Total_Num_Cells,1);
    Chirp_Qual_length_ksdensity_grid = 1000; % 100 coule use 1e3
    Chirp_Qual_ksdensity_bdwth       = 5*1e-2; % Tom said between 1e-2 and 1e-1 is best
    
    Chirp_Trig_Fig      = figure;
    Chirp_Trig_Dist_Fig = figure;
    Chirp_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(Chirp_DataFiles{i});
        %load Chirp.mat;
        %load 03_12_2020_ZF_Ph01_Chirp_1_100_Cells_Corr.mat;
        %load 03_12_2020_ZF_Ph01_Chirp_1_All_Cells.mat;
        
        % Spike times
        if Cell_Choice == 1 % All cells
            Chirp_spike_times_arr{i}   = spiketimestamps;
            %Chirp_Cell_Choice_vec      = 1:1:size(spiketimestamps,2);
            %Chirp_True_cell_index_vec = Chirp_Cell_Choice_vec(any(~isnan(Chirp_spike_times_arr{i}))); % Record indices of responsive cells
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define Chirp version)
            Chirp_spike_times_arr{i}   = spiketimestamps(:,Cell_Choice_vec);
            %Chirp_True_cell_index_vec = Cell_Choice_vec(any(~isnan(Chirp_spike_times_arr{i}))); % Record indices of responsive cells
        end
        
        % Remove NaN columns (non-responsive cells)?
        %Chirp_spike_times_arr{i}(:,~any(~isnan(Chirp_spike_times_arr{i}))) = [];
        %cl_var.Chirp_True_Num_Cells = size(Chirp_spike_times_arr{i},2);
        
        % Triggers
        Chirp_trigCh_vec            = Ch_new.trigger_ch;
        if isa(Chirp_trigCh_vec,'cell')
            Chirp_trigCh_vec = cell2mat(Chirp_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        Chirp_min_trigCh_vec        = min(Chirp_trigCh_vec);
        Chirp_max_trigCh_vec        = max(Chirp_trigCh_vec);
        Chirp_trigThreshFac         = 0.05;
        Chirp_trigHigh_vec          = double(Chirp_trigCh_vec > Chirp_min_trigCh_vec + Chirp_trigThreshFac*(Chirp_max_trigCh_vec-Chirp_min_trigCh_vec));
        [~,Chirp_trig_index_vec]    = findpeaks(Chirp_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        Chirp_trig_index_vec_length = length(Chirp_trig_index_vec);
        
        Chirp_NumTrigRep = Chirp_trig_index_vec_length/Chirp_NumTrigPerRep;
        
        Chirp_sampling_freq      = Ch_new.SamplingFrequency;
        Chirp_sampling_int       = 1/Chirp_sampling_freq ;
        Chirp_trig_t_vec         = (0:1:length(Chirp_trigCh_vec)-1)*Chirp_sampling_int;
        
        Chirp_trig_times_vec_temp   = Chirp_trig_t_vec(Chirp_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        Chirp_trig_times_mat{i} = Chirp_trig_times_vec_temp - Chirp_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(Chirp_Trig_Fig);
        subplot(Chirp_Plot_Row,Chirp_Plot_Col,i);
        plot(Chirp_trig_t_vec,Chirp_trigCh_vec); hold on;
        plot(Chirp_trig_times_vec_temp,4085*ones(1,Chirp_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % chirp trigger channel
        if (i/Chirp_Plot_Col > Chirp_Plot_Row-1) || ((i/Chirp_Plot_Col > Chirp_Plot_Row-2)&&(i/Chirp_Plot_Col < Chirp_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp_Plot_Col)~=0)&&(mod(i,Chirp_Plot_Col)>mod(p.Num_data_sets,Chirp_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,Chirp_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        
        % Find and plot the distribution of trigger intervals (not so relevant here
        % as trig times deliberately nonuniform)
        Chirp_actual_stim_int = diff(Chirp_trig_times_mat{i});
        Chirp_mean_stim_int   = mean(Chirp_actual_stim_int);
        Chirp_stddev_stim_int = std(Chirp_actual_stim_int);
        figure(Chirp_Trig_Dist_Fig);
        subplot(Chirp_Plot_Row,Chirp_Plot_Col,i);
        histogram(Chirp_actual_stim_int,100); hold on;
        Chirp_v1 = vline(Chirp_mean_stim_int,'r');
        Chirp_v2 = vline(Chirp_mean_stim_int-Chirp_stddev_stim_int,'g');
        Chirp_v3 = vline(Chirp_mean_stim_int+Chirp_stddev_stim_int,'g');
        if (i/Chirp_Plot_Col > Chirp_Plot_Row-1) || ((i/Chirp_Plot_Col > Chirp_Plot_Row-2)&&(i/Chirp_Plot_Col < Chirp_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp_Plot_Col)~=0)&&(mod(i,Chirp_Plot_Col)>mod(p.Num_data_sets,Chirp_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,Chirp_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % Gratings 400px trigger int. dist.
        set(Chirp_v1,'LineWidth',1.5);
        set(Chirp_v2,'LineWidth',1.5);
        set(Chirp_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        
        % I find 6 triggers - so 3 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cell's spikes into increasing time order.
        
        % 1. remove spikes before and after the stimulus;
        if Chirp_trig_index_vec_length>Chirp_NumTrigPerRep
            Chirp_last_stim_int_vec(i) = Chirp_trig_times_mat{i}(3)-Chirp_trig_times_mat{i}(2);
        else % Chirp_trig_index_vec_length<=Chirp_NumTrigPerRep
            Chirp_last_stim_int_vec(i) = 72; % sec
        end
        Chirp_spike_times_arr{i}(Chirp_spike_times_arr{i}<=Chirp_trig_times_mat{i}(1)) = NaN;
        Chirp_spike_times_arr{i}(Chirp_spike_times_arr{i}>=(Chirp_trig_times_mat{i}(end) + Chirp_last_stim_int_vec(i))) = NaN;
        
        
        if Chirp_trig_index_vec_length>Chirp_NumTrigPerRep
            Chirp_stim_end_time_vec(i) = Chirp_trig_times_mat{i}(Chirp_NumTrigPerRep+1); % period of a single repeat
        else % Chirp_trig_index_vec_length<=Chirp_NumTrigPerRep
            Chirp_stim_end_time_vec(i) = 77; % sec
        end
        
        % 1.5 Calculate Quality Indices 
        Chirp_Qual_ksdensity_grid  = linspace(0,Chirp_stim_end_time_vec(i),Chirp_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(Chirp_Qual_length_ksdensity_grid,Chirp_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = Chirp_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/Chirp_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:Chirp_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),Chirp_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*Chirp_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,Chirp_Qual_ksdensity_grid,'Bandwidth',Chirp_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                Chirp_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                Chirp_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % Chirp_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(Chirp_QI_vec);
        %figure;
        %hist(Chirp_QI_vec);
        
        % Plot the distribution of QIs
        figure(Chirp_QI_Dist_Fig);
        subplot(Chirp_Plot_Row,Chirp_Plot_Col,i);
        if i == 1
            Plot_vec_loop = Chirp_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = Chirp_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        Chirp_QI_v1 = vline(Chirp_QI_thresh,'r','QI thresh.');
        Chirp_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/Chirp_Plot_Col > Chirp_Plot_Row-1) || ((i/Chirp_Plot_Col > Chirp_Plot_Row-2)&&(i/Chirp_Plot_Col < Chirp_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp_Plot_Col)~=0)&&(mod(i,Chirp_Plot_Col)>mod(p.Num_data_sets,Chirp_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,Chirp_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(Chirp_QI_v1,'LineWidth',1.5);
        set(Chirp_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if Chirp_trig_index_vec_length>Chirp_NumTrigPerRep
            %Chirp_stim_end_time_vec(i) = Chirp_trig_times_mat{i}(Chirp_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            Chirp_spike_times_arr{i} = mod(Chirp_spike_times_arr{i},Chirp_stim_end_time_vec(i));
        %else % Chirp_trig_index_vec_length<=Chirp_NumTrigPerRep % --> commented as now use above
           %Chirp_stim_end_time_vec(i) = 77; % sec % --> commented as now use above
        end
        
        % 3. sort each cell's spikes into increasing time order.
        Chirp_spike_times_arr{i} = sort(Chirp_spike_times_arr{i},1);
        
    end
    
    figure(Chirp_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Chirp_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Chirp_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(Chirp_Trig_Fig,'color','w');
    set(Chirp_Trig_Dist_Fig,'color','w');
    set(Chirp_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(Chirp_QI_vec); hold on;
    Chirp_QI_v1 = vline(Chirp_QI_thresh,'r','QI thresh.');
    Chirp_QI_v2 = vline(max(Chirp_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('Chirp QI - All Cells');
    set(Chirp_QI_v1,'LineWidth',1.5);
    set(Chirp_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % Calculate Chirp QI logical vec
    Chirp_QI_logical_vec = Chirp_QI_vec > Chirp_QI_thresh;
    Num_Chirp_Qual_Cells = sum(Chirp_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_Chirp_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,Chirp_QI_logical_vec],[],2);
    
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1 ---> replace Chirp_mean_stim_int_vec(j) with
    %%% corrected as diff trig ints
    for j = 2:p.Num_data_sets
        for i = 1:Chirp_NumTrigPerRep
            if i < Chirp_NumTrigPerRep
                indices_loop = find((Chirp_spike_times_arr{j}>=Chirp_trig_times_mat{j}(i))&(Chirp_spike_times_arr{j}<Chirp_trig_times_mat{j}(i+1)));
                b_loop       = Chirp_trig_times_mat{1}(i+1);
                b_dash_loop  = Chirp_trig_times_mat{j}(i+1);
            else %i == Chirp_NumTrigPerRep
                if length(Chirp_trig_times_mat{j}) > Chirp_NumTrigPerRep
                    indices_loop = find((Chirp_spike_times_arr{j}>=Chirp_trig_times_mat{j}(i))&(Chirp_spike_times_arr{j}<=Chirp_trig_times_mat{j}(i+1)));
                    b_dash_loop  = Chirp_trig_times_mat{j}(i+1);
                else
                    indices_loop = find((Chirp_spike_times_arr{j}>=Chirp_trig_times_mat{j}(i))&(Chirp_spike_times_arr{j}<=Chirp_trig_times_mat{j}(i) + Chirp_last_stim_int_vec(j)));
                    b_dash_loop  = Chirp_trig_times_mat{j}(i) + Chirp_last_stim_int_vec(j);
                end
                if length(Chirp_trig_times_mat{1}) > Chirp_NumTrigPerRep % was Chirp_trig_times_mat{j} --> changed 08,04,2021
                    b_loop       = Chirp_trig_times_mat{1}(i+1);
                else
                    b_loop       = Chirp_trig_times_mat{1}(i) + Chirp_last_stim_int_vec(1);
                end
            end
            a_loop       = Chirp_trig_times_mat{1}(i);
            a_dash_loop  = Chirp_trig_times_mat{j}(i);
            Chirp_spike_times_arr{j}(indices_loop) = a_loop + (Chirp_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    Chirp_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        Chirp_num_spikes_vec(i) = size(Chirp_spike_times_arr{i},1);
    end
    Chirp_max_num_spikes = max(Chirp_num_spikes_vec);
    
    cl_var.Chirp_spike_times_mat = NaN(Chirp_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.Chirp_spike_times_mat(1:Chirp_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = Chirp_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
    %    cl_var.Chirp_spike_times_mat = ;                         % Just above.
    %    cl_var.Chirp_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
    cl_var.Chirp_trig_times_vec  = Chirp_trig_times_mat{1};    % As map subsequent data set times back to first data set times
    cl_var.Chirp_stim_end_time   = Chirp_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
    
end



%%% FFF_Noise
if p.Obj_clust_vec(3) == 1 || p.Obj_plot_vec(5) == 1
    
    %FFF_Noise_NumTrigPerRep     = #; % User defined
    %FFF_Noise_spike_times_arr   = cell(p.Num_data_sets,1);
    %FFF_Noise_trig_times_mat    = cell(p.Num_data_sets,1);
    %FFF_Noise_stim_end_time_vec = NaN(p.Num_data_sets,1);
    %FFF_Noise_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    
    FFF_Noise_STA_arr            = cell(p.Num_data_sets,1);
    FFF_Noise_STA_arr_norm       = cell(p.Num_data_sets,1); % For normalisation code (10,06,2021)
    FFF_Noise_STA_arr_scaled     = cell(p.Num_data_sets,1); % For normalisation code (10,06,2021)
    FFF_Noise_Full_STA_arr       = cell(p.Num_data_sets,1);
    FFF_Noise_Full_STA_arr_scaled = cell(p.Num_data_sets,1); % For normalisation code (10,06,2021)
    FFF_Noise_Full_t_mat         = cell(p.Num_data_sets,1); % Not sure if vecs are diff lengths for diff data sets (same for first two).
    FFF_Noise_Full_t_vec_int_vec = NaN(p.Num_data_sets,1);
    FFF_Noise_Full_t_vec_int_vec_norm = NaN(p.Num_data_sets,1); % For normalisation code (10,06,2021)
    FFF_Noise_Full_t_vec_end_vec = NaN(p.Num_data_sets,1);
    FFF_Noise_SD_mat             = NaN(Total_Num_Cells,4); % For normalisation code (10,06,2021) --> If want to keep track of which colours are significant
    
    FFF_Noise_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    FFF_Noise_Plot_Row = ceil(p.Num_data_sets/FFF_Noise_Plot_Dim);
    FFF_Noise_Plot_Col = FFF_Noise_Plot_Dim;
    
    FFF_Noise_QI_vec = NaN(Total_Num_Cells,1);
    %FFF_Noise_norm_mean_mat = NaN(Total_Num_Cells,4); % For normalisation code (10,06,2021) --> assume 4 colours
    %FFF_Noise_norm_SD_mat   = NaN(Total_Num_Cells,4); % For normalisation code (10,06,2021) --> assume 4 colours
    
    load('data_FFFNoiseStim.mat');
    
    FFF_Noise_Trig_Fig      = figure;
    FFF_Noise_Trig_Dist_Fig = figure;
    FFF_Noise_STA_Line      = figure;
    FFF_Noise_STA_HeatMap   = figure;
    FFF_Noise_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(FFF_Noise_DataFiles{i});
        %load FFF_Noise.mat;
        %load 03_12_2020_ZF_Ph01_FFF_Noise_100_Cells_Corr.mat;
        %load 03_12_2020_ZF_Ph01_FFF_Noise_All_Cells.mat;
        
        % Spike times
        if Cell_Choice == 1 % All cells
            FFF_Noise_spike_times_mat     = spiketimestamps;
            %FFF_Noise_Cell_Choice_vec     = 1:1:size(spiketimestamps,2);
            %FFF_Noise_True_cell_index_vec = FFF_Noise_Cell_Choice_vec(any(~isnan(FFF_Noise_spike_times_mat))); % Record indices of responsive cells
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define FFF Noise version)
            FFF_Noise_spike_times_mat     = spiketimestamps(:,Cell_Choice_vec);
            %FFF_Noise_True_cell_index_vec = Cell_Choice_vec(any(~isnan(FFF_Noise_spike_times_mat))); % Record indices of responsive cells
        end
        
        % Remove NaN columns (non-responsive cells)?
        %FFF_Noise_spike_times_mat(:,~any(~isnan(FFF_Noise_spike_times_mat))) = [];
        %FFF_Noise_True_Num_Cells = size(FFF_Noise_spike_times_mat,2);
        
        % Triggers
        FFF_Noise_trigCh_vec            = Ch_new.trigger_ch;
        if isa(FFF_Noise_trigCh_vec,'cell')
            FFF_Noise_trigCh_vec            = cell2mat(FFF_Noise_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        FFF_Noise_min_trigCh_vec        = min(FFF_Noise_trigCh_vec);
        FFF_Noise_max_trigCh_vec        = max(FFF_Noise_trigCh_vec);
        FFF_Noise_trigThreshFac         = 0.05;
        FFF_Noise_trigHigh_vec          = double(FFF_Noise_trigCh_vec > FFF_Noise_min_trigCh_vec + FFF_Noise_trigThreshFac*(FFF_Noise_max_trigCh_vec-FFF_Noise_min_trigCh_vec));
        [~,FFF_Noise_trig_index_vec]    = findpeaks(FFF_Noise_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        FFF_Noise_trig_index_vec_length = length(FFF_Noise_trig_index_vec);
        
        % FFF_Noise_NumTrigRep = FFF_Noise_trig_index_vec_length (no repetitions, just a series of triggers);
        
        FFF_Noise_sampling_freq      = Ch_new.SamplingFrequency;
        FFF_Noise_sampling_int       = 1/FFF_Noise_sampling_freq;
        FFF_Noise_trig_t_vec         = (0:1:length(FFF_Noise_trigCh_vec)-1)*FFF_Noise_sampling_int ;
        
        FFF_Noise_trig_times_vec_temp = FFF_Noise_trig_t_vec(FFF_Noise_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        FFF_Noise_trig_times_vec      = FFF_Noise_trig_times_vec_temp - FFF_Noise_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(FFF_Noise_Trig_Fig);
        subplot(FFF_Noise_Plot_Row,FFF_Noise_Plot_Col,i);
        plot(FFF_Noise_trig_t_vec,FFF_Noise_trigCh_vec); hold on;
        plot(FFF_Noise_trig_times_vec_temp,4085*ones(1,FFF_Noise_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % FFF Noise trigger channel
        if (i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-1) || ((i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-2)&&(i/FFF_Noise_Plot_Col < FFF_Noise_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Noise_Plot_Col)~=0)&&(mod(i,FFF_Noise_Plot_Col)>mod(p.Num_data_sets,FFF_Noise_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,FFF_Noise_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        % Find and plot the distribution of trigger intervals
        FFF_Noise_actual_stim_int = diff(FFF_Noise_trig_times_vec);
        FFF_Noise_mean_stim_int   = mean(FFF_Noise_actual_stim_int);
        FFF_Noise_stddev_stim_int = std(FFF_Noise_actual_stim_int);
        figure(FFF_Noise_Trig_Dist_Fig);
        subplot(FFF_Noise_Plot_Row,FFF_Noise_Plot_Col,i);
        histogram(FFF_Noise_actual_stim_int); hold on;
        FFF_Noise_v1 = vline(FFF_Noise_mean_stim_int,'r');
        FFF_Noise_v2 = vline(FFF_Noise_mean_stim_int-FFF_Noise_stddev_stim_int,'g');
        FFF_Noise_v3 = vline(FFF_Noise_mean_stim_int+FFF_Noise_stddev_stim_int,'g');
        if (i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-1) || ((i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-2)&&(i/FFF_Noise_Plot_Col < FFF_Noise_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Noise_Plot_Col)~=0)&&(mod(i,FFF_Noise_Plot_Col)>mod(p.Num_data_sets,FFF_Noise_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,FFF_Noise_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % FFF Noise trigger int. dist.
        set(FFF_Noise_v1,'LineWidth',1.5);
        set(FFF_Noise_v2,'LineWidth',1.5);
        set(FFF_Noise_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        % I find 3994 triggers - no repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. N/A (map spike in repeats into the first repeat);
        % 3. N/A (sort each cell's spikes into increasing time order).
        
        % 1. remove spikes before and after the stimulus;
        % Actually don't do this for FFF Noise as will remove spikes later when
        % finding the STAs and will want some of the spikes that are removed at the end here
        % FFF_Noise_spike_times_mat(FFF_Noise_spike_times_mat<=FFF_Noise_trig_times_vec(1)) = NaN;
        % FFF_Noise_spike_times_mat(FFF_Noise_spike_times_mat>=(FFF_Noise_trig_times_vec(end) + FFF_Noise_mean_stim_int)) = NaN;
        
        %%%
        
        % Stimulus data
        
        % A pop-up will ask for the number of columns and rows
        FFF_Noise_trig_per_frozen = FFF_Noise_trig_index_vec_length; % FFF_Noise_trig_index_vec_length; (tried 1200 when num was too big)
        % FFF_Noise_stimulus_arr    =
        % load_noise_from_hdf5_new(string('Noise_restored_logical1.h5'),true,1,double(FFF_Noise_trig_per_frozen)); % it turns out we should have been using the old data
        % FFF_Noise_stimulus_arr    = load_noise_from_hdf5(string('Noise.h5'),false,1,double(FFF_Noise_trig_per_frozen)); % was true, with false takes 1 box automatically given '1' chosen
        % FFF_Noise_stimulus_arr    = squeeze(FFF_Noise_stimulus_arr);
        % load('data_FFFNoiseStim.mat'); --> load above so just load once.
        FFF_Noise_stimulus_arr = [FFFNoiseStim_Red(1:double(FFF_Noise_trig_per_frozen)),FFFNoiseStim_Green(1:double(FFF_Noise_trig_per_frozen)),FFFNoiseStim_Blue(1:double(FFF_Noise_trig_per_frozen)),FFFNoiseStim_UV(1:double(FFF_Noise_trig_per_frozen))];
        % FFF_Noise_nr_boxes        = int64(sqrt(size(FFF_Noise_stimulus_arr,2)));
        % FFF_Noise_nr_colours      = int64(size(FFF_Noise_stimulus_arr,3));
        % FFF_Noise_stimulus_arr    = permute(FFF_Noise_stimulus_arr,[2 1 3]); % PAR Mod 01,02,2021
        % FFF_Noise_stimulus_arr    = reshape(FFF_Noise_stimulus_arr,[nr_boxes,nr_boxes,trig_per_frozen,nr_colours]);
        
        p.FFF_Noise_Spectral_Dim = size(FFF_Noise_stimulus_arr,2);
        
        p.stim_rows    = 1;
        p.stim_columns = 1;
        p.stim_frames  = FFF_Noise_trig_per_frozen;        % (size(stimulus_arr,3)) already defined as: FFF_Noise_trig_per_frozen
        p.Num_trigs    = FFF_Noise_trig_index_vec_length; % (length(trig_times_vec)) already defined as: FFF_Noise_trig_index_vec_length
        
        
        
        
        if p.Num_trigs > p.stim_frames % Frozen noise repeated case
            
            Num_FNoise_rep        = p.Num_trigs/p.stim_frames;
            p.Num_FNoise_rep_ceil = ceil(Num_FNoise_rep);
            max_trig_int          = max(diff(FFF_Noise_trig_times_vec(1:p.stim_frames)));
            inter_noise_int_vec   = NaN(p.Num_FNoise_rep_ceil-1,1);
            for k = 1:p.Num_FNoise_rep_ceil-1
                inter_noise_int_vec(k) = FFF_Noise_trig_times_vec(k*p.stim_frames+1)-FFF_Noise_trig_times_vec(k*p.stim_frames);
            end
            
            % Define variable to note if we are in the gaps or w/o gaps case
            % 1. w/o gaps;
            % 2. with gaps.
            if sum(inter_noise_int_vec <= max_trig_int) == p.Num_FNoise_rep_ceil-1
                p.Gap_Ind = 1;
            else % sum(inter_noise_int_vec <= max_trig_int) < p.Num_FNoise_rep_ceil-1 (anticipate = 0)
                p.Gap_Ind = 2;
            end
            
        else % p.Num_trigs <= p.stim_frames
            
            p.Gap_Ind = 2; % If there are no repeats then treat as if there are gaps
            
        end
        
        if p.Num_trigs < p.stim_frames
            p.Num_FNoise_rep_ceil = 1;
        end
        
        p.noise_length                 = min(p.Num_trigs,p.stim_frames);
        FFF_Noise_trig_times_vec_trunc = FFF_Noise_trig_times_vec(1:p.noise_length);
        Num_Trig_Final_NoiseChunk      = mod(p.Num_trigs,p.stim_frames);
        
        % Stimulus ppties
        stim_freq  = 20;          % Hz
        stim_int   = 1/stim_freq; % sec
        p.stim_int = stim_int;
        
        if p.Num_trigs >= p.Num_STE_bins+1 % PAR Mod 09,10,2020 (contents aren't new, just this outer if statement wrapping around it)
            if p.Time_Choice == 1 % stimulus frames
                
                p.stim_timesample_vec = linspace(-stim_int*p.Num_STE_bins,-stim_int,p.Num_STE_bins);
                
                % Mimimum time intervals after first and last triggers in which to consider
                % spikes (to ensure spikes occur following stimulus presentation across
                % full STE window)
                min_start_int = FFF_Noise_trig_times_vec(p.Num_STE_bins+1);  % sec (was 'stim_int*p.Num_STE_bins', this gave 0 indices when calc STA/STE_Full)
                min_end_int   = 2*stim_int;               % sec
                
            elseif p.Time_Choice == 2 % define own time grid
                
                % Pre-spike stimulus segment time vector
                p.stim_timesample_vec = linspace(-STE_int,0,p.Num_STE_bins);
                
                % Mimimum time intervals after first and last triggers in which to consider
                % spikes (to ensure spikes occur following stimulus presentation across
                % full STE window)
                min_start_int = STE_int;  % sec
                min_end_int   = stim_int; % sec
                
            end
        end
        
        % Eqivalent statement to above, but for normalisation code
        %if p.Num_trigs >= -STE_norm_bin_start+1 % STE_norm_bin_start only
        %occurs for p.Time_Choice == 1
            if p.Time_Choice == 1 % stimulus frames
                
                p.stim_timesample_vec_norm = linspace(stim_int*STE_norm_bin_start,stim_int*STE_norm_bin_end,p.Num_STE_bins_norm);
                
                % Mimimum time intervals after first and last triggers in which to consider
                % spikes (to ensure spikes occur following stimulus presentation across
                % full STE window)
                min_start_int_norm = FFF_Noise_trig_times_vec(-STE_norm_bin_start+1);
                min_end_int_norm   = 2*stim_int;               % sec
                
            elseif p.Time_Choice == 2 % define own time grid
                
                % Pre-spike stimulus segment time vector
                p.stim_timesample_vec_norm = linspace(STE_norm_t_start,STE_norm_t_end,p.Num_STE_bins_norm);
                
                % Mimimum time intervals after first and last triggers in which to consider
                % spikes (to ensure spikes occur following stimulus presentation across
                % full STE window)
                min_start_int_norm = -STE_norm_t_start;  % sec
                min_end_int_norm   = stim_int; % sec
                
            end
       % end
        
        
        FFF_Noise_first_trig_time      = FFF_Noise_trig_times_vec(1);
        FFF_Noise_last_trig_time       = FFF_Noise_trig_times_vec(end);
        FFF_Noise_last_trig_time_trunc = FFF_Noise_trig_times_vec_trunc(end);
        
        if p.STA_Choice  == 2 % subtract mean raw stim
            if p.Time_Choice == 1 || p.Mean_Stim_Choice  == 1 % stimulus frames || Calc from stim frames
                p.Num_Raw_Stim    = p.stim_frames - p.Num_STE_bins + 1;
                p.Num_Raw_Stim_norm = p.stim_frames - p.Num_STE_bins_norm + 1; % For normalisation code (10,06,2021)
                % Total number of raw stimuli displayed
            elseif p.Time_Choice == 2      % define own time grid
                p.Num_Raw_Stim    = ceil((FFF_Noise_last_trig_time_trunc - FFF_Noise_first_trig_time + stim_int - STE_int)/(STE_int/p.Num_STE_bins)) + 1;
                p.Num_Raw_Stim_norm = ceil((FFF_Noise_last_trig_time_trunc - FFF_Noise_first_trig_time + stim_int + STE_norm_t_start)/(STE_int_norm/p.Num_STE_bins_norm)) + 1; % For normalisation code (10,06,2021)
                % This is the total time over which stimuli are played 'last_trig_time-first_trig_time+stim_int'
                % minus the length of an STE 'STE_int' divided by the time between time
                % samples 'STE_int/p.Num_STE_bins' rounded upwards to the nearest integer, can times 10 to try to
                % make sure don't miss any stimulus combinations and plus 1 to count
                % the first stimulus before shifting the window.
                % 111.208542 seconds for model 1 without *10 so do this o/w takes too
                % long.
                p.Num_Raw_Stim_t_vec    = linspace(FFF_Noise_first_trig_time + STE_int,FFF_Noise_last_trig_time_trunc + stim_int,p.Num_Raw_Stim);
                p.Num_Raw_Stim_t_vec_norm = linspace(FFF_Noise_first_trig_time - STE_norm_t_start,FFF_Noise_last_trig_time_trunc + stim_int,p.Num_Raw_Stim_norm); % For normalisation code (10,06,2021)
            end
        end
        
        %p.Num_Raw_Stim
        %p.Num_FNoise_rep_ceil
        %p.Num_Raw_Stim_t_vec
        %p.stim_timesample_vec
        %p.Gap_Ind
        %p.length_spike_times
        %p.Spectral_Dim
        %p.noise_length
        %p.stim_int
        
        if p.STA_Choice == 2
            mean_raw_stim_arr      = NaN(p.Num_STE_bins,p.FFF_Noise_Spectral_Dim);
            mean_raw_stim_arr_norm = NaN(p.Num_STE_bins_norm,p.FFF_Noise_Spectral_Dim); % For normalisation code (10,06,2021)
            for k = 1:p.FFF_Noise_Spectral_Dim
                mean_raw_stim_arr(:,k)      = mean_raw_stim_FFF_Noise_fn(FFF_Noise_stimulus_arr(:,k),FFF_Noise_trig_times_vec_trunc,p);
                mean_raw_stim_arr_norm(:,k) = mean_raw_stim_FFF_Noise_norm_fn(FFF_Noise_stimulus_arr(:,k),FFF_Noise_trig_times_vec_trunc,p); % For normalisation code (10,06,2021)
            end
        end
        
        FFF_Noise_Full_t_vec_int_vec(i)      = p.stim_timesample_vec(2)-p.stim_timesample_vec(1);
        FFF_Noise_Full_t_vec_int_vec_norm(i) = p.stim_timesample_vec_norm(2)-p.stim_timesample_vec_norm(1); % For normalisation code (10,06,2021)
        if p.Time_Choice == 1
            Num_t_bins = p.Num_STE_bins*p.FFF_Noise_Spectral_Dim;
            FFF_Noise_Full_t_vec_end_vec(i)      = (Num_t_bins-1)*FFF_Noise_Full_t_vec_int_vec(i);
            FFF_Noise_Full_t_vec_end_vec_norm(i) = (Num_t_bins-1)*FFF_Noise_Full_t_vec_int_vec_norm(i); % For normalisation code (10,06,2021)
        else % p.Time_Choice == 2
            FFF_Noise_Full_t_vec_end_vec(i)      = -4*(p.stim_timesample_vec(1)) + 3*FFF_Noise_Full_t_vec_int_vec(i);
            FFF_Noise_Full_t_vec_end_vec_norm(i) = -4*(p.stim_timesample_vec_norm(1)) + 3*FFF_Noise_Full_t_vec_int_vec_norm(i); % For normalisation code (10,06,2021)
        end
        FFF_Noise_Full_t_mat{i}      = 0:FFF_Noise_Full_t_vec_int_vec(i):FFF_Noise_Full_t_vec_int_vec(i)*round(FFF_Noise_Full_t_vec_end_vec(i)/FFF_Noise_Full_t_vec_int_vec(i)); % this way of rounding will likely avoind the errors of the version above
        FFF_Noise_Full_t_mat_norm{i} = 0:FFF_Noise_Full_t_vec_int_vec_norm(i):FFF_Noise_Full_t_vec_int_vec_norm(i)*round(FFF_Noise_Full_t_vec_end_vec_norm(i)/FFF_Noise_Full_t_vec_int_vec_norm(i));  % For normalisation code (10,06,2021)
        FFF_Noise_STA_arr{i}         = NaN(Num_Cell_Per_Exp_vec(i),p.Num_STE_bins,p.FFF_Noise_Spectral_Dim); % NaN(FFF_Noise_True_Num_Cells,p.Num_STE_bins,p.FFF_Noise_Spectral_Dim)
        FFF_Noise_STA_arr_norm{i}    = NaN(Num_Cell_Per_Exp_vec(i),p.Num_STE_bins_norm,p.FFF_Noise_Spectral_Dim);
        
        tic
%         if p.Parpool == 1 % parallel processing - on
%             
%             parpool(p.Num_Cores);
%             
%             for k = 1:Num_Cell_Per_Exp_vec(i) % 1:FFF_Noise_True_Num_Cells
%                 
%                 spike_times_vec_loop = FFF_Noise_spike_times_mat(:,k);
%                 spike_times_vec_loop(isnan(spike_times_vec_loop)) = [];
%                 spike_times_vec_loop = spike_times_vec_loop(spike_times_vec_loop>(FFF_Noise_first_trig_time + min_start_int));
%                 if (Num_Trig_Final_NoiseChunk >= p.Num_STE_bins + 1) || (Num_Trig_Final_NoiseChunk == 0)
%                     spike_times_vec_loop     = spike_times_vec_loop(spike_times_vec_loop<(FFF_Noise_last_trig_time + min_end_int));
%                 else % Num_Trig_Final_NoiseChunk < p.Num_STE_bins + 1
%                     if p.Gap_Ind==1 % w/o gaps
%                         spike_times_vec_loop = spike_times_vec_loop(spike_times_vec_loop<(FFF_Noise_last_trig_time + min_end_int));
%                     else % p.Gap_Ind==2 with gaps
%                         spike_times_vec_loop = spike_times_vec_loop(spike_times_vec_loop<(FFF_Noise_trig_times_vec(p.stim_frames*(p.Num_FNoise_rep_ceil-1)) + min_end_int));
%                     end
%                 end
%                 length_spike_times_loop = length(spike_times_vec_loop);
%                 
%                 parfor j = 1:p.FFF_Noise_Spectral_Dim
%                     %STE_Full = NaN(p.Num_STE_bins,p.length_spike_times);
%                     STE_Full_loop               = STE_Full_FFF_Noise_fn(FFF_Noise_stimulus_arr(:,j),FFF_Noise_trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p);
%                     FFF_Noise_STA_arr{i}(k,:,j) = STA_FFF_Noise_fn(STE_Full_loop,mean_raw_stim_arr(:,j),p);%STA_fn_v2(STE_Full(:,:,j),mean_raw_stim_arr(:,j),p)
%                 end
%             end
%             
%         else % p.Parpool == 2 % parallel processing - off
            
            for k = 1:Num_Cell_Per_Exp_vec(i) % 1:FFF_Noise_True_Num_Cells
                
                spike_times_vec_loop = FFF_Noise_spike_times_mat(:,k);
                spike_times_vec_loop(isnan(spike_times_vec_loop)) = [];
                spike_times_vec_loop_norm = spike_times_vec_loop; % For normalisation code (10,06,2021)
                spike_times_vec_loop      = spike_times_vec_loop(spike_times_vec_loop>(FFF_Noise_first_trig_time + min_start_int));
                spike_times_vec_loop_norm = spike_times_vec_loop_norm(spike_times_vec_loop_norm>(FFF_Noise_first_trig_time + min_start_int_norm)); % For normalisation code (10,06,2021)
                if (Num_Trig_Final_NoiseChunk >= p.Num_STE_bins + 1) || (Num_Trig_Final_NoiseChunk == 0)
                    spike_times_vec_loop     = spike_times_vec_loop(spike_times_vec_loop<(FFF_Noise_last_trig_time + min_end_int));
                else % Num_Trig_Final_NoiseChunk < p.Num_STE_bins + 1
                    if p.Gap_Ind==1 % w/o gaps
                        spike_times_vec_loop = spike_times_vec_loop(spike_times_vec_loop<(FFF_Noise_last_trig_time + min_end_int));
                    else % p.Gap_Ind==2 with gaps
                        spike_times_vec_loop = spike_times_vec_loop(spike_times_vec_loop<(FFF_Noise_trig_times_vec(p.stim_frames*(p.Num_FNoise_rep_ceil-1)) + min_end_int));
                    end
                end
                length_spike_times_loop = length(spike_times_vec_loop);
                if (Num_Trig_Final_NoiseChunk >= p.Num_STE_bins_norm + 1) || (Num_Trig_Final_NoiseChunk == 0) % For normalisation code (10,06,2021) - whole if statement
                    spike_times_vec_loop_norm     = spike_times_vec_loop_norm(spike_times_vec_loop_norm<(FFF_Noise_last_trig_time + min_end_int_norm));
                else % Num_Trig_Final_NoiseChunk < p.Num_STE_bins_norm + 1
                    if p.Gap_Ind==1 % w/o gaps
                        spike_times_vec_loop_norm = spike_times_vec_loop_norm(spike_times_vec_loop_norm<(FFF_Noise_last_trig_time + min_end_int_norm));
                    else % p.Gap_Ind==2 with gaps
                        spike_times_vec_loop_norm = spike_times_vec_loop_norm(spike_times_vec_loop_norm<(FFF_Noise_trig_times_vec(p.stim_frames*(p.Num_FNoise_rep_ceil-1)) + min_end_int_norm));
                    end
                end
                length_spike_times_loop_norm = length(spike_times_vec_loop_norm); % For normalisation code (10,06,2021)
                
                for j = 1:p.FFF_Noise_Spectral_Dim
                    %STE_Full = NaN(p.Num_STE_bins,p.length_spike_times);
                    STE_Full_loop                    = STE_Full_FFF_Noise_fn(FFF_Noise_stimulus_arr(:,j),FFF_Noise_trig_times_vec_trunc,spike_times_vec_loop,length_spike_times_loop,p);
                    FFF_Noise_STA_arr{i}(k,:,j)      = STA_FFF_Noise_fn(STE_Full_loop,mean_raw_stim_arr(:,j),p);%STA_fn_v2(STE_Full(:,:,j),mean_raw_stim_arr(:,j),p)
                    STE_Full_loop_norm               = STE_Full_FFF_Noise_norm_fn(FFF_Noise_stimulus_arr(:,j),FFF_Noise_trig_times_vec_trunc,spike_times_vec_loop_norm,length_spike_times_loop_norm,p); % For normalisation code (10,06,2021)
                    FFF_Noise_STA_arr_norm{i}(k,:,j) = STA_FFF_Noise_fn(STE_Full_loop_norm,mean_raw_stim_arr_norm(:,j),p); % For normalisation code (10,06,2021)
                end
            end
            
        %end
        toc
        
        % Plot STA (kernels) for a single cell (lineplot)
        cell_choice_temp = 2;
        figure(FFF_Noise_STA_Line);
        subplot(FFF_Noise_Plot_Row,FFF_Noise_Plot_Col,i);
        plot(p.stim_timesample_vec,FFF_Noise_STA_arr{i}(cell_choice_temp,:,1),'r','LineWidth',1.5); hold on;
        plot(p.stim_timesample_vec,FFF_Noise_STA_arr{i}(cell_choice_temp,:,2),'g','LineWidth',1.5);
        plot(p.stim_timesample_vec,FFF_Noise_STA_arr{i}(cell_choice_temp,:,3),'b','LineWidth',1.5);
        plot(p.stim_timesample_vec,FFF_Noise_STA_arr{i}(cell_choice_temp,:,4),'m','LineWidth',1.5);
        xlim([p.stim_timesample_vec(1) 0]);
        if (i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-1) || ((i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-2)&&(i/FFF_Noise_Plot_Col < FFF_Noise_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Noise_Plot_Col)~=0)&&(mod(i,FFF_Noise_Plot_Col)>mod(p.Num_data_sets,FFF_Noise_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,FFF_Noise_Plot_Col) == 1
            ylabel('stimulus intensity (a.u.)');
        end
        title(Gen_Data_Names_vec{i}); % FFF Noise - STA
        set(gca,'FontSize',12);
        
        % Plot STAs (kernels) across all cells (heatmap)
       
        FFF_Noise_Full_STA_arr{i}        = reshape(FFF_Noise_STA_arr{i},[Num_Cell_Per_Exp_vec(i),p.Num_STE_bins*p.FFF_Noise_Spectral_Dim]); % reshape(FFF_Noise_STA_arr{i},[FFF_Noise_True_Num_Cells,p.Num_STE_bins*p.FFF_Noise_Spectral_Dim])
        
        figure(FFF_Noise_STA_HeatMap);
        subplot(FFF_Noise_Plot_Row,FFF_Noise_Plot_Col,i);
        imagesc(FFF_Noise_Full_t_mat{i}+0.5*FFF_Noise_Full_t_vec_int_vec(i),[],FFF_Noise_Full_STA_arr{i}); hold on; % cl_var.FFF_Noise_Full_t_vec+0.5*FFF_Noise_Full_t_vec_int_vec(i)
        colormap(FFF_Noise_STA_HeatMap,bluewhitered(256)); % gray(256), parula(256)
        colorbar;
        FFF_Noise_STA_v1 = vline(0.25*(FFF_Noise_Full_t_vec_end_vec(i)+FFF_Noise_Full_t_vec_int_vec(i)),'k');
        FFF_Noise_STA_v2 = vline(0.5*(FFF_Noise_Full_t_vec_end_vec(i)+FFF_Noise_Full_t_vec_int_vec(i)),'k');
        FFF_Noise_STA_v3 = vline(0.75*(FFF_Noise_Full_t_vec_end_vec(i)+FFF_Noise_Full_t_vec_int_vec(i)),'k');
        if (i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-1) || ((i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-2)&&(i/FFF_Noise_Plot_Col < FFF_Noise_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Noise_Plot_Col)~=0)&&(mod(i,FFF_Noise_Plot_Col)>mod(p.Num_data_sets,FFF_Noise_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,FFF_Noise_Plot_Col) == 1
            ylabel('cell');
        end
        title(Gen_Data_Names_vec{i}); % FFF Noise - STA
        set(FFF_Noise_STA_v1,'LineWidth',1.5);
        set(FFF_Noise_STA_v2,'LineWidth',1.5);
        set(FFF_Noise_STA_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        % Find means and SDs of norm STAs
%         FFF_Noise_norm_mean_mat_loop = NaN(Total_Num_Cells,4); % For normalisation code (10,06,2021)
%         FFF_Noise_norm_SD_mat_loop   = NaN(Total_Num_Cells,4); % For normalisation code (10,06,2021)
        %FFF_Noise_STA_arr_scaled     = cell(p.Num_data_sets,1);
        
        FFF_Noise_norm_mean_mat_loop = squeeze(mean(FFF_Noise_STA_arr_norm{i}(1:Num_Cell_Per_Exp_vec(i),:,:),2));
        FFF_Noise_norm_SD_mat_loop   = squeeze(std(FFF_Noise_STA_arr_norm{i}(1:Num_Cell_Per_Exp_vec(i),:,:),0,2));
        %FFF_Noise_STA_arr_scaled{i}  = FFF_Noise_STA_arr{i}; % subtract mean and divide SD for each cell for each colour
        FFF_Noise_STA_arr_scaled{i} = NaN(Num_Cell_Per_Exp_vec(i),p.Num_STE_bins,p.FFF_Noise_Spectral_Dim);
        for j=1:Num_Cell_Per_Exp_vec(i)
            for k = 1:p.FFF_Noise_Spectral_Dim
                FFF_Noise_STA_arr_scaled{i}(j,:,k) = (FFF_Noise_STA_arr{i}(j,:,k) - FFF_Noise_norm_mean_mat_loop(j,k))/FFF_Noise_norm_SD_mat_loop(j,k);
            end
        end
        
        FFF_Noise_Full_STA_arr_scaled{i} = reshape(FFF_Noise_STA_arr_scaled{i},[Num_Cell_Per_Exp_vec(i),p.Num_STE_bins*p.FFF_Noise_Spectral_Dim]); % For normalisation code (10,06,2021)
        
        %FFF_Noise_SD_mat             = NaN(Total_Num_Cells,4);
        %for j = 1:Num_Cell_Per_Exp_vec(i)   FFF_Noise_STA_arr_norm
        if i == 1
            % Old QI code
            %                 FFF_Noise_QI_vec(1:Num_Cell_Per_Exp_vec(1)) = max([std(FFF_Noise_STA_arr{1}(1:Num_Cell_Per_Exp_vec(1),:,1),0,2),...
            %                                                                    std(FFF_Noise_STA_arr{1}(1:Num_Cell_Per_Exp_vec(1),:,2),0,2),...
            %                                                                    std(FFF_Noise_STA_arr{1}(1:Num_Cell_Per_Exp_vec(1),:,3),0,2),...
            %                                                                    std(FFF_Noise_STA_arr{1}(1:Num_Cell_Per_Exp_vec(1),:,4),0,2)],[],2); % std(FFF_Noise_Full_STA_arr{i}(j,:))
            SD_mat_temp = squeeze(std(FFF_Noise_STA_arr_scaled{1},0,2));         % For normalisation code (10,06,2021)
            FFF_Noise_SD_mat(1:Num_Cell_Per_Exp_vec(1),:) = SD_mat_temp;
            FFF_Noise_QI_vec(1:Num_Cell_Per_Exp_vec(1)) = max(SD_mat_temp,[],2); % For normalisation code (10,06,2021)
        else % i>1
            % Old QI code
            %                 FFF_Noise_QI_vec(sum(Num_Cell_Per_Exp_vec(1:(i-1)))+1:sum(Num_Cell_Per_Exp_vec(1:i))) = max([std(FFF_Noise_STA_arr{i}(1:Num_Cell_Per_Exp_vec(i),:,1),0,2),...
            %                                                                                                              std(FFF_Noise_STA_arr{i}(1:Num_Cell_Per_Exp_vec(i),:,2),0,2),...
            %                                                                                                              std(FFF_Noise_STA_arr{i}(1:Num_Cell_Per_Exp_vec(i),:,3),0,2),...
            %                                                                                                              std(FFF_Noise_STA_arr{i}(1:Num_Cell_Per_Exp_vec(i),:,4),0,2)],[],2);
            %                                                                                                          %max([std(FFF_Noise_STA_arr{i}(sum(Num_Cell_Per_Exp_vec(1:(i-1)))+1:sum(Num_Cell_Per_Exp_vec(1:i)),:,1),0,2),...
            %                                                                                                              %std(FFF_Noise_STA_arr{i}(sum(Num_Cell_Per_Exp_vec(1:(i-1)))+1:sum(Num_Cell_Per_Exp_vec(1:i)),:,2),0,2),...
            %                                                                                                              %std(FFF_Noise_STA_arr{i}(sum(Num_Cell_Per_Exp_vec(1:(i-1)))+1:sum(Num_Cell_Per_Exp_vec(1:i)),:,3),0,2),...
            %                                                                                                              %std(FFF_Noise_STA_arr{i}(sum(Num_Cell_Per_Exp_vec(1:(i-1)))+1:sum(Num_Cell_Per_Exp_vec(1:i)),:,4),0,2)],[],2);
            SD_mat_temp = squeeze(std(FFF_Noise_STA_arr_scaled{i},0,2));         % For normalisation code (10,06,2021)
            FFF_Noise_SD_mat(sum(Num_Cell_Per_Exp_vec(1:(i-1)))+1:sum(Num_Cell_Per_Exp_vec(1:i)),:) = SD_mat_temp;
            FFF_Noise_QI_vec(sum(Num_Cell_Per_Exp_vec(1:(i-1)))+1:sum(Num_Cell_Per_Exp_vec(1:i)))   = max(SD_mat_temp,[],2); % For normalisation code (10,06,2021)
        end
        %end
        %FFF_Noise_Full_STA_arr{i} % Num_Cell_Per_Exp_vec(i),p.Num_STE_bins*p.FFF_Noise_Spectral_Dim
        %FFF_Noise_STA_arr = cell(p.Num_data_sets,1); FFF_Noise_STA_arr{i} = NaN(Num_Cell_Per_Exp_vec(i),p.Num_STE_bins,p.FFF_Noise_Spectral_Dim);
        
        
        % Plot the distribution of QIs
        figure(FFF_Noise_QI_Dist_Fig);
        subplot(FFF_Noise_Plot_Row,FFF_Noise_Plot_Col,i);
        if i == 1
            Plot_vec_loop = FFF_Noise_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = FFF_Noise_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        FFF_Noise_QI_v1 = vline(FFF_Noise_QI_thresh,'r','QI thresh.');
        FFF_Noise_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-1) || ((i/FFF_Noise_Plot_Col > FFF_Noise_Plot_Row-2)&&(i/FFF_Noise_Plot_Col < FFF_Noise_Plot_Row-1)&&(mod(p.Num_data_sets,FFF_Noise_Plot_Col)~=0)&&(mod(i,FFF_Noise_Plot_Col)>mod(p.Num_data_sets,FFF_Noise_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,FFF_Noise_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(FFF_Noise_QI_v1,'LineWidth',1.5);
        set(FFF_Noise_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);

        
    end
    
    %FFF_Noise_STA_Line
    %FFF_Noise_STA_HeatMap
    figure(FFF_Noise_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF Noise trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF_Noise_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF Noise trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF_Noise_STA_Line);
    annotation('textbox',[.5 .975 0 0],'String','FFF Noise - STA','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF_Noise_STA_HeatMap);
    annotation('textbox',[.5 .975 0 0],'String','FFF Noise - STA','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF_Noise_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF Noise QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(FFF_Noise_Trig_Fig,'color','w');
    set(FFF_Noise_Trig_Dist_Fig,'color','w');
    set(FFF_Noise_STA_Line,'color','w');
    set(FFF_Noise_STA_HeatMap,'color','w');
    set(FFF_Noise_QI_Dist_Fig,'color','w');
    
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(FFF_Noise_QI_vec); hold on;
    FFF_Noise_QI_v1 = vline(FFF_Noise_QI_thresh,'r','QI thresh.');
    FFF_Noise_QI_v2 = vline(max(FFF_Noise_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('FFF Noise QI - All Cells');
    set(FFF_Noise_QI_v1,'LineWidth',1.5);
    set(FFF_Noise_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % Calculate FFF Noise QI logical vec
    FFF_Noise_QI_logical_vec = FFF_Noise_QI_vec > FFF_Noise_QI_thresh;
    FFF_Noise_SD_logical_mat = FFF_Noise_SD_mat > FFF_Noise_QI_thresh; % For normalisation code (10,06,2021)
    Num_FFF_Noise_Qual_Cells = sum(FFF_Noise_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_FFF_Noise_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,FFF_Noise_QI_logical_vec],[],2);
    
    cl_var.FFF_Noise_STA             = NaN(Total_Num_Cells,p.Num_STE_bins,p.FFF_Noise_Spectral_Dim);
    cl_var.FFF_Noise_Full_STA        = NaN(Total_Num_Cells,p.Num_STE_bins*p.FFF_Noise_Spectral_Dim);
    cl_var.FFF_Noise_STA_scaled      = NaN(Total_Num_Cells,p.Num_STE_bins,p.FFF_Noise_Spectral_Dim); % For normalisation code (10,06,2021)
    cl_var.FFF_Noise_Full_STA_scaled = NaN(Total_Num_Cells,p.Num_STE_bins*p.FFF_Noise_Spectral_Dim); % For normalisation code (10,06,2021)
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.FFF_Noise_STA(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i),:,:)           = FFF_Noise_STA_arr{i};
        cl_var.FFF_Noise_Full_STA(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i),:)        = FFF_Noise_Full_STA_arr{i};
        cl_var.FFF_Noise_STA_scaled(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i),:,:)    = FFF_Noise_STA_arr_scaled{i};  % For normalisation code (10,06,2021)
        cl_var.FFF_Noise_Full_STA_scaled(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i),:) = FFF_Noise_Full_STA_arr_scaled{i}; % For normalisation code (10,06,2021)
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
    %cl_var.FFF_Noise_STA           = ; % Defined above
    %cl_var.FFF_Noise_Full_STA      = ; % Defined above
    cl_var.FFF_Noise_Full_t_vec     = FFF_Noise_Full_t_mat{1};
    cl_var.FFF_Noise_Full_t_vec_int = FFF_Noise_Full_t_vec_int_vec(1);
    cl_var.FFF_Noise_Full_t_vec_end = FFF_Noise_Full_t_vec_end_vec(1);
    
    clear FFFNoiseStim_Red FFFNoiseStim_Green FFFNoiseStim_Blue FFFNoiseStim_UV; % Remove original data files as no longer needed.
    
end



%%% Gratings_400px
if p.Obj_clust_vec(4) == 1 || p.Obj_plot_vec(6) == 1
    
    Gratings_400px_NumTrigPerRep     = 8; % User defined
    Gratings_400px_spike_times_arr   = cell(p.Num_data_sets,1);
    Gratings_400px_trig_times_mat    = cell(p.Num_data_sets,1);
    Gratings_400px_stim_end_time_vec = NaN(p.Num_data_sets,1);
    Gratings_400px_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    
    Gratings_400px_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    Gratings_400px_Plot_Row = ceil(p.Num_data_sets/Gratings_400px_Plot_Dim);
    Gratings_400px_Plot_Col = Gratings_400px_Plot_Dim;
    
    Gratings_400px_QI_vec                     = NaN(Total_Num_Cells,1);
    Gratings_400px_Qual_length_ksdensity_grid = 1000; % 100 coule use 1e3
    Gratings_400px_Qual_ksdensity_bdwth       = 1e-2; % Tom said 1e-2 is best (of options presented)
    
    Gratings_400px_Trig_Fig      = figure;
    Gratings_400px_Trig_Dist_Fig = figure;
    Gratings_400px_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(Gratings_400px_DataFiles{i});
        %load Gratings_400px.mat;
        %load 03_12_2020_ZF_Ph01_Gratings_100_Cells_Corr.mat;
        %load 03_12_2020_ZF_Ph01_Gratings_All_Cells.mat;
        
        % Spike times
        if Cell_Choice == 1 % All cells
            Gratings_400px_spike_times_arr{i} = spiketimestamps;
            %Gratings_400px_Cell_Choice_vec        = 1:1:size(spiketimestamps,2);
            %Gratings_400px_True_cell_index_vec    = Gratings_400px_Cell_Choice_vec(any(~isnan(cl_var.Gratings_400px_spike_times_mat))); % Record indices of responsive cells
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define Gratings_400px version)
            Gratings_400px_spike_times_arr{i} = spiketimestamps(:,Cell_Choice_vec);
            %Gratings_400px_True_cell_index_vec    = Cell_Choice_vec(any(~isnan(cl_var.Gratings_400px_spike_times_mat))); % Record indices of responsive cells
        end
        
        % Remove NaN columns (non-responsive cells)?
        %cl_var.Gratings_400px_spike_times_mat(:,~any(~isnan(cl_var.Gratings_400px_spike_times_mat))) = [];
        %cl_var.Gratings_400px_True_Num_Cells = size(cl_var.Gratings_400px_spike_times_mat,2);
        
        % Triggers
        Gratings_400px_trigCh_vec            = Ch_new.trigger_ch;
        if isa(Gratings_400px_trigCh_vec,'cell')
            Gratings_400px_trigCh_vec            = cell2mat(Gratings_400px_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        Gratings_400px_min_trigCh_vec        = min(Gratings_400px_trigCh_vec);
        Gratings_400px_max_trigCh_vec        = max(Gratings_400px_trigCh_vec);
        Gratings_400px_trigThreshFac         = 0.05;
        Gratings_400px_trigHigh_vec          = double(Gratings_400px_trigCh_vec > Gratings_400px_min_trigCh_vec + Gratings_400px_trigThreshFac*(Gratings_400px_max_trigCh_vec-Gratings_400px_min_trigCh_vec));
        [~,Gratings_400px_trig_index_vec]    = findpeaks(Gratings_400px_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        Gratings_400px_trig_index_vec_length = length(Gratings_400px_trig_index_vec);
        
        Gratings_400px_NumTrigRep = Gratings_400px_trig_index_vec_length/Gratings_400px_NumTrigPerRep;
        
        Gratings_400px_sampling_freq      = Ch_new.SamplingFrequency;
        Gratings_400px_sampling_int       = 1/Gratings_400px_sampling_freq;
        Gratings_400px_trig_t_vec         = (0:1:length(Gratings_400px_trigCh_vec)-1)*Gratings_400px_sampling_int ;
        
        Gratings_400px_trig_times_vec_temp = Gratings_400px_trig_t_vec(Gratings_400px_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        Gratings_400px_trig_times_mat{i}   = Gratings_400px_trig_times_vec_temp - Gratings_400px_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(Gratings_400px_Trig_Fig);
        subplot(Gratings_400px_Plot_Row,Gratings_400px_Plot_Col,i);
        plot(Gratings_400px_trig_t_vec,Gratings_400px_trigCh_vec); hold on;
        plot(Gratings_400px_trig_times_vec_temp,4085*ones(1,Gratings_400px_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % Gratings 400px trigger channel
        if (i/Gratings_400px_Plot_Col > Gratings_400px_Plot_Row-1) || ((i/Gratings_400px_Plot_Col > Gratings_400px_Plot_Row-2)&&(i/Gratings_400px_Plot_Col < Gratings_400px_Plot_Row-1)&&(mod(p.Num_data_sets,Gratings_400px_Plot_Col)~=0)&&(mod(i,Gratings_400px_Plot_Col)>mod(p.Num_data_sets,Gratings_400px_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,Gratings_400px_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        % Find and plot the distribution of trigger intervals
        Gratings_400px_actual_stim_int      = diff(Gratings_400px_trig_times_mat{i});
        Gratings_400px_mean_stim_int_vec(i) = mean(Gratings_400px_actual_stim_int);
        Gratings_400px_stddev_stim_int      = std(Gratings_400px_actual_stim_int);
        figure(Gratings_400px_Trig_Dist_Fig);
        subplot(Gratings_400px_Plot_Row,Gratings_400px_Plot_Col,i);
        histogram(Gratings_400px_actual_stim_int); hold on;
        Gratings_400px_v1 = vline(Gratings_400px_mean_stim_int_vec(i),'r');
        Gratings_400px_v2 = vline(Gratings_400px_mean_stim_int_vec(i)-Gratings_400px_stddev_stim_int,'g');
        Gratings_400px_v3 = vline(Gratings_400px_mean_stim_int_vec(i)+Gratings_400px_stddev_stim_int,'g');
        if (i/Gratings_400px_Plot_Col > Gratings_400px_Plot_Row-1) || ((i/Gratings_400px_Plot_Col > Gratings_400px_Plot_Row-2)&&(i/Gratings_400px_Plot_Col < Gratings_400px_Plot_Row-1)&&(mod(p.Num_data_sets,Gratings_400px_Plot_Col)~=0)&&(mod(i,Gratings_400px_Plot_Col)>mod(p.Num_data_sets,Gratings_400px_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,Gratings_400px_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % Gratings 400px trigger int. dist.
        set(Gratings_400px_v1,'LineWidth',1.5);
        set(Gratings_400px_v2,'LineWidth',1.5);
        set(Gratings_400px_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        % I find 48 triggers - so 6 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cells spike into increasing time order.
        
        % 1. remove spikes before and after the stimulus;
        Gratings_400px_spike_times_arr{i}(Gratings_400px_spike_times_arr{i}<=Gratings_400px_trig_times_mat{i}(1)) = NaN;
        Gratings_400px_spike_times_arr{i}(Gratings_400px_spike_times_arr{i}>=(Gratings_400px_trig_times_mat{i}(end) + Gratings_400px_mean_stim_int_vec(i))) = NaN;
        
        
        if Gratings_400px_trig_index_vec_length>Gratings_400px_NumTrigPerRep
            Gratings_400px_stim_end_time_vec(i) = Gratings_400px_trig_times_mat{i}(Gratings_400px_NumTrigPerRep+1); % period of a single repeat
        else % Gratings_400px_trig_index_vec_length<=Gratings_400px_NumTrigPerRep
            Gratings_400px_stim_end_time_vec(i) = 44; % sec
        end
        
        
        % 1.5 Calculate Quality Indices 
        Gratings_400px_Qual_ksdensity_grid  = linspace(0,Gratings_400px_stim_end_time_vec(i),Gratings_400px_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(Gratings_400px_Qual_length_ksdensity_grid,Gratings_400px_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = Gratings_400px_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/Gratings_400px_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:Gratings_400px_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),Gratings_400px_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*Gratings_400px_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,Gratings_400px_Qual_ksdensity_grid,'Bandwidth',Gratings_400px_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                Gratings_400px_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                Gratings_400px_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % Gratings_400px_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(Gratings_400px_QI_vec);
        %figure;
        %hist(Gratings_400px_QI_vec);
        
        % Plot the distribution of QIs
        figure(Gratings_400px_QI_Dist_Fig);
        subplot(Gratings_400px_Plot_Row,Gratings_400px_Plot_Col,i);
        if i == 1
            Plot_vec_loop = Gratings_400px_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = Gratings_400px_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        Gratings_400px_QI_v1 = vline(Gratings_400px_QI_thresh,'r','QI thresh.');
        Gratings_400px_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/Gratings_400px_Plot_Col > Gratings_400px_Plot_Row-1) || ((i/Gratings_400px_Plot_Col > Gratings_400px_Plot_Row-2)&&(i/Gratings_400px_Plot_Col < Gratings_400px_Plot_Row-1)&&(mod(p.Num_data_sets,Gratings_400px_Plot_Col)~=0)&&(mod(i,Gratings_400px_Plot_Col)>mod(p.Num_data_sets,Gratings_400px_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,Gratings_400px_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(Gratings_400px_QI_v1,'LineWidth',1.5);
        set(Gratings_400px_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if Gratings_400px_trig_index_vec_length>Gratings_400px_NumTrigPerRep
            %Gratings_400px_stim_end_time_vec(i) = Gratings_400px_trig_times_mat{i}(Gratings_400px_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            Gratings_400px_spike_times_arr{i}   = mod(Gratings_400px_spike_times_arr{i},Gratings_400px_stim_end_time_vec(i));
        %else % Gratings_400px_trig_index_vec_length<=Gratings_400px_NumTrigPerRep % --> commented as now use above
            %Gratings_400px_stim_end_time_vec(i) = 44; % sec % --> commented as now use above
        end
        
        % 3. sort each cells spike into increasing time order.
        Gratings_400px_spike_times_arr{i} = sort(Gratings_400px_spike_times_arr{i},1);
        
    end
    
    figure(Gratings_400px_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Gratings 400px trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Gratings_400px_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Gratings 400px trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Gratings_400px_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Gratings 400px QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(Gratings_400px_Trig_Fig,'color','w');
    set(Gratings_400px_Trig_Dist_Fig,'color','w');
    set(Gratings_400px_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(Gratings_400px_QI_vec); hold on;
    Gratings_400px_QI_v1 = vline(Gratings_400px_QI_thresh,'r','QI thresh.');
    Gratings_400px_QI_v2 = vline(max(Gratings_400px_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('Gratings 400px QI - All Cells');
    set(Gratings_400px_QI_v1,'LineWidth',1.5);
    set(Gratings_400px_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % Calculate Gratings_400px QI logical vec
    Gratings_400px_QI_logical_vec = Gratings_400px_QI_vec > Gratings_400px_QI_thresh;
    Num_Gratings_400px_Qual_Cells = sum(Gratings_400px_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_Gratings_400px_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,Gratings_400px_QI_logical_vec],[],2);
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1
    for j = 2:p.Num_data_sets
        for i = 1:Gratings_400px_NumTrigPerRep
            if i < Gratings_400px_NumTrigPerRep
                indices_loop = find((Gratings_400px_spike_times_arr{j}>=Gratings_400px_trig_times_mat{j}(i))&(Gratings_400px_spike_times_arr{j}<Gratings_400px_trig_times_mat{j}(i+1)));
                b_loop       = Gratings_400px_trig_times_mat{1}(i+1);
                b_dash_loop  = Gratings_400px_trig_times_mat{j}(i+1);
            else %i == Gratings_400px_NumTrigPerRep
                if length(Gratings_400px_trig_times_mat{j}) > Gratings_400px_NumTrigPerRep
                    indices_loop = find((Gratings_400px_spike_times_arr{j}>=Gratings_400px_trig_times_mat{j}(i))&(Gratings_400px_spike_times_arr{j}<=Gratings_400px_trig_times_mat{j}(i+1)));
                    b_dash_loop  = Gratings_400px_trig_times_mat{j}(i+1);
                else
                    indices_loop = find((Gratings_400px_spike_times_arr{j}>=Gratings_400px_trig_times_mat{j}(i))&(Gratings_400px_spike_times_arr{j}<=Gratings_400px_trig_times_mat{j}(i) + Gratings_400px_mean_stim_int_vec(j)));
                    b_dash_loop  = Gratings_400px_trig_times_mat{j}(i) + Gratings_400px_mean_stim_int_vec(j);
                end
                if length(Gratings_400px_trig_times_mat{1}) > Gratings_400px_NumTrigPerRep % was Gratings_400px_trig_times_mat{{j} --> changed 08,04,2021
                    b_loop       = Gratings_400px_trig_times_mat{1}(i+1);
                else
                    b_loop       = Gratings_400px_trig_times_mat{1}(i) + Gratings_400px_mean_stim_int_vec(1);
                end
            end
            a_loop       = Gratings_400px_trig_times_mat{1}(i);
            a_dash_loop  = Gratings_400px_trig_times_mat{j}(i);
            Gratings_400px_spike_times_arr{j}(indices_loop) = a_loop + (Gratings_400px_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    Gratings_400px_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        Gratings_400px_num_spikes_vec(i) = size(Gratings_400px_spike_times_arr{i},1);
    end
    Gratings_400px_max_num_spikes = max(Gratings_400px_num_spikes_vec);
    
    cl_var.Gratings_400px_spike_times_mat = NaN(Gratings_400px_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.Gratings_400px_spike_times_mat(1:Gratings_400px_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = Gratings_400px_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
%    cl_var.Gratings_400px_spike_times_mat = ;                         % Just above.
%    cl_var.Gratings_400px_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
     cl_var.Gratings_400px_trig_times_vec  = Gratings_400px_trig_times_mat{1};    % As map subsequent data set times back to first data set times
     cl_var.Gratings_400px_stim_end_time   = Gratings_400px_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
 
end


%%% CNoise: Full RF Size
if p.Obj_clust_vec(5) == 1 || p.Obj_plot_vec(7) == 1
    
    
    Full_RF_Size_mat                       = NaN(p.Num_data_sets,max(Num_Cell_Per_Exp_vec));
    %cl_var.Full_RF_Size_hist_BinEdges_mat = cell(p.Num_data_sets);
    %cl_var.Full_RF_Size_hist_NumBins_vec  = NaN(p.Num_data_sets,1);
    cl_var.Full_RF_Size_hist_BinWidth_vec  = NaN(p.Num_data_sets,1);
    cl_var.Full_RF_Size_hist_BinLimits_mat = NaN(p.Num_data_sets,2);
    
    RF_Size_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    RF_Size_Plot_Row = ceil(p.Num_data_sets/RF_Size_Plot_Dim);
    RF_Size_Plot_Col = RF_Size_Plot_Dim;
    
    RF_Size_Fig = figure;
    
    for i = 1:p.Num_data_sets
        
        % Load data
        if RF_Ident_Method_RF_Size == 1     % STA-SD
            %load('data_CNoise_20px_Paul_zf_RF_test_20_AllStat.mat','Num_STASD_ASP_FullRF_Pixels_vec');
            load(CNoise_20px_DataFiles{i},'Num_STASD_ASP_FullRF_Pixels_vec','True_cell_index_vec');%,'True_Num_Cells'
            Full_RF_Size_mat(i,True_cell_index_vec) = Num_STASD_ASP_FullRF_Pixels_vec;%Full_RF_Size_mat(i,1:Num_Cell_Per_Exp_vec(i)), int32(True_cell_index_vec)
            clear Num_STASD_ASP_FullRF_Pixels_vec True_cell_index_vec;
        else % RF_Ident_Method_RF_Size == 2 % SC
            %load('data_CNoise_20px_Paul_zf_RF_test_20_AllStat.mat','Num_SC_ASP_FullRF_Pixels_vec');
            load(CNoise_20px_DataFiles{i},'Num_SC_ASP_FullRF_Pixels_vec','True_cell_index_vec');
            Full_RF_Size_mat(i,1:True_cell_index_vec) = Num_SC_ASP_FullRF_Pixels_vec; % Full_RF_Size_mat(i,1:Num_Cell_Per_Exp_vec(i))
            clear Num_SC_ASP_FullRF_Pixels_vec True_cell_index_vec;
        end
        
        % Plot data
        figure(RF_Size_Fig);
        subplot(RF_Size_Plot_Row,RF_Size_Plot_Col,i);
        RF_Size_hist = histogram(Full_RF_Size_mat(i,1:Num_Cell_Per_Exp_vec(i)));
        xlim([RF_Size_hist.BinEdges(1) RF_Size_hist.BinEdges(end)]);
        ylim([0 max(RF_Size_hist.Values)]);
        if (i/RF_Size_Plot_Col > RF_Size_Plot_Row-1) || ((i/RF_Size_Plot_Col > RF_Size_Plot_Row-2)&&(i/RF_Size_Plot_Col < RF_Size_Plot_Row-1)&&(mod(p.Num_data_sets,RF_Size_Plot_Col)~=0)&&(mod(i,RF_Size_Plot_Col)>mod(p.Num_data_sets,RF_Size_Plot_Col)))
            xlabel('num. pixels');
        end
        if mod(i,RF_Size_Plot_Col) == 1
            ylabel('num. cells');
        end
        title(Gen_Data_Names_vec{i}); % RF size
        set(gca,'FontSize',12);
        
        %cl_var.Full_RF_Size_hist_BinEdges_mat{i}   = RF_Size_hist.BinEdges;
        %cl_var.Full_RF_Size_hist_NumBins_vec(i)    = RF_Size_hist.NumBins;
        cl_var.Full_RF_Size_hist_BinWidth_vec(i)    = RF_Size_hist.BinWidth;
        cl_var.Full_RF_Size_hist_BinLimits_mat(i,:) = RF_Size_hist.BinLimits;
        %cl_var.RF_Size_hist_MaxValue               = max(RF_Size_hist.Values); % Not required as calc max for each cluster
        
    end
    
    figure(RF_Size_Fig);
    annotation('textbox',[.5 .975 0 0],'String','RF size','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    set(RF_Size_Fig,'color','w');
    
    cl_var.Full_RF_Size_vec = NaN(Total_Num_Cells,1);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.Full_RF_Size_vec(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = Full_RF_Size_mat(i,1:Num_Cell_Per_Exp_vec(i));
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
    % Calculate Full RF Size QI logical vec
    Full_RF_Size_QI_logical_vec = cl_var.Full_RF_Size_vec > CNoise_QI_thresh;
    Num_Full_RF_Size_Qual_Cells = sum(Full_RF_Size_QI_logical_vec);
    
    % Update Overall Quality Vector
    if Use_CNoise_QI == 1
    Overall_Quality_vec = max([Overall_Quality_vec,Full_RF_Size_QI_logical_vec],[],2);  
    end
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_Full_RF_Size_Qual_Cells];
    
end


%%% CNoise: Full RF Ellipticity
if p.Obj_clust_vec(6) == 1 || p.Obj_plot_vec(8) == 1
    
    Full_RF_Ellipticity_mat                       = NaN(p.Num_data_sets,max(Num_Cell_Per_Exp_vec));
    %cl_var.Full_RF_Ellipticity_hist_BinEdges_mat = cell(p.Num_data_sets);
    %cl_var.Full_RF_Ellipticity_hist_NumBins_vec  = NaN(p.Num_data_sets,1);
    cl_var.Full_RF_Ellipticity_hist_BinWidth_vec  = NaN(p.Num_data_sets,1);
    cl_var.Full_RF_Ellipticity_hist_BinLimits_mat = NaN(p.Num_data_sets,2);
    
    RF_Ellipticity_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    RF_Ellipticity_Plot_Row = ceil(p.Num_data_sets/RF_Ellipticity_Plot_Dim);
    RF_Ellipticity_Plot_Col = RF_Ellipticity_Plot_Dim;
    
    RF_Ellipticity_Fig = figure;
    
    for i = 1:p.Num_data_sets
        
        % Load data
        if RF_Ident_Method_RF_Ellipticity == 1     % STA-SD
            %load('data_CNoise_20px_Paul_zf_RF_test_20_AllStat.mat','STASD_Gaus_FullRF_Axis_eval_Ratio_vec');
            load(CNoise_20px_DataFiles{i},'STASD_Gaus_FullRF_Axis_eval_Ratio_vec','True_cell_index_vec');
            Full_RF_Ellipticity_mat(i,True_cell_index_vec) = STASD_Gaus_FullRF_Axis_eval_Ratio_vec; %Full_RF_Ellipticity_mat(i,1:Num_Cell_Per_Exp_vec(i))
            clear STASD_Gaus_FullRF_Axis_eval_Ratio_vec True_cell_index_vec;
        else % RF_Ident_Method_RF_Ellipticity == 2 % SC
            %load('data_CNoise_20px_Paul_zf_RF_test_20_AllStat.mat','SC_Gaus_FullRF_Axis_eval_Ratio_vec');
            load(CNoise_20px_DataFiles{i},'SC_Gaus_FullRF_Axis_eval_Ratio_vec','True_cell_index_vec');
            Full_RF_Ellipticity_mat(i,True_cell_index_vec) = SC_Gaus_FullRF_Axis_eval_Ratio_vec; %Full_RF_Ellipticity_mat(i,1:Num_Cell_Per_Exp_vec(i))
            clear SC_Gaus_FullRF_Axis_eval_Ratio_vec True_cell_index_vec;
        end
        %STASD_Gaus_FullRF_Axis_eval_Ratio_vec(STASD_Gaus_FullRF_Axis_eval_Ratio_vec<=20) % I have removed RF ellipticities > 20
        %SC_Gaus_FullRF_Axis_eval_Ratio_vec(SC_Gaus_FullRF_Axis_eval_Ratio_vec<=20)       % I have removed RF ellipticities > 20
        
        if p.RF_Ellip_thresh_choice == 1
            Full_RF_Ellipticity_mat(i,Full_RF_Ellipticity_mat(i,1:Num_Cell_Per_Exp_vec(i))>p.RF_Ellip_thresh) = p.RF_Ellip_thresh;
            %cl_var.Full_RF_Ellipticity_vec(cl_var.Full_RF_Ellipticity_vec>p.RF_Ellip_thresh) = p.RF_Ellip_thresh;
        end
        
        % Plot data
        figure(RF_Ellipticity_Fig);
        subplot(RF_Ellipticity_Plot_Row,RF_Ellipticity_Plot_Col,i);
        RF_Ellip_hist = histogram(Full_RF_Ellipticity_mat(i,1:Num_Cell_Per_Exp_vec(i)));
        xlim([RF_Ellip_hist.BinEdges(1) RF_Ellip_hist.BinEdges(end)]);
        ylim([0 max(RF_Ellip_hist.Values)]);
        if (i/RF_Ellipticity_Plot_Col > RF_Ellipticity_Plot_Row-1) || ((i/RF_Ellipticity_Plot_Col > RF_Ellipticity_Plot_Row-2)&&(i/RF_Ellipticity_Plot_Col < RF_Ellipticity_Plot_Row-1)&&(mod(p.Num_data_sets,RF_Ellipticity_Plot_Col)~=0)&&(mod(i,RF_Ellipticity_Plot_Col)>mod(p.Num_data_sets,RF_Ellipticity_Plot_Col)))
            xlabel('RF ellipticity');
        end
        if mod(i,RF_Ellipticity_Plot_Col) == 1
            ylabel('num. cells');
        end
        title(Gen_Data_Names_vec{i}); % RF ellipticity
        set(gca,'FontSize',12);
        
        %cl_var.Full_RF_Ellipticity_hist_BinEdges_mat{i}   = RF_Ellip_hist.BinEdges;
        %cl_var.Full_RF_Ellipticity_hist_NumBins_vec(i)    = RF_Ellip_hist.NumBins;
        cl_var.Full_RF_Ellipticity_hist_BinWidth_vec(i)    = RF_Ellip_hist.BinWidth;
        cl_var.Full_RF_Ellipticity_hist_BinLimits_mat(i,:) = RF_Ellip_hist.BinLimits;
        %cl_var.RF_Ellip_hist_MaxValue = max(RF_Ellip_hist.Values); % Not required as calc max for each cluster
        
    end
    
    figure(RF_Ellipticity_Fig);
    annotation('textbox',[.5 .975 0 0],'String','RF ellipticity','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    set(RF_Ellipticity_Fig,'color','w');
    
    cl_var.Full_RF_Ellipticity_vec = NaN(Total_Num_Cells,1);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.Full_RF_Ellipticity_vec(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = Full_RF_Ellipticity_mat(i,1:Num_Cell_Per_Exp_vec(i));
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
%     if p.RF_Ellip_thresh_choice == 1 % Do this above to avoid wrong
%     histogram ppties.
%         cl_var.Full_RF_Ellipticity_vec(cl_var.Full_RF_Ellipticity_vec>p.RF_Ellip_thresh) = p.RF_Ellip_thresh;
%     end
    
end


%%% CNoise: Full RF Dominant Axis Angle
if p.Obj_clust_vec(7) == 1 || p.Obj_plot_vec(9) == 1
    
    Full_RF_Dom_Ax_Ang_mat                       = NaN(p.Num_data_sets,max(Num_Cell_Per_Exp_vec));
    cl_var.Full_RF_Dom_Ax_Ang_hist_BinWidth_vec  = NaN(p.Num_data_sets,1);
    cl_var.Full_RF_Dom_Ax_Ang_hist_BinLimits_mat = NaN(p.Num_data_sets,2);
    
    RF_Dom_Ax_Ang_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    RF_Dom_Ax_Ang_Plot_Row = ceil(p.Num_data_sets/RF_Dom_Ax_Ang_Plot_Dim);
    RF_Dom_Ax_Ang_Plot_Col = RF_Dom_Ax_Ang_Plot_Dim;
    
    RF_Dom_Ax_Ang_Fig = figure;
    
    for i = 1:p.Num_data_sets
        
        % Load data
        if RF_Ident_Method_RF_Dom_Ax_Ang == 1     % STA-SD
            load(CNoise_20px_DataFiles{i},'STASD_Gaus_FullRF_Angle_Major_Axis_vec','STASD_Gaus_FullRF_Axis_eval_Ratio_vec');
            Full_RF_Dom_Ax_Ang_mat(i,STASD_Gaus_FullRF_Axis_eval_Ratio_vec~=1) = STASD_Gaus_FullRF_Angle_Major_Axis_vec(STASD_Gaus_FullRF_Axis_eval_Ratio_vec~=1); % Full_RF_Dom_Ax_Ang_mat(i,1:Num_Cell_Per_Exp_vec(i))
            clear STASD_Gaus_FullRF_Angle_Major_Axis_vec STASD_Gaus_FullRF_Axis_eval_Ratio_vec;
        else % RF_Ident_Method_RF_Dom_Ax_Ang == 2 % SC
            load(CNoise_20px_DataFiles{i},'SC_Gaus_RF_Angle_Major_Axis_mat','SC_Gaus_FullRF_Axis_eval_Ratio_vec');
            Full_RF_Dom_Ax_Ang_mat(i,SC_Gaus_FullRF_Axis_eval_Ratio_vec~=1) = SC_Gaus_RF_Angle_Major_Axis_mat(SC_Gaus_FullRF_Axis_eval_Ratio_vec~=1); % Full_RF_Dom_Ax_Ang_mat(i,1:Num_Cell_Per_Exp_vec(i))
            clear SC_Gaus_RF_Angle_Major_Axis_mat SC_Gaus_FullRF_Axis_eval_Ratio_vec;
        end
        
        % Plot data
        figure(RF_Dom_Ax_Ang_Fig);
        subplot(RF_Dom_Ax_Ang_Plot_Row,RF_Dom_Ax_Ang_Plot_Col,i);
        %RF_Dom_Ax_Ang_hist = histogram(Full_RF_Dom_Ax_Ang_mat(i,1:Num_Cell_Per_Exp_vec(i)));
        RF_Dom_Ax_Ang_hist = polarhistogram(Full_RF_Dom_Ax_Ang_mat(i,1:Num_Cell_Per_Exp_vec(i)),(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]);
        %xlim([RF_Dom_Ax_Ang_hist.BinEdges(1) RF_Dom_Ax_Ang_hist.BinEdges(end)]);
        %ylim([0 max(RF_Dom_Ax_Ang_hist.Values)]);
        %         if (i/RF_Dom_Ax_Ang_Plot_Col > RF_Dom_Ax_Ang_Plot_Row-1) || ((i/RF_Dom_Ax_Ang_Plot_Col > RF_Dom_Ax_Ang_Plot_Row-2)&&(i/RF_Dom_Ax_Ang_Plot_Col < RF_Dom_Ax_Ang_Plot_Row-1)&&(mod(p.Num_data_sets,RF_Dom_Ax_Ang_Plot_Col)~=0)&&(mod(i,RF_Dom_Ax_Ang_Plot_Col)>mod(p.Num_data_sets,RF_Dom_Ax_Ang_Plot_Col)))
        %             xlabel('dom. axis angle');
        %         end
        %         if mod(i,RF_Dom_Ax_Ang_Plot_Col) == 1
        %             ylabel('num. cells');
        %         end
        set(gca,'ThetaTick',[0 45 90 270 315]);
        set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
        title(Gen_Data_Names_vec{i}); % RF size
        set(gca,'FontSize',12);
        
        %cl_var.Full_RF_Dom_Ax_Ang_hist_BinWidth_vec(i)    = RF_Dom_Ax_Ang_hist.BinWidth; % --> not required.
        %cl_var.Full_RF_Dom_Ax_Ang_hist_BinLimits_mat(i,:) = RF_Dom_Ax_Ang_hist.BinLimits; % --> not required.
        
    end
    
    figure(RF_Dom_Ax_Ang_Fig);
    annotation('textbox',[.5 .975 0 0],'String','RF dom. axis angle','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    set(RF_Dom_Ax_Ang_Fig,'color','w');
    
    cl_var.Full_RF_Dom_Ax_Ang_vec = NaN(Total_Num_Cells,1);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.Full_RF_Dom_Ax_Ang_vec(loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = Full_RF_Dom_Ax_Ang_mat(i,1:Num_Cell_Per_Exp_vec(i));
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
end



%%% FFF2
if p.Obj_clust_vec(8) == 1 || p.Obj_plot_vec(10) == 1
    
    FFF2_NumTrigPerRep     = 6; % User defined
    FFF2_spike_times_arr   = cell(p.Num_data_sets,1);
    FFF2_trig_times_mat    = cell(p.Num_data_sets,1);
    FFF2_stim_end_time_vec = NaN(p.Num_data_sets,1);
    FFF2_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    
    FFF2_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    FFF2_Plot_Row = ceil(p.Num_data_sets/FFF2_Plot_Dim);
    FFF2_Plot_Col = FFF2_Plot_Dim;
    
    FFF2_QI_vec                     = NaN(Total_Num_Cells,1);
    FFF2_Qual_length_ksdensity_grid = 1000; % 100 coule use 1e3 (size(FFF2_spike_times_arr{i},1))
    FFF2_Qual_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best
    
    FFF2_Trig_Fig      = figure;
    FFF2_Trig_Dist_Fig = figure;
    FFF2_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(FFF2_DataFiles{i});
        
        % Spike times
        if Cell_Choice == 1 % All cells
            FFF2_spike_times_arr{i} = spiketimestamps;
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define FFF2 version)
            FFF2_spike_times_arr{i} = spiketimestamps(:,Cell_Choice_vec);
        end
        
        % Triggers
        FFF2_trigCh_vec            = Ch_new.trigger_ch;
        if isa(FFF2_trigCh_vec,'cell')
            FFF2_trigCh_vec = cell2mat(FFF2_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        FFF2_min_trigCh_vec        = min(FFF2_trigCh_vec);
        FFF2_max_trigCh_vec        = max(FFF2_trigCh_vec);
        FFF2_trigThreshFac         = 0.05;
        FFF2_trigHigh_vec          = double(FFF2_trigCh_vec > FFF2_min_trigCh_vec + FFF2_trigThreshFac*(FFF2_max_trigCh_vec-FFF2_min_trigCh_vec));
        [~,FFF2_trig_index_vec]    = findpeaks(FFF2_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        FFF2_trig_index_vec_length = length(FFF2_trig_index_vec);
        
        FFF2_NumTrigRep_mod = mod(FFF2_trig_index_vec_length,FFF2_NumTrigPerRep); % New: 28,05,2021
        
        if FFF2_NumTrigRep_mod ~= 0 % New: 28,05,2021 (whole if statment (last line was already there w/o if statement))
            FFF2_trig_index_vec = FFF2_trig_index_vec(1:end-FFF2_NumTrigRep_mod);
            FFF2_trig_index_vec_length = length(FFF2_trig_index_vec);    
        end
        FFF2_NumTrigRep = FFF2_trig_index_vec_length/FFF2_NumTrigPerRep;
        
        FFF2_sampling_freq = Ch_new.SamplingFrequency;
        FFF2_sampling_int  = 1/FFF2_sampling_freq;
        FFF2_trig_t_vec    = (0:1:length(FFF2_trigCh_vec)-1)*FFF2_sampling_int;
        
        FFF2_trig_times_vec_temp = FFF2_trig_t_vec(FFF2_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        FFF2_trig_times_mat{i}   = FFF2_trig_times_vec_temp - FFF2_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(FFF2_Trig_Fig);
        subplot(FFF2_Plot_Row,FFF2_Plot_Col,i);
        plot(FFF2_trig_t_vec,FFF2_trigCh_vec); hold on;
        plot(FFF2_trig_times_vec_temp,4085*ones(1,FFF2_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % FFF2 trigger channel
        if (i/FFF2_Plot_Col > FFF2_Plot_Row-1) || ((i/FFF2_Plot_Col > FFF2_Plot_Row-2)&&(i/FFF2_Plot_Col < FFF2_Plot_Row-1)&&(mod(p.Num_data_sets,FFF2_Plot_Col)~=0)&&(mod(i,FFF2_Plot_Col)>mod(p.Num_data_sets,FFF2_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,FFF2_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        % Find and plot the distribution of trigger intervals
        FFF2_actual_stim_int      = diff(FFF2_trig_times_mat{i});
        FFF2_mean_stim_int_vec(i) = mean(FFF2_actual_stim_int);
        FFF2_stddev_stim_int      = std(FFF2_actual_stim_int);
        figure(FFF2_Trig_Dist_Fig);
        subplot(FFF2_Plot_Row,FFF2_Plot_Col,i);
        histogram(FFF2_actual_stim_int); hold on;
        FFF2_v1 = vline(FFF2_mean_stim_int_vec(i),'r');
        FFF2_v2 = vline(FFF2_mean_stim_int_vec(i)-FFF2_stddev_stim_int,'g');
        FFF2_v3 = vline(FFF2_mean_stim_int_vec(i)+FFF2_stddev_stim_int,'g');
        if (i/FFF2_Plot_Col > FFF2_Plot_Row-1) || ((i/FFF2_Plot_Col > FFF2_Plot_Row-2)&&(i/FFF2_Plot_Col < FFF2_Plot_Row-1)&&(mod(p.Num_data_sets,FFF2_Plot_Col)~=0)&&(mod(i,FFF2_Plot_Col)>mod(p.Num_data_sets,FFF2_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,FFF2_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % FFF2 trigger int. dist.
        set(FFF2_v1,'LineWidth',1.5);
        set(FFF2_v2,'LineWidth',1.5);
        set(FFF2_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        
        % I find 40 triggers - so 10 repeats.
        % Trig 2s red 2s dark, Trig 2s green 2s dark, Trig 2s blue 2s dark, Trig 2s UV 2s dark
        % Hence 4 sec between triggers, 16 sec for 1 repeat and 160 sec for 10 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cell's spikes into increasing time order.
        % Check:
        % Paul = [1 2;3 4];
        % Paul(Paul<3) = NaN
        % Paul = [Paul; 1 2]
        % Paul = sort(Paul,1)
        
        % 1. remove spikes before and after the stimulus;
        % Test Before
        % min(FFF2_spike_times_arr{i},[],'all')
        % max(FFF2_spike_times_arr{i},[],'all')
        FFF2_spike_times_arr{i}(FFF2_spike_times_arr{i}<=FFF2_trig_times_mat{i}(1)) = NaN;
        FFF2_spike_times_arr{i}(FFF2_spike_times_arr{i}>=(FFF2_trig_times_mat{i}(end) + FFF2_mean_stim_int_vec(i))) = NaN;
        % Test After
        % min(FFF2_spike_times_arr{i},[],'all')
        % max(FFF2_spike_times_arr{i},[],'all')
        
        
        if FFF2_trig_index_vec_length>FFF2_NumTrigPerRep
            FFF2_stim_end_time_vec(i) = FFF2_trig_times_mat{i}(FFF2_NumTrigPerRep+1); % period of a single repeat (Time of start of 2nd repeat (hence length of first repeat))
        else % FFF2_trig_index_vec_length<=FFF2_NumTrigPerRep
            FFF2_stim_end_time_vec(i) = 24; % sec (6 trigs, 4 sec between each)
        end
        
        % 1.5 Calculate Quality Indices 
        FFF2_Qual_ksdensity_grid  = linspace(0,FFF2_stim_end_time_vec(i),FFF2_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(FFF2_Qual_length_ksdensity_grid,FFF2_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = FFF2_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/FFF2_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:FFF2_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),FFF2_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*FFF2_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,FFF2_Qual_ksdensity_grid,'Bandwidth',FFF2_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                FFF2_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                FFF2_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % FFF2_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(FFF2_QI_vec);
        %figure;
        %hist(FFF2_QI_vec);
        
        % Plot the distribution of QIs
        figure(FFF2_QI_Dist_Fig);
        subplot(FFF2_Plot_Row,FFF2_Plot_Col,i);
        if i == 1
            Plot_vec_loop = FFF2_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = FFF2_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        FFF2_QI_v1 = vline(FFF2_QI_thresh,'r','QI thresh.');
        FFF2_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/FFF2_Plot_Col > FFF2_Plot_Row-1) || ((i/FFF2_Plot_Col > FFF2_Plot_Row-2)&&(i/FFF2_Plot_Col < FFF2_Plot_Row-1)&&(mod(p.Num_data_sets,FFF2_Plot_Col)~=0)&&(mod(i,FFF2_Plot_Col)>mod(p.Num_data_sets,FFF2_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,FFF2_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(FFF2_QI_v1,'LineWidth',1.5);
        set(FFF2_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if FFF2_trig_index_vec_length>FFF2_NumTrigPerRep
           % FFF2_stim_end_time_vec(i) =
           % FFF2_trig_times_mat{i}(FFF2_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            FFF2_spike_times_arr{i}   = mod(FFF2_spike_times_arr{i},FFF2_stim_end_time_vec(i));
       % else % FFF2_trig_index_vec_length<=FFF2_NumTrigPerRep % --> commented as now use above
            %FFF2_stim_end_time_vec(i) = 16; % sec % --> commented as now use above
        end
        % FFF2_NumTrigPerRep
        % FFF2_NumTrigRep
        
        % 3. sort each cell's spikes into increasing time order.
        FFF2_spike_times_arr{i} = sort(FFF2_spike_times_arr{i},1);
        
    end
    
    figure(FFF2_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF2 trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF2_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF2 trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(FFF2_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','FFF2 QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(FFF2_Trig_Fig,'color','w');
    set(FFF2_Trig_Dist_Fig,'color','w');
    set(FFF2_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(FFF2_QI_vec); hold on;
    FFF2_QI_v1 = vline(FFF2_QI_thresh,'r','QI thresh.');
    FFF2_QI_v2 = vline(max(FFF2_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('FFF2 QI - All Cells');
    set(FFF2_QI_v1,'LineWidth',1.5);
    set(FFF2_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    
    % Calculate FFF2 QI logical vec
    FFF2_QI_logical_vec = FFF2_QI_vec > FFF2_QI_thresh;
    Num_FFF2_Qual_Cells = sum(FFF2_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_FFF2_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,FFF2_QI_logical_vec],[],2);
    
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1
    for j = 2:p.Num_data_sets
        for i = 1:FFF2_NumTrigPerRep
            if i < FFF2_NumTrigPerRep
                indices_loop = find((FFF2_spike_times_arr{j}>=FFF2_trig_times_mat{j}(i))&(FFF2_spike_times_arr{j}<FFF2_trig_times_mat{j}(i+1)));
                b_loop       = FFF2_trig_times_mat{1}(i+1);
                b_dash_loop  = FFF2_trig_times_mat{j}(i+1);
            else %i == FFF2_NumTrigPerRep
                if length(FFF2_trig_times_mat{j}) > FFF2_NumTrigPerRep
                    indices_loop = find((FFF2_spike_times_arr{j}>=FFF2_trig_times_mat{j}(i))&(FFF2_spike_times_arr{j}<=FFF2_trig_times_mat{j}(i+1)));
                    b_dash_loop  = FFF2_trig_times_mat{j}(i+1);
                else
                    indices_loop = find((FFF2_spike_times_arr{j}>=FFF2_trig_times_mat{j}(i))&(FFF2_spike_times_arr{j}<=FFF2_trig_times_mat{j}(i) + FFF2_mean_stim_int_vec(j)));
                    b_dash_loop  = FFF2_trig_times_mat{j}(i) + FFF2_mean_stim_int_vec(j);
                end
                if length(FFF2_trig_times_mat{1}) > FFF2_NumTrigPerRep % was FFF2_trig_times_mat{j} --> changed 08,04,2021
                    b_loop       = FFF2_trig_times_mat{1}(i+1);
                else
                    b_loop       = FFF2_trig_times_mat{1}(i) + FFF2_mean_stim_int_vec(1);
                end
            end
            a_loop       = FFF2_trig_times_mat{1}(i);
            a_dash_loop  = FFF2_trig_times_mat{j}(i);
            FFF2_spike_times_arr{j}(indices_loop) = a_loop + (FFF2_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    FFF2_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        FFF2_num_spikes_vec(i) = size(FFF2_spike_times_arr{i},1);
    end
    FFF2_max_num_spikes = max(FFF2_num_spikes_vec);
    
    cl_var.FFF2_spike_times_mat = NaN(FFF2_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.FFF2_spike_times_mat(1:FFF2_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = FFF2_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
%    cl_var.FFF2_spike_times_mat = ;                         % Just above.
%    cl_var.FFF2_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
     cl_var.FFF2_trig_times_vec  = FFF2_trig_times_mat{1};    % As map subsequent data set times back to first data set times
     cl_var.FFF2_stim_end_time   = FFF2_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
    
end



%%% SSub
if p.Obj_clust_vec(10) == 1 || p.Obj_plot_vec(12) == 1
    
    p.SSub_NumTrigPerRep   = 18; % User defined
    SSub_spike_times_arr   = cell(p.Num_data_sets,1);
    SSub_trig_times_mat    = cell(p.Num_data_sets,1);
    SSub_stim_end_time_vec = NaN(p.Num_data_sets,1);
    SSub_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    
    SSub_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    SSub_Plot_Row = ceil(p.Num_data_sets/SSub_Plot_Dim);
    SSub_Plot_Col = SSub_Plot_Dim;
    
    SSub_QI_vec                     = NaN(Total_Num_Cells,1);
    SSub_Qual_length_ksdensity_grid = 1000; % 100 could use 1e3 (size(SSub_spike_times_arr{i},1))
    SSub_Qual_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best for FFF
    
    SSub_Trig_Fig      = figure;
    SSub_Trig_Dist_Fig = figure;
    SSub_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(SSub_DataFiles{i});
        
        % Spike times
        if Cell_Choice == 1 % All cells
            SSub_spike_times_arr{i} = spiketimestamps;
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define SSub version)
            SSub_spike_times_arr{i} = spiketimestamps(:,Cell_Choice_vec);
        end
        
        % Triggers
        SSub_trigCh_vec            = Ch_new.trigger_ch;
        if isa(SSub_trigCh_vec,'cell')
            SSub_trigCh_vec = cell2mat(SSub_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        SSub_min_trigCh_vec        = min(SSub_trigCh_vec);
        SSub_max_trigCh_vec        = max(SSub_trigCh_vec);
        SSub_trigThreshFac         = 0.05;
        SSub_trigHigh_vec          = double(SSub_trigCh_vec > SSub_min_trigCh_vec + SSub_trigThreshFac*(SSub_max_trigCh_vec-SSub_min_trigCh_vec));
        [~,SSub_trig_index_vec]    = findpeaks(SSub_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        SSub_trig_index_vec_length = length(SSub_trig_index_vec);
        
        SSub_NumTrigRep_mod = mod(SSub_trig_index_vec_length,p.SSub_NumTrigPerRep); % New: 28,05,2021
        
        if SSub_NumTrigRep_mod ~= 0 % New: 28,05,2021 (whole if statment (last line was already there w/o if statement))
            SSub_trig_index_vec = SSub_trig_index_vec(1:end-SSub_NumTrigRep_mod);
            SSub_trig_index_vec_length = length(SSub_trig_index_vec);
        end
        SSub_NumTrigRep = SSub_trig_index_vec_length/p.SSub_NumTrigPerRep;
        
        SSub_sampling_freq = Ch_new.SamplingFrequency;
        SSub_sampling_int  = 1/SSub_sampling_freq;
        SSub_trig_t_vec    = (0:1:length(SSub_trigCh_vec)-1)*SSub_sampling_int;
        
        SSub_trig_times_vec_temp = SSub_trig_t_vec(SSub_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        SSub_trig_times_mat{i}   = SSub_trig_times_vec_temp - SSub_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(SSub_Trig_Fig);
        subplot(SSub_Plot_Row,SSub_Plot_Col,i);
        plot(SSub_trig_t_vec,SSub_trigCh_vec); hold on;
        plot(SSub_trig_times_vec_temp,4085*ones(1,SSub_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % SSub trigger channel
        if (i/SSub_Plot_Col > SSub_Plot_Row-1) || ((i/SSub_Plot_Col > SSub_Plot_Row-2)&&(i/SSub_Plot_Col < SSub_Plot_Row-1)&&(mod(p.Num_data_sets,SSub_Plot_Col)~=0)&&(mod(i,SSub_Plot_Col)>mod(p.Num_data_sets,SSub_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,SSub_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        % Find and plot the distribution of trigger intervals
        SSub_actual_stim_int      = diff(SSub_trig_times_mat{i});
        SSub_mean_stim_int_vec(i) = mean(SSub_actual_stim_int);
        SSub_stddev_stim_int      = std(SSub_actual_stim_int);
        figure(SSub_Trig_Dist_Fig);
        subplot(SSub_Plot_Row,SSub_Plot_Col,i);
        histogram(SSub_actual_stim_int); hold on;
        SSub_v1 = vline(SSub_mean_stim_int_vec(i),'r');
        SSub_v2 = vline(SSub_mean_stim_int_vec(i)-SSub_stddev_stim_int,'g');
        SSub_v3 = vline(SSub_mean_stim_int_vec(i)+SSub_stddev_stim_int,'g');
        if (i/SSub_Plot_Col > SSub_Plot_Row-1) || ((i/SSub_Plot_Col > SSub_Plot_Row-2)&&(i/SSub_Plot_Col < SSub_Plot_Row-1)&&(mod(p.Num_data_sets,SSub_Plot_Col)~=0)&&(mod(i,SSub_Plot_Col)>mod(p.Num_data_sets,SSub_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,SSub_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % SSub trigger int. dist.
        set(SSub_v1,'LineWidth',1.5);
        set(SSub_v2,'LineWidth',1.5);
        set(SSub_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        
        % I find 40 triggers - so 10 repeats.
        % Trig 2s red 2s dark, Trig 2s green 2s dark, Trig 2s blue 2s dark, Trig 2s UV 2s dark
        % Hence 4 sec between triggers, 16 sec for 1 repeat and 160 sec for 10 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cell's spikes into increasing time order.
        % Check:
        % Paul = [1 2;3 4];
        % Paul(Paul<3) = NaN
        % Paul = [Paul; 1 2]
        % Paul = sort(Paul,1)
        
        % 1. remove spikes before and after the stimulus;
        % Test Before
        % min(SSub_spike_times_arr{i},[],'all')
        % max(SSub_spike_times_arr{i},[],'all')
        SSub_spike_times_arr{i}(SSub_spike_times_arr{i}<=SSub_trig_times_mat{i}(1)) = NaN;
        SSub_spike_times_arr{i}(SSub_spike_times_arr{i}>=(SSub_trig_times_mat{i}(end) + SSub_mean_stim_int_vec(i))) = NaN;
        % Test After
        % min(SSub_spike_times_arr{i},[],'all')
        % max(SSub_spike_times_arr{i},[],'all')
        
        
        if SSub_trig_index_vec_length>p.SSub_NumTrigPerRep
            SSub_stim_end_time_vec(i) = SSub_trig_times_mat{i}(p.SSub_NumTrigPerRep+1); % period of a single repeat (Time of start of 2nd repeat (hence length of first repeat))
        else % SSub_trig_index_vec_length<=p.SSub_NumTrigPerRep
            SSub_stim_end_time_vec(i) = 36; % sec (18 trigs, 2 sec between each)
        end
        
        % 1.5 Calculate Quality Indices 
        SSub_Qual_ksdensity_grid  = linspace(0,SSub_stim_end_time_vec(i),SSub_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(SSub_Qual_length_ksdensity_grid,SSub_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = SSub_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/SSub_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:SSub_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),SSub_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*SSub_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,SSub_Qual_ksdensity_grid,'Bandwidth',SSub_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                SSub_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                SSub_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % SSub_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(SSub_QI_vec);
        %figure;
        %hist(SSub_QI_vec);
        
        % Plot the distribution of QIs
        figure(SSub_QI_Dist_Fig);
        subplot(SSub_Plot_Row,SSub_Plot_Col,i);
        if i == 1
            Plot_vec_loop = SSub_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = SSub_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        SSub_QI_v1 = vline(SSub_QI_thresh,'r','QI thresh.');
        SSub_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/SSub_Plot_Col > SSub_Plot_Row-1) || ((i/SSub_Plot_Col > SSub_Plot_Row-2)&&(i/SSub_Plot_Col < SSub_Plot_Row-1)&&(mod(p.Num_data_sets,SSub_Plot_Col)~=0)&&(mod(i,SSub_Plot_Col)>mod(p.Num_data_sets,SSub_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,SSub_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(SSub_QI_v1,'LineWidth',1.5);
        set(SSub_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if SSub_trig_index_vec_length>p.SSub_NumTrigPerRep
           % SSub_stim_end_time_vec(i) =
           % SSub_trig_times_mat{i}(p.SSub_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            SSub_spike_times_arr{i}   = mod(SSub_spike_times_arr{i},SSub_stim_end_time_vec(i));
       % else % SSub_trig_index_vec_length<=p.SSub_NumTrigPerRep % --> commented as now use above
            %SSub_stim_end_time_vec(i) = 16; % sec % --> commented as now use above
        end
        % p.SSub_NumTrigPerRep
        % SSub_NumTrigRep
        
        % 3. sort each cell's spikes into increasing time order.
        SSub_spike_times_arr{i} = sort(SSub_spike_times_arr{i},1);
        
    end
    
    figure(SSub_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','SSub trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(SSub_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','SSub trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(SSub_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','SSub QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(SSub_Trig_Fig,'color','w');
    set(SSub_Trig_Dist_Fig,'color','w');
    set(SSub_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(SSub_QI_vec); hold on;
    SSub_QI_v1 = vline(SSub_QI_thresh,'r','QI thresh.');
    SSub_QI_v2 = vline(max(SSub_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('SSub QI - All Cells');
    set(SSub_QI_v1,'LineWidth',1.5);
    set(SSub_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    
    % Calculate SSub QI logical vec
    SSub_QI_logical_vec = SSub_QI_vec > SSub_QI_thresh;
    Num_SSub_Qual_Cells = sum(SSub_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_SSub_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,SSub_QI_logical_vec],[],2);
    
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1
    for j = 2:p.Num_data_sets
        for i = 1:p.SSub_NumTrigPerRep
            if i < p.SSub_NumTrigPerRep
                indices_loop = find((SSub_spike_times_arr{j}>=SSub_trig_times_mat{j}(i))&(SSub_spike_times_arr{j}<SSub_trig_times_mat{j}(i+1)));
                b_loop       = SSub_trig_times_mat{1}(i+1);
                b_dash_loop  = SSub_trig_times_mat{j}(i+1);
            else %i == p.SSub_NumTrigPerRep
                if length(SSub_trig_times_mat{j}) > p.SSub_NumTrigPerRep
                    indices_loop = find((SSub_spike_times_arr{j}>=SSub_trig_times_mat{j}(i))&(SSub_spike_times_arr{j}<=SSub_trig_times_mat{j}(i+1)));
                    b_dash_loop  = SSub_trig_times_mat{j}(i+1);
                else
                    indices_loop = find((SSub_spike_times_arr{j}>=SSub_trig_times_mat{j}(i))&(SSub_spike_times_arr{j}<=SSub_trig_times_mat{j}(i) + SSub_mean_stim_int_vec(j)));
                    b_dash_loop  = SSub_trig_times_mat{j}(i) + SSub_mean_stim_int_vec(j);
                end
                if length(SSub_trig_times_mat{1}) > p.SSub_NumTrigPerRep % was SSub_trig_times_mat{j} --> changed 08,04,2021
                    b_loop       = SSub_trig_times_mat{1}(i+1);
                else
                    b_loop       = SSub_trig_times_mat{1}(i) + SSub_mean_stim_int_vec(1);
                end
            end
            a_loop       = SSub_trig_times_mat{1}(i);
            a_dash_loop  = SSub_trig_times_mat{j}(i);
            SSub_spike_times_arr{j}(indices_loop) = a_loop + (SSub_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    SSub_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        SSub_num_spikes_vec(i) = size(SSub_spike_times_arr{i},1);
    end
    SSub_max_num_spikes = max(SSub_num_spikes_vec);
    
    cl_var.SSub_spike_times_mat = NaN(SSub_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.SSub_spike_times_mat(1:SSub_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = SSub_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
%    cl_var.SSub_spike_times_mat = ;                         % Just above.
%    cl_var.SSub_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
     cl_var.SSub_trig_times_vec  = SSub_trig_times_mat{1};    % As map subsequent data set times back to first data set times
     cl_var.SSub_stim_end_time   = SSub_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
    
end


%%% CSteps
if p.Obj_clust_vec(11) == 1 || p.Obj_plot_vec(13) == 1
    
    p.CSteps_NumTrigPerRep   = 10; % User defined
    CSteps_spike_times_arr   = cell(p.Num_data_sets,1);
    CSteps_trig_times_mat    = cell(p.Num_data_sets,1);
    CSteps_stim_end_time_vec = NaN(p.Num_data_sets,1);
    CSteps_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    
    CSteps_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    CSteps_Plot_Row = ceil(p.Num_data_sets/CSteps_Plot_Dim);
    CSteps_Plot_Col = CSteps_Plot_Dim;
    
    CSteps_QI_vec                     = NaN(Total_Num_Cells,1);
    CSteps_Qual_length_ksdensity_grid = 1000; % 100 could use 1e3 (size(CSteps_spike_times_arr{i},1))
    CSteps_Qual_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best for FFF
    
    CSteps_Trig_Fig      = figure;
    CSteps_Trig_Dist_Fig = figure;
    CSteps_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(CSteps_DataFiles{i});
        
        % Spike times
        if Cell_Choice == 1 % All cells
            CSteps_spike_times_arr{i} = spiketimestamps;
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define CSteps version)
            CSteps_spike_times_arr{i} = spiketimestamps(:,Cell_Choice_vec);
        end
        
        % Triggers
        CSteps_trigCh_vec            = Ch_new.trigger_ch;
        if isa(CSteps_trigCh_vec,'cell')
            CSteps_trigCh_vec = cell2mat(CSteps_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        CSteps_min_trigCh_vec        = min(CSteps_trigCh_vec);
        CSteps_max_trigCh_vec        = max(CSteps_trigCh_vec);
        CSteps_trigThreshFac         = 0.05;
        CSteps_trigHigh_vec          = double(CSteps_trigCh_vec > CSteps_min_trigCh_vec + CSteps_trigThreshFac*(CSteps_max_trigCh_vec-CSteps_min_trigCh_vec));
        [~,CSteps_trig_index_vec]    = findpeaks(CSteps_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        CSteps_trig_index_vec_length = length(CSteps_trig_index_vec);
        
        CSteps_NumTrigRep_mod = mod(CSteps_trig_index_vec_length,p.CSteps_NumTrigPerRep); % New: 28,05,2021
        
        if CSteps_NumTrigRep_mod ~= 0 % New: 28,05,2021 (whole if statment (last line was already there w/o if statement))
            CSteps_trig_index_vec = CSteps_trig_index_vec(1:end-CSteps_NumTrigRep_mod);
            CSteps_trig_index_vec_length = length(CSteps_trig_index_vec);
        end
        CSteps_NumTrigRep = CSteps_trig_index_vec_length/p.CSteps_NumTrigPerRep;
        
        CSteps_sampling_freq = Ch_new.SamplingFrequency;
        CSteps_sampling_int  = 1/CSteps_sampling_freq;
        CSteps_trig_t_vec    = (0:1:length(CSteps_trigCh_vec)-1)*CSteps_sampling_int;
        
        CSteps_trig_times_vec_temp = CSteps_trig_t_vec(CSteps_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        CSteps_trig_times_mat{i}   = CSteps_trig_times_vec_temp - CSteps_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(CSteps_Trig_Fig);
        subplot(CSteps_Plot_Row,CSteps_Plot_Col,i);
        plot(CSteps_trig_t_vec,CSteps_trigCh_vec); hold on;
        plot(CSteps_trig_times_vec_temp,4085*ones(1,CSteps_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % CSteps trigger channel
        if (i/CSteps_Plot_Col > CSteps_Plot_Row-1) || ((i/CSteps_Plot_Col > CSteps_Plot_Row-2)&&(i/CSteps_Plot_Col < CSteps_Plot_Row-1)&&(mod(p.Num_data_sets,CSteps_Plot_Col)~=0)&&(mod(i,CSteps_Plot_Col)>mod(p.Num_data_sets,CSteps_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,CSteps_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        % Find and plot the distribution of trigger intervals
        CSteps_actual_stim_int      = diff(CSteps_trig_times_mat{i});
        CSteps_mean_stim_int_vec(i) = mean(CSteps_actual_stim_int);
        CSteps_stddev_stim_int      = std(CSteps_actual_stim_int);
        figure(CSteps_Trig_Dist_Fig);
        subplot(CSteps_Plot_Row,CSteps_Plot_Col,i);
        histogram(CSteps_actual_stim_int); hold on;
        CSteps_v1 = vline(CSteps_mean_stim_int_vec(i),'r');
        CSteps_v2 = vline(CSteps_mean_stim_int_vec(i)-CSteps_stddev_stim_int,'g');
        CSteps_v3 = vline(CSteps_mean_stim_int_vec(i)+CSteps_stddev_stim_int,'g');
        if (i/CSteps_Plot_Col > CSteps_Plot_Row-1) || ((i/CSteps_Plot_Col > CSteps_Plot_Row-2)&&(i/CSteps_Plot_Col < CSteps_Plot_Row-1)&&(mod(p.Num_data_sets,CSteps_Plot_Col)~=0)&&(mod(i,CSteps_Plot_Col)>mod(p.Num_data_sets,CSteps_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,CSteps_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % CSteps trigger int. dist.
        set(CSteps_v1,'LineWidth',1.5);
        set(CSteps_v2,'LineWidth',1.5);
        set(CSteps_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        
        % I find 40 triggers - so 10 repeats.
        % Trig 2s red 2s dark, Trig 2s green 2s dark, Trig 2s blue 2s dark, Trig 2s UV 2s dark
        % Hence 4 sec between triggers, 16 sec for 1 repeat and 160 sec for 10 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cell's spikes into increasing time order.
        % Check:
        % Paul = [1 2;3 4];
        % Paul(Paul<3) = NaN
        % Paul = [Paul; 1 2]
        % Paul = sort(Paul,1)
        
        % 1. remove spikes before and after the stimulus;
        % Test Before
        % min(CSteps_spike_times_arr{i},[],'all')
        % max(CSteps_spike_times_arr{i},[],'all')
        CSteps_spike_times_arr{i}(CSteps_spike_times_arr{i}<=CSteps_trig_times_mat{i}(1)) = NaN;
        CSteps_spike_times_arr{i}(CSteps_spike_times_arr{i}>=(CSteps_trig_times_mat{i}(end) + CSteps_mean_stim_int_vec(i))) = NaN;
        % Test After
        % min(CSteps_spike_times_arr{i},[],'all')
        % max(CSteps_spike_times_arr{i},[],'all')
        
        
        if CSteps_trig_index_vec_length>p.CSteps_NumTrigPerRep
            CSteps_stim_end_time_vec(i) = CSteps_trig_times_mat{i}(p.CSteps_NumTrigPerRep+1); % period of a single repeat (Time of start of 2nd repeat (hence length of first repeat))
        else % CSteps_trig_index_vec_length<=p.CSteps_NumTrigPerRep
            CSteps_stim_end_time_vec(i) = 40; % sec (10 trigs, 4 sec between each)
        end
        
        % 1.5 Calculate Quality Indices 
        CSteps_Qual_ksdensity_grid  = linspace(0,CSteps_stim_end_time_vec(i),CSteps_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(CSteps_Qual_length_ksdensity_grid,CSteps_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = CSteps_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/CSteps_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:CSteps_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),CSteps_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*CSteps_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,CSteps_Qual_ksdensity_grid,'Bandwidth',CSteps_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                CSteps_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                CSteps_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % CSteps_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(CSteps_QI_vec);
        %figure;
        %hist(CSteps_QI_vec);
        
        % Plot the distribution of QIs
        figure(CSteps_QI_Dist_Fig);
        subplot(CSteps_Plot_Row,CSteps_Plot_Col,i);
        if i == 1
            Plot_vec_loop = CSteps_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = CSteps_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        CSteps_QI_v1 = vline(CSteps_QI_thresh,'r','QI thresh.');
        CSteps_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/CSteps_Plot_Col > CSteps_Plot_Row-1) || ((i/CSteps_Plot_Col > CSteps_Plot_Row-2)&&(i/CSteps_Plot_Col < CSteps_Plot_Row-1)&&(mod(p.Num_data_sets,CSteps_Plot_Col)~=0)&&(mod(i,CSteps_Plot_Col)>mod(p.Num_data_sets,CSteps_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,CSteps_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(CSteps_QI_v1,'LineWidth',1.5);
        set(CSteps_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if CSteps_trig_index_vec_length>p.CSteps_NumTrigPerRep
           % CSteps_stim_end_time_vec(i) =
           % CSteps_trig_times_mat{i}(p.CSteps_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            CSteps_spike_times_arr{i}   = mod(CSteps_spike_times_arr{i},CSteps_stim_end_time_vec(i));
       % else % CSteps_trig_index_vec_length<=p.CSteps_NumTrigPerRep % --> commented as now use above
            %CSteps_stim_end_time_vec(i) = 16; % sec % --> commented as now use above
        end
        % p.CSteps_NumTrigPerRep
        % CSteps_NumTrigRep
        
        % 3. sort each cell's spikes into increasing time order.
        CSteps_spike_times_arr{i} = sort(CSteps_spike_times_arr{i},1);
        
    end
    
    figure(CSteps_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','CSteps trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(CSteps_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','CSteps trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(CSteps_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','CSteps QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(CSteps_Trig_Fig,'color','w');
    set(CSteps_Trig_Dist_Fig,'color','w');
    set(CSteps_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(CSteps_QI_vec); hold on;
    CSteps_QI_v1 = vline(CSteps_QI_thresh,'r','QI thresh.');
    CSteps_QI_v2 = vline(max(CSteps_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('CSteps QI - All Cells');
    set(CSteps_QI_v1,'LineWidth',1.5);
    set(CSteps_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    
    % Calculate CSteps QI logical vec
    CSteps_QI_logical_vec = CSteps_QI_vec > CSteps_QI_thresh;
    Num_CSteps_Qual_Cells = sum(CSteps_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_CSteps_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,CSteps_QI_logical_vec],[],2);
    
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1
    for j = 2:p.Num_data_sets
        for i = 1:p.CSteps_NumTrigPerRep
            if i < p.CSteps_NumTrigPerRep
                indices_loop = find((CSteps_spike_times_arr{j}>=CSteps_trig_times_mat{j}(i))&(CSteps_spike_times_arr{j}<CSteps_trig_times_mat{j}(i+1)));
                b_loop       = CSteps_trig_times_mat{1}(i+1);
                b_dash_loop  = CSteps_trig_times_mat{j}(i+1);
            else %i == p.CSteps_NumTrigPerRep
                if length(CSteps_trig_times_mat{j}) > p.CSteps_NumTrigPerRep
                    indices_loop = find((CSteps_spike_times_arr{j}>=CSteps_trig_times_mat{j}(i))&(CSteps_spike_times_arr{j}<=CSteps_trig_times_mat{j}(i+1)));
                    b_dash_loop  = CSteps_trig_times_mat{j}(i+1);
                else
                    indices_loop = find((CSteps_spike_times_arr{j}>=CSteps_trig_times_mat{j}(i))&(CSteps_spike_times_arr{j}<=CSteps_trig_times_mat{j}(i) + CSteps_mean_stim_int_vec(j)));
                    b_dash_loop  = CSteps_trig_times_mat{j}(i) + CSteps_mean_stim_int_vec(j);
                end
                if length(CSteps_trig_times_mat{1}) > p.CSteps_NumTrigPerRep % was CSteps_trig_times_mat{j} --> changed 08,04,2021
                    b_loop       = CSteps_trig_times_mat{1}(i+1);
                else
                    b_loop       = CSteps_trig_times_mat{1}(i) + CSteps_mean_stim_int_vec(1);
                end
            end
            a_loop       = CSteps_trig_times_mat{1}(i);
            a_dash_loop  = CSteps_trig_times_mat{j}(i);
            CSteps_spike_times_arr{j}(indices_loop) = a_loop + (CSteps_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    CSteps_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        CSteps_num_spikes_vec(i) = size(CSteps_spike_times_arr{i},1);
    end
    CSteps_max_num_spikes = max(CSteps_num_spikes_vec);
    
    cl_var.CSteps_spike_times_mat = NaN(CSteps_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.CSteps_spike_times_mat(1:CSteps_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = CSteps_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
%    cl_var.CSteps_spike_times_mat = ;                         % Just above.
%    cl_var.CSteps_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
     cl_var.CSteps_trig_times_vec  = CSteps_trig_times_mat{1};    % As map subsequent data set times back to first data set times
     cl_var.CSteps_stim_end_time   = CSteps_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
    
end



%%% Chirp2
if p.Obj_clust_vec(9) == 1 || p.Obj_plot_vec(11) == 1
    
    Chirp2_NumTrigPerRep     = 3; % User defined
    Chirp2_spike_times_arr   = cell(p.Num_data_sets,1);
    Chirp2_trig_times_mat    = cell(p.Num_data_sets,1);
    Chirp2_stim_end_time_vec = NaN(p.Num_data_sets,1);
    %Chirp2_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    Chirp2_last_stim_int_vec = NaN(p.Num_data_sets,1);
    
    Chirp2_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    Chirp2_Plot_Row = ceil(p.Num_data_sets/Chirp2_Plot_Dim);
    Chirp2_Plot_Col = Chirp2_Plot_Dim;
    
    Chirp2_QI_vec                     = NaN(Total_Num_Cells,1);
    Chirp2_Qual_length_ksdensity_grid = 1000; % 100 coule use 1e3
    Chirp2_Qual_ksdensity_bdwth       = 5*1e-2; % Tom said between 1e-2 and 1e-1 is best for Orig Chirp
    
    Chirp2_Trig_Fig      = figure;
    Chirp2_Trig_Dist_Fig = figure;
    Chirp2_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(Chirp2_DataFiles{i});
        
        % Spike times
        if Cell_Choice == 1 % All cells
            Chirp2_spike_times_arr{i}   = spiketimestamps;
            %Chirp2_Cell_Choice_vec      = 1:1:size(spiketimestamps,2);
            %Chirp2_True_cell_index_vec = Chirp2_Cell_Choice_vec(any(~isnan(Chirp2_spike_times_arr{i}))); % Record indices of responsive cells
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define Chirp2 version)
            Chirp2_spike_times_arr{i}   = spiketimestamps(:,Cell_Choice_vec);
            %Chirp2_True_cell_index_vec = Cell_Choice_vec(any(~isnan(Chirp2_spike_times_arr{i}))); % Record indices of responsive cells
        end
        
        % Remove NaN columns (non-responsive cells)?
        %Chirp2_spike_times_arr{i}(:,~any(~isnan(Chirp2_spike_times_arr{i}))) = [];
        %cl_var.Chirp2_True_Num_Cells = size(Chirp2_spike_times_arr{i},2);
        
        % Triggers
        Chirp2_trigCh_vec            = Ch_new.trigger_ch;
        if isa(Chirp2_trigCh_vec,'cell')
            Chirp2_trigCh_vec = cell2mat(Chirp2_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        Chirp2_min_trigCh_vec        = min(Chirp2_trigCh_vec);
        Chirp2_max_trigCh_vec        = max(Chirp2_trigCh_vec);
        Chirp2_trigThreshFac         = 0.05;
        Chirp2_trigHigh_vec          = double(Chirp2_trigCh_vec > Chirp2_min_trigCh_vec + Chirp2_trigThreshFac*(Chirp2_max_trigCh_vec-Chirp2_min_trigCh_vec));
        [~,Chirp2_trig_index_vec]    = findpeaks(Chirp2_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        Chirp2_trig_index_vec_length = length(Chirp2_trig_index_vec);
        
        Chirp2_NumTrigRep_mod = mod(Chirp2_trig_index_vec_length,Chirp2_NumTrigPerRep); % New: 28,05,2021

        if Chirp2_NumTrigRep_mod ~= 0 % New: 28,05,2021 (whole if statment (last line was already there w/o if statement))
            Chirp2_trig_index_vec = Chirp2_trig_index_vec(1:end-Chirp2_NumTrigRep_mod);
            Chirp2_trig_index_vec_length = length(Chirp2_trig_index_vec);
        end
        Chirp2_NumTrigRep = Chirp2_trig_index_vec_length/Chirp2_NumTrigPerRep;
        
        Chirp2_sampling_freq      = Ch_new.SamplingFrequency;
        Chirp2_sampling_int       = 1/Chirp2_sampling_freq ;
        Chirp2_trig_t_vec         = (0:1:length(Chirp2_trigCh_vec)-1)*Chirp2_sampling_int;
        
        Chirp2_trig_times_vec_temp   = Chirp2_trig_t_vec(Chirp2_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        Chirp2_trig_times_mat{i} = Chirp2_trig_times_vec_temp - Chirp2_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(Chirp2_Trig_Fig);
        subplot(Chirp2_Plot_Row,Chirp2_Plot_Col,i);
        plot(Chirp2_trig_t_vec,Chirp2_trigCh_vec); hold on;
        plot(Chirp2_trig_times_vec_temp,4085*ones(1,Chirp2_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % chirp trigger channel
        if (i/Chirp2_Plot_Col > Chirp2_Plot_Row-1) || ((i/Chirp2_Plot_Col > Chirp2_Plot_Row-2)&&(i/Chirp2_Plot_Col < Chirp2_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp2_Plot_Col)~=0)&&(mod(i,Chirp2_Plot_Col)>mod(p.Num_data_sets,Chirp2_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,Chirp2_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        
        % Find and plot the distribution of trigger intervals (not so relevant here
        % as trig times deliberately nonuniform)
        Chirp2_actual_stim_int = diff(Chirp2_trig_times_mat{i});
        Chirp2_mean_stim_int   = mean(Chirp2_actual_stim_int);
        Chirp2_stddev_stim_int = std(Chirp2_actual_stim_int);
        figure(Chirp2_Trig_Dist_Fig);
        subplot(Chirp2_Plot_Row,Chirp2_Plot_Col,i);
        histogram(Chirp2_actual_stim_int,100); hold on;
        Chirp2_v1 = vline(Chirp2_mean_stim_int,'r');
        Chirp2_v2 = vline(Chirp2_mean_stim_int-Chirp2_stddev_stim_int,'g');
        Chirp2_v3 = vline(Chirp2_mean_stim_int+Chirp2_stddev_stim_int,'g');
        if (i/Chirp2_Plot_Col > Chirp2_Plot_Row-1) || ((i/Chirp2_Plot_Col > Chirp2_Plot_Row-2)&&(i/Chirp2_Plot_Col < Chirp2_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp2_Plot_Col)~=0)&&(mod(i,Chirp2_Plot_Col)>mod(p.Num_data_sets,Chirp2_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,Chirp2_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % Gratings 400px trigger int. dist.
        set(Chirp2_v1,'LineWidth',1.5);
        set(Chirp2_v2,'LineWidth',1.5);
        set(Chirp2_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        
        % I find 6 triggers - so 3 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cell's spikes into increasing time order.
        
        % 1. remove spikes before and after the stimulus;
        % This if statement not necessary as have trigger at end of each
        % stim repeat and remove incomplete repeats.
%         if Chirp2_trig_index_vec_length>Chirp2_NumTrigPerRep
%             Chirp2_last_stim_int_vec(i) = Chirp2_trig_times_mat{i}(3)-Chirp2_trig_times_mat{i}(2);
%         else % Chirp2_trig_index_vec_length<=Chirp2_NumTrigPerRep
%             Chirp2_last_stim_int_vec(i) = 67; % sec (67 sec between 2nd and 3rd trigs)
%         end
        Chirp2_spike_times_arr{i}(Chirp2_spike_times_arr{i}<=Chirp2_trig_times_mat{i}(1)) = NaN;
        Chirp2_spike_times_arr{i}(Chirp2_spike_times_arr{i}>=(Chirp2_trig_times_mat{i}(end))) = NaN;
        
        
        if Chirp2_trig_index_vec_length>Chirp2_NumTrigPerRep
            Chirp2_stim_end_time_vec(i) = Chirp2_trig_times_mat{i}(Chirp2_NumTrigPerRep); % period of a single repeat (was Chirp2_NumTrigPerRep + 1 in other cases as no trig at end of stim usually)
        else % Chirp2_trig_index_vec_length<=Chirp2_NumTrigPerRep
            Chirp2_stim_end_time_vec(i) = 102; % sec
        end
        
        % 1.5 Calculate Quality Indices 
        Chirp2_Qual_ksdensity_grid  = linspace(0,Chirp2_stim_end_time_vec(i),Chirp2_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(Chirp2_Qual_length_ksdensity_grid,Chirp2_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = Chirp2_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/Chirp2_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:Chirp2_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),Chirp2_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*Chirp2_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,Chirp2_Qual_ksdensity_grid,'Bandwidth',Chirp2_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                Chirp2_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                Chirp2_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % Chirp2_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(Chirp2_QI_vec);
        %figure;
        %hist(Chirp2_QI_vec);
        
        % Plot the distribution of QIs
        figure(Chirp2_QI_Dist_Fig);
        subplot(Chirp2_Plot_Row,Chirp2_Plot_Col,i);
        if i == 1
            Plot_vec_loop = Chirp2_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = Chirp2_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        Chirp2_QI_v1 = vline(Chirp2_QI_thresh,'r','QI thresh.');
        Chirp2_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/Chirp2_Plot_Col > Chirp2_Plot_Row-1) || ((i/Chirp2_Plot_Col > Chirp2_Plot_Row-2)&&(i/Chirp2_Plot_Col < Chirp2_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp2_Plot_Col)~=0)&&(mod(i,Chirp2_Plot_Col)>mod(p.Num_data_sets,Chirp2_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,Chirp2_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(Chirp2_QI_v1,'LineWidth',1.5);
        set(Chirp2_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if Chirp2_trig_index_vec_length>Chirp2_NumTrigPerRep
            %Chirp2_stim_end_time_vec(i) = Chirp2_trig_times_mat{i}(Chirp2_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            Chirp2_spike_times_arr{i} = mod(Chirp2_spike_times_arr{i},Chirp2_stim_end_time_vec(i));
        %else % Chirp2_trig_index_vec_length<=Chirp2_NumTrigPerRep % --> commented as now use above
           %Chirp2_stim_end_time_vec(i) = 77; % sec % --> commented as now use above
        end
        
        % 3. sort each cell's spikes into increasing time order.
        Chirp2_spike_times_arr{i} = sort(Chirp2_spike_times_arr{i},1);
        
    end
    
    figure(Chirp2_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp2 trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Chirp2_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp2 trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Chirp2_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp2 QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(Chirp2_Trig_Fig,'color','w');
    set(Chirp2_Trig_Dist_Fig,'color','w');
    set(Chirp2_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(Chirp2_QI_vec); hold on;
    Chirp2_QI_v1 = vline(Chirp2_QI_thresh,'r','QI thresh.');
    Chirp2_QI_v2 = vline(max(Chirp2_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('Chirp2 QI - All Cells');
    set(Chirp2_QI_v1,'LineWidth',1.5);
    set(Chirp2_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % Calculate Chirp2 QI logical vec
    Chirp2_QI_logical_vec = Chirp2_QI_vec > Chirp2_QI_thresh;
    Num_Chirp2_Qual_Cells = sum(Chirp2_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_Chirp2_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,Chirp2_QI_logical_vec],[],2);
    
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1 ---> replace Chirp2_mean_stim_int_vec(j) with
    %%% corrected as diff trig ints
    for j = 2:p.Num_data_sets
        for i = 1:Chirp2_NumTrigPerRep
%             if i < Chirp2_NumTrigPerRep % This is always true when last
%             trigger is at end of stim
                indices_loop = find((Chirp2_spike_times_arr{j}>=Chirp2_trig_times_mat{j}(i))&(Chirp2_spike_times_arr{j}<Chirp2_trig_times_mat{j}(i+1)));
                b_loop       = Chirp2_trig_times_mat{1}(i+1);
                b_dash_loop  = Chirp2_trig_times_mat{j}(i+1);
%             else %i == Chirp2_NumTrigPerRep
%                 if length(Chirp2_trig_times_mat{j}) > Chirp2_NumTrigPerRep
%                     indices_loop = find((Chirp2_spike_times_arr{j}>=Chirp2_trig_times_mat{j}(i))&(Chirp2_spike_times_arr{j}<=Chirp2_trig_times_mat{j}(i+1)));
%                     b_dash_loop  = Chirp2_trig_times_mat{j}(i+1);
%                 else
%                     indices_loop = find((Chirp2_spike_times_arr{j}>=Chirp2_trig_times_mat{j}(i))&(Chirp2_spike_times_arr{j}<=Chirp2_trig_times_mat{j}(i) + Chirp2_last_stim_int_vec(j)));
%                     b_dash_loop  = Chirp2_trig_times_mat{j}(i) + Chirp2_last_stim_int_vec(j);
%                 end
%                 if length(Chirp2_trig_times_mat{1}) > Chirp2_NumTrigPerRep % was Chirp2_trig_times_mat{j} --> changed 08,04,2021
%                     b_loop       = Chirp2_trig_times_mat{1}(i+1);
%                 else
%                     b_loop       = Chirp2_trig_times_mat{1}(i) + Chirp2_last_stim_int_vec(1);
%                 end
%             end
            a_loop       = Chirp2_trig_times_mat{1}(i);
            a_dash_loop  = Chirp2_trig_times_mat{j}(i);
            Chirp2_spike_times_arr{j}(indices_loop) = a_loop + (Chirp2_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    Chirp2_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        Chirp2_num_spikes_vec(i) = size(Chirp2_spike_times_arr{i},1);
    end
    Chirp2_max_num_spikes = max(Chirp2_num_spikes_vec);
    
    cl_var.Chirp2_spike_times_mat = NaN(Chirp2_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.Chirp2_spike_times_mat(1:Chirp2_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = Chirp2_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
    %    cl_var.Chirp2_spike_times_mat = ;                         % Just above.
    %    cl_var.Chirp2_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
    cl_var.Chirp2_trig_times_vec  = Chirp2_trig_times_mat{1};    % As map subsequent data set times back to first data set times
    cl_var.Chirp2_stim_end_time   = Chirp2_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
    
end


%%% Chirp3
if p.Obj_clust_vec(12) == 1 || p.Obj_plot_vec(14) == 1
    
    Chirp3_NumTrigPerRep     = 4; % User defined
    Chirp3_spike_times_arr   = cell(p.Num_data_sets,1);
    Chirp3_trig_times_mat    = cell(p.Num_data_sets,1);
    Chirp3_stim_end_time_vec = NaN(p.Num_data_sets,1);
    %Chirp3_mean_stim_int_vec = NaN(p.Num_data_sets,1);
    Chirp3_last_stim_int_vec = NaN(p.Num_data_sets,1);
    
    Chirp3_Plot_Dim = ceil(sqrt(p.Num_data_sets));
    Chirp3_Plot_Row = ceil(p.Num_data_sets/Chirp3_Plot_Dim);
    Chirp3_Plot_Col = Chirp3_Plot_Dim;
    
    Chirp3_QI_vec                     = NaN(Total_Num_Cells,1);
    Chirp3_Qual_length_ksdensity_grid = 1000; % 100 coule use 1e3
    Chirp3_Qual_ksdensity_bdwth       = 5*1e-2; % Tom said between 1e-2 and 1e-1 is best for Orig Chirp
    
    Chirp3_Trig_Fig      = figure;
    Chirp3_Trig_Dist_Fig = figure;
    Chirp3_QI_Dist_Fig   = figure;
    
    for i = 1:p.Num_data_sets
        
        load(Chirp3_DataFiles{i});
        
        % Spike times
        if Cell_Choice == 1 % All cells
            Chirp3_spike_times_arr{i}   = spiketimestamps;
            %Chirp3_Cell_Choice_vec      = 1:1:size(spiketimestamps,2);
            %Chirp3_True_cell_index_vec = Chirp3_Cell_Choice_vec(any(~isnan(Chirp3_spike_times_arr{i}))); % Record indices of responsive cells
        else % Cell_Choice == 2 % Subset of cells (use 'Cell_Choice_vec' from above - don't define Chirp3 version)
            Chirp3_spike_times_arr{i}   = spiketimestamps(:,Cell_Choice_vec);
            %Chirp3_True_cell_index_vec = Cell_Choice_vec(any(~isnan(Chirp3_spike_times_arr{i}))); % Record indices of responsive cells
        end
        
        % Remove NaN columns (non-responsive cells)?
        %Chirp3_spike_times_arr{i}(:,~any(~isnan(Chirp3_spike_times_arr{i}))) = [];
        %cl_var.Chirp3_True_Num_Cells = size(Chirp3_spike_times_arr{i},2);
        
        % Triggers
        Chirp3_trigCh_vec            = Ch_new.trigger_ch;
        if isa(Chirp3_trigCh_vec,'cell')
            Chirp3_trigCh_vec = cell2mat(Chirp3_trigCh_vec); % 12,08,2021 Extra line when in cell form
        end
        Chirp3_min_trigCh_vec        = min(Chirp3_trigCh_vec);
        Chirp3_max_trigCh_vec        = max(Chirp3_trigCh_vec);
        Chirp3_trigThreshFac         = 0.05;
        Chirp3_trigHigh_vec          = double(Chirp3_trigCh_vec > Chirp3_min_trigCh_vec + Chirp3_trigThreshFac*(Chirp3_max_trigCh_vec-Chirp3_min_trigCh_vec));
        [~,Chirp3_trig_index_vec]    = findpeaks(Chirp3_trigHigh_vec); % 'MinPeakProminence',#,'MinPeakDistance',#
        Chirp3_trig_index_vec_length = length(Chirp3_trig_index_vec);
        
        Chirp3_NumTrigRep_mod = mod(Chirp3_trig_index_vec_length,Chirp3_NumTrigPerRep); % New: 28,05,2021

        if Chirp3_NumTrigRep_mod ~= 0 % New: 28,05,2021 (whole if statment (last line was already there w/o if statement))
            Chirp3_trig_index_vec = Chirp3_trig_index_vec(1:end-Chirp3_NumTrigRep_mod);
            Chirp3_trig_index_vec_length = length(Chirp3_trig_index_vec);
        end
        Chirp3_NumTrigRep = Chirp3_trig_index_vec_length/Chirp3_NumTrigPerRep;
        
        Chirp3_sampling_freq      = Ch_new.SamplingFrequency;
        Chirp3_sampling_int       = 1/Chirp3_sampling_freq ;
        Chirp3_trig_t_vec         = (0:1:length(Chirp3_trigCh_vec)-1)*Chirp3_sampling_int;
        
        Chirp3_trig_times_vec_temp   = Chirp3_trig_t_vec(Chirp3_trig_index_vec); % NB: this gives the same answer as the way I calculate it for RF ident.
        Chirp3_trig_times_mat{i} = Chirp3_trig_times_vec_temp - Chirp3_trig_times_vec_temp(1);
        
        % Plot trigger channel and located triggers
        figure(Chirp3_Trig_Fig);
        subplot(Chirp3_Plot_Row,Chirp3_Plot_Col,i);
        plot(Chirp3_trig_t_vec,Chirp3_trigCh_vec); hold on;
        plot(Chirp3_trig_times_vec_temp,4085*ones(1,Chirp3_trig_index_vec_length),'ro');
        title(Gen_Data_Names_vec{i}); % chirp trigger channel
        if (i/Chirp3_Plot_Col > Chirp3_Plot_Row-1) || ((i/Chirp3_Plot_Col > Chirp3_Plot_Row-2)&&(i/Chirp3_Plot_Col < Chirp3_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp3_Plot_Col)~=0)&&(mod(i,Chirp3_Plot_Col)>mod(p.Num_data_sets,Chirp3_Plot_Col)))
            xlabel('time (sec)');
        end
        if mod(i,Chirp3_Plot_Col) == 1
            ylabel('trigger channel');
        end
        set(gca,'FontSize',12);
        
        
        % Find and plot the distribution of trigger intervals (not so relevant here
        % as trig times deliberately nonuniform)
        Chirp3_actual_stim_int = diff(Chirp3_trig_times_mat{i});
        Chirp3_mean_stim_int   = mean(Chirp3_actual_stim_int);
        Chirp3_stddev_stim_int = std(Chirp3_actual_stim_int);
        figure(Chirp3_Trig_Dist_Fig);
        subplot(Chirp3_Plot_Row,Chirp3_Plot_Col,i);
        histogram(Chirp3_actual_stim_int,100); hold on;
        Chirp3_v1 = vline(Chirp3_mean_stim_int,'r');
        Chirp3_v2 = vline(Chirp3_mean_stim_int-Chirp3_stddev_stim_int,'g');
        Chirp3_v3 = vline(Chirp3_mean_stim_int+Chirp3_stddev_stim_int,'g');
        if (i/Chirp3_Plot_Col > Chirp3_Plot_Row-1) || ((i/Chirp3_Plot_Col > Chirp3_Plot_Row-2)&&(i/Chirp3_Plot_Col < Chirp3_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp3_Plot_Col)~=0)&&(mod(i,Chirp3_Plot_Col)>mod(p.Num_data_sets,Chirp3_Plot_Col)))
            xlabel('trigger int.');
        end
        if mod(i,Chirp3_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i}); % Gratings 400px trigger int. dist.
        set(Chirp3_v1,'LineWidth',1.5);
        set(Chirp3_v2,'LineWidth',1.5);
        set(Chirp3_v3,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        clear spiketimestamps Ch_new; % Remove original data files as no longer needed.
        
        
        % I find 12 triggers - so 3 repeats.
        
        % Across the matrix:
        % 1. remove spikes before and after the stimulus;
        % 2. map spike in repeats into the first repeat;
        % 3. sort each cell's spikes into increasing time order.
        
        % 1. remove spikes before and after the stimulus;
        % This if statement not necessary as have trigger at end of each
        % stim repeat and remove incomplete repeats.
%         if Chirp3_trig_index_vec_length>Chirp3_NumTrigPerRep
%             Chirp3_last_stim_int_vec(i) = Chirp3_trig_times_mat{i}(3)-Chirp3_trig_times_mat{i}(2);
%         else % Chirp3_trig_index_vec_length<=Chirp3_NumTrigPerRep
%             Chirp3_last_stim_int_vec(i) = 67; % sec (67 sec between 2nd and 3rd trigs)
%         end
        Chirp3_spike_times_arr{i}(Chirp3_spike_times_arr{i}<=Chirp3_trig_times_mat{i}(1)) = NaN;
        Chirp3_spike_times_arr{i}(Chirp3_spike_times_arr{i}>=(Chirp3_trig_times_mat{i}(end))) = NaN;
        
        
        if Chirp3_trig_index_vec_length>Chirp3_NumTrigPerRep
            Chirp3_stim_end_time_vec(i) = Chirp3_trig_times_mat{i}(Chirp3_NumTrigPerRep); % period of a single repeat (was Chirp3_NumTrigPerRep + 1 in other cases as no trig at end of stim usually)
        else % Chirp3_trig_index_vec_length<=Chirp3_NumTrigPerRep
            Chirp3_stim_end_time_vec(i) = 102; % sec
        end
        
        % 1.5 Calculate Quality Indices 
        Chirp3_Qual_ksdensity_grid  = linspace(0,Chirp3_stim_end_time_vec(i),Chirp3_Qual_length_ksdensity_grid);
        for j = 1:Num_Cell_Per_Exp_vec(i)
            response_mat         = zeros(Chirp3_Qual_length_ksdensity_grid,Chirp3_NumTrigRep); % time samples x stimulus repetitions
            Spike_times_vec_loop = Chirp3_spike_times_arr{i}(:,j);                          % Spike times for the jth cell of the ith experiment
            index_loop           = ceil(Spike_times_vec_loop/Chirp3_stim_end_time_vec(i));  % Which repetition of the stimulus
            for k = 1:Chirp3_NumTrigRep    % Loop over stimulus repetitions
                if sum(index_loop==k)>0 % Check there are spikes in the given repetition
                    stim_times_loop       = mod(Spike_times_vec_loop(index_loop==k),Chirp3_stim_end_time_vec(i));                % Find stimulus times for this repetition and map onto initial time window
                    vec_loop              = [-stim_times_loop;stim_times_loop;(2*Chirp3_stim_end_time_vec(i)-stim_times_loop)];  % Create vector with reflected spike pattern to left and right
                    [response_mat(:,k),~] = ksdensity(vec_loop,Chirp3_Qual_ksdensity_grid,'Bandwidth',Chirp3_Qual_ksdensity_bdwth); % Calculate ksdensiy smoothed reposnse for each repeat
                end
            end
            if i == 1
                Chirp3_QI_vec(j)                                  = var(mean(response_mat,2))/mean(var(response_mat,0,1));
            else % i>1
                Chirp3_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+j) = var(mean(response_mat,2))/mean(var(response_mat,0,1)); % Chirp3_QI_vec(Num_Cell_Per_Exp_vec(i-1)+j)
            end
        end
        % Check output
        %figure;
        %plot(response_mat); % response_mat(:,1:2)
        %figure;
        %plot(Chirp3_QI_vec);
        %figure;
        %hist(Chirp3_QI_vec);
        
        % Plot the distribution of QIs
        figure(Chirp3_QI_Dist_Fig);
        subplot(Chirp3_Plot_Row,Chirp3_Plot_Col,i);
        if i == 1
            Plot_vec_loop = Chirp3_QI_vec(1:Num_Cell_Per_Exp_vec(1));
        else
            Plot_vec_loop = Chirp3_QI_vec(sum(Num_Cell_Per_Exp_vec(1:i-1))+1:sum(Num_Cell_Per_Exp_vec(1:i)));
        end
        histogram(Plot_vec_loop); hold on;
        Chirp3_QI_v1 = vline(Chirp3_QI_thresh,'r','QI thresh.');
        Chirp3_QI_v2 = vline(max(Plot_vec_loop),'g','max QI');
        if (i/Chirp3_Plot_Col > Chirp3_Plot_Row-1) || ((i/Chirp3_Plot_Col > Chirp3_Plot_Row-2)&&(i/Chirp3_Plot_Col < Chirp3_Plot_Row-1)&&(mod(p.Num_data_sets,Chirp3_Plot_Col)~=0)&&(mod(i,Chirp3_Plot_Col)>mod(p.Num_data_sets,Chirp3_Plot_Col)))
            xlabel('QI');
        end
        if mod(i,Chirp3_Plot_Col) == 1
            ylabel('freq.');
        end
        title(Gen_Data_Names_vec{i});
        set(Chirp3_QI_v1,'LineWidth',1.5);
        set(Chirp3_QI_v2,'LineWidth',1.5);
        set(gca,'FontSize',12);
        
        
        
        
        % 2. map spike in repeats into the first repeat;
        % Use modular arithmetic to map spikes in repeats onto first repeat
        if Chirp3_trig_index_vec_length>Chirp3_NumTrigPerRep
            %Chirp3_stim_end_time_vec(i) = Chirp3_trig_times_mat{i}(Chirp3_NumTrigPerRep+1); % period of a single repeat % --> commented as now use above
            Chirp3_spike_times_arr{i} = mod(Chirp3_spike_times_arr{i},Chirp3_stim_end_time_vec(i));
        %else % Chirp3_trig_index_vec_length<=Chirp3_NumTrigPerRep % --> commented as now use above
           %Chirp3_stim_end_time_vec(i) = 77; % sec % --> commented as now use above
        end
        
        % 3. sort each cell's spikes into increasing time order.
        Chirp3_spike_times_arr{i} = sort(Chirp3_spike_times_arr{i},1);
        
    end
    
    figure(Chirp3_Trig_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp3 trigger channel','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Chirp3_Trig_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp3 trigger int. dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    figure(Chirp3_QI_Dist_Fig);
    annotation('textbox',[.5 .975 0 0],'String','Chirp3 QI dist.','FitBoxToText','on','FontSize',12,'FontWeight','bold','LineStyle','none','HorizontalAlignment','center','VerticalAlignment','middle');
    % By default, the units are normalized to the figure. The lower left corner of the figure maps to
    % (0,0) and the upper right corner maps to (1,1). dim - [x y w h].
    set(Chirp3_Trig_Fig,'color','w');
    set(Chirp3_Trig_Dist_Fig,'color','w');
    set(Chirp3_QI_Dist_Fig,'color','w');
    
    % Plot the distribution of QIs for the full data set
    figure;
    histogram(Chirp3_QI_vec); hold on;
    Chirp3_QI_v1 = vline(Chirp3_QI_thresh,'r','QI thresh.');
    Chirp3_QI_v2 = vline(max(Chirp3_QI_vec),'g','max QI');
    xlabel('QI');
    ylabel('freq.');
    title('Chirp3 QI - All Cells');
    set(Chirp3_QI_v1,'LineWidth',1.5);
    set(Chirp3_QI_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % Calculate Chirp3 QI logical vec
    Chirp3_QI_logical_vec = Chirp3_QI_vec > Chirp3_QI_thresh;
    Num_Chirp3_Qual_Cells = sum(Chirp3_QI_logical_vec);
    
    % Update Final number of cells per stimulus vector (after poor quality data removed)
    cl_var.Num_Cell_Per_Stim_vec = [cl_var.Num_Cell_Per_Stim_vec,Num_Chirp3_Qual_Cells];
    
    % Update Overall Quality Vector
    Overall_Quality_vec = max([Overall_Quality_vec,Chirp3_QI_logical_vec],[],2);
    
    
    %%% Map spikes in data sets 2 onwards onto the same inter-trigger intervals
    %%% as in data set 1 ---> replace Chirp3_mean_stim_int_vec(j) with
    %%% corrected as diff trig ints
    for j = 2:p.Num_data_sets
        for i = 1:Chirp3_NumTrigPerRep
%             if i < Chirp3_NumTrigPerRep % This is always true when last
%             trigger is at end of stim
                indices_loop = find((Chirp3_spike_times_arr{j}>=Chirp3_trig_times_mat{j}(i))&(Chirp3_spike_times_arr{j}<Chirp3_trig_times_mat{j}(i+1)));
                b_loop       = Chirp3_trig_times_mat{1}(i+1);
                b_dash_loop  = Chirp3_trig_times_mat{j}(i+1);
%             else %i == Chirp3_NumTrigPerRep
%                 if length(Chirp3_trig_times_mat{j}) > Chirp3_NumTrigPerRep
%                     indices_loop = find((Chirp3_spike_times_arr{j}>=Chirp3_trig_times_mat{j}(i))&(Chirp3_spike_times_arr{j}<=Chirp3_trig_times_mat{j}(i+1)));
%                     b_dash_loop  = Chirp3_trig_times_mat{j}(i+1);
%                 else
%                     indices_loop = find((Chirp3_spike_times_arr{j}>=Chirp3_trig_times_mat{j}(i))&(Chirp3_spike_times_arr{j}<=Chirp3_trig_times_mat{j}(i) + Chirp3_last_stim_int_vec(j)));
%                     b_dash_loop  = Chirp3_trig_times_mat{j}(i) + Chirp3_last_stim_int_vec(j);
%                 end
%                 if length(Chirp3_trig_times_mat{1}) > Chirp3_NumTrigPerRep % was Chirp3_trig_times_mat{j} --> changed 08,04,2021
%                     b_loop       = Chirp3_trig_times_mat{1}(i+1);
%                 else
%                     b_loop       = Chirp3_trig_times_mat{1}(i) + Chirp3_last_stim_int_vec(1);
%                 end
%             end
            a_loop       = Chirp3_trig_times_mat{1}(i);
            a_dash_loop  = Chirp3_trig_times_mat{j}(i);
            Chirp3_spike_times_arr{j}(indices_loop) = a_loop + (Chirp3_spike_times_arr{j}(indices_loop) - a_dash_loop)*(b_loop - a_loop)/(b_dash_loop - a_dash_loop);
        end
    end
    
    %%% Combine spike times array entries into a single matrix
    Chirp3_num_spikes_vec = NaN(p.Num_data_sets,1);
    for i = 1:p.Num_data_sets
        Chirp3_num_spikes_vec(i) = size(Chirp3_spike_times_arr{i},1);
    end
    Chirp3_max_num_spikes = max(Chirp3_num_spikes_vec);
    
    cl_var.Chirp3_spike_times_mat = NaN(Chirp3_max_num_spikes,Total_Num_Cells);
    loop_var = 0;
    for i = 1:p.Num_data_sets
        cl_var.Chirp3_spike_times_mat(1:Chirp3_num_spikes_vec(i),loop_var+1:loop_var+Num_Cell_Per_Exp_vec(i)) = Chirp3_spike_times_arr{i};
        loop_var = loop_var + Num_Cell_Per_Exp_vec(i);
    end
    
    %    cl_var.Chirp3_spike_times_mat = ;                         % Just above.
    %    cl_var.Chirp3_True_Num_Cells  = ;                         % Won't use, will use same cells across all stimuli
    cl_var.Chirp3_trig_times_vec  = Chirp3_trig_times_mat{1};    % As map subsequent data set times back to first data set times
    cl_var.Chirp3_stim_end_time   = Chirp3_stim_end_time_vec(1); % As map subsequent data set times back to first data set times
    
end


% Save before applying QC

%save data_DataPrep4_Chick_07_09_2021_9_NoQC;
%save data_DataPrep5_Chick_24_09_2021_1_NoQC;
%save data_DataPrep5_Chick_24_09_2021_3_NoQC;
%save data_DataPrep5_Chick_18_11_2021_PreQC_1;
%save('data_FFF2_Cstep_Chirp3_19_11_2021_PreQC_1.mat','-v7.3');

% %%% Find the number of true cells (which spike) % ---> Commented out below when got new QIs 28,04,2021
% 
% % Find the indices of cells which don't spike for each stimulus
% No_Spike_Indices_mat = NaN(Num_Stim_Loaded,Total_Num_Cells); % Each row is for a different stimulus, the indices of cells unresponsive to that stimulus being listed.
% Total_Cell_Choice_vec = 1:1:Total_Num_Cells;
% loop_var = 1;
% if p.Obj_clust_vec(1) == 1 || p.Obj_plot_vec(3) == 1 % FFF
%     loop_vec =  Total_Cell_Choice_vec(~any(~isnan(cl_var.FFF_spike_times_mat))); % Record indices of unresponsive cells
%     No_Spike_Indices_mat(loop_var,1:length(loop_vec)) = loop_vec;
%     loop_var = loop_var + 1;
% end
% if p.Obj_clust_vec(2) == 1 || p.Obj_plot_vec(4) == 1 % Chirp
%     loop_vec =  Total_Cell_Choice_vec(~any(~isnan(cl_var.Chirp_spike_times_mat))); % Record indices of unresponsive cells
%     No_Spike_Indices_mat(loop_var,1:length(loop_vec)) = loop_vec;
%     loop_var = loop_var + 1;
% end
% if p.Obj_clust_vec(3) == 1 || p.Obj_plot_vec(5) == 1 % FFF Noise
%     loop_vec =  Total_Cell_Choice_vec(~any(~isnan(cl_var.FFF_Noise_STA(:,:,1)),2)); % Record indices of unresponsive cells
%     No_Spike_Indices_mat(loop_var,1:length(loop_vec)) = loop_vec;
%     loop_var = loop_var + 1;
% end
% if p.Obj_clust_vec(4) == 1 || p.Obj_plot_vec(6) == 1 % Gratings
%     loop_vec =  Total_Cell_Choice_vec(~any(~isnan(cl_var.Gratings_400px_spike_times_mat))); % Record indices of unresponsive cells
%     No_Spike_Indices_mat(loop_var,1:length(loop_vec)) = loop_vec;
%     loop_var = loop_var + 1;
% end
% % if sum(p.Obj_clust_vec(5:6)) > 0 || sum(p.Obj_plot_vec([1,2,7,8])) > 0 % CNoise
% %     loop_vec =  Total_Cell_Choice_vec(~any(~isnan(cl_var.CNoise_20px_spike_times_mat))); % Record indices of unresponsive cells
% %     No_Spike_Indices_mat(loop_var,1:length(loop_vec)) = loop_vec;
% % end
% 
% %No_Spike_Indices_Total_vec = unique(No_Spike_Indices_mat);
% %No_Spike_Indices_Total_vec(isnan(No_Spike_Indices_Total_vec)) = [];
% % Num_Unresponsive_Cells = length(No_Spike_Indices_Total_vec);
% Fully_Unresponsive_Cells = No_Spike_Indices_mat(1,:);
% if (Num_Stim_Loaded > 1) && (~isempty(Fully_Unresponsive_Cells))
%     for i = 2:Num_Stim_Loaded
%         if ~isempty(Fully_Unresponsive_Cells) % Fully_Unresponsive_Cells changes on each loop and can stop once empty
%             Fully_Unresponsive_Cells = intersect(Fully_Unresponsive_Cells,No_Spike_Indices_mat(i,:)); % was No_Spike_Indices_mat(2,:) --> changed 08,04,2021
%         end
%     end   
% end
% Num_Unresponsive_Cells = length(Fully_Unresponsive_Cells);
% 
% % Eliminate cells that don't spike for any one stimulus (or which don't spike for all stim) cl_var.spike_times_mat's,
% % cl_var.Exp_Ident_vec and cl_var.True_Num_Cells (also Num_Cell_Per_Exp_vec?)
% if ~isempty(Fully_Unresponsive_Cells) % Remove any cell missing from all stim, no cells removed if Fully_Unresponsive_Cells is empty.
%     
%     if p.Obj_clust_vec(1) == 1 || p.Obj_plot_vec(3) == 1 % FFF
%         %cl_var.FFF_spike_times_mat(:,No_Spike_Indices_Total_vec) = [];
%         cl_var.FFF_spike_times_mat(:,Fully_Unresponsive_Cells) = [];
%     end
%     if p.Obj_clust_vec(2) == 1 || p.Obj_plot_vec(4) == 1 % Chirp
%         %cl_var.Chirp_spike_times_mat(:,No_Spike_Indices_Total_vec) = [];
%         cl_var.Chirp_spike_times_mat(:,Fully_Unresponsive_Cells) = [];
%     end
%     if p.Obj_clust_vec(3) == 1 || p.Obj_plot_vec(5) == 1 % FFF Noise
%         cl_var.FFF_Noise_STA(Fully_Unresponsive_Cells,:,:)    = [];
%         cl_var.FFF_Noise_Full_STA(Fully_Unresponsive_Cells,:) = [];
%     end
%     if p.Obj_clust_vec(4) == 1 || p.Obj_plot_vec(6) == 1 % Gratings
%         %cl_var.Gratings_400px_spike_times_mat(:,No_Spike_Indices_Total_vec) = [];
%         cl_var.Gratings_400px_spike_times_mat(:,Fully_Unresponsive_Cells) = [];
%     end
%     if sum(p.Obj_clust_vec(5:6)) > 0 || sum(p.Obj_plot_vec([1,2,7,8])) > 0 % CNoise
%         %cl_var.CNoise_20px_spike_times_mat(:,No_Spike_Indices_Total_vec) = [];
%         % Prob diff here as kernels
% 
%     end
%     if p.Obj_clust_vec(5) == 1 || p.Obj_plot_vec(7) == 1 % CNoise: Full RF Size
%         %cl_var.Full_RF_Size_vec(No_Spike_Indices_Total_vec) = [];
%         cl_var.Full_RF_Size_vec(Fully_Unresponsive_Cells) = [];
%     end
%     if p.Obj_clust_vec(6) == 1 || p.Obj_plot_vec(8) == 1 % CNoise: Full RF Ellipticity
%         %cl_var.Full_RF_Ellipticity_vec(No_Spike_Indices_Total_vec) = [];
%         cl_var.Full_RF_Ellipticity_vec(Fully_Unresponsive_Cells) = [];
%     end
%     
%     %cl_var.Exp_Ident_vec(No_Spike_Indices_Total_vec) = [];
%     cl_var.Exp_Ident_vec(Fully_Unresponsive_Cells) = [];
%     
% end
% 
% cl_var.True_Num_Cells = Total_Num_Cells - Num_Unresponsive_Cells;  % ---> Commented out above when got new QIs 28,04,2021

%%% Remove low quality cells (with no above threshold responses for any clustered sitmulus)
% Also remove low quality cells from all plotted stimuli

if p.Obj_plot_vec(1) == 1 % Cell positions
    cl_var.Full_RF_Pos_mat                = cl_var.Full_RF_Pos_mat(Overall_Quality_vec,:);
end
if p.Obj_clust_vec(1) == 1 || p.Obj_plot_vec(3) == 1 % FFF
    cl_var.FFF_spike_times_mat            = cl_var.FFF_spike_times_mat(:,Overall_Quality_vec);
end
if p.Obj_clust_vec(2) == 1 || p.Obj_plot_vec(4) == 1 % Chirp
    cl_var.Chirp_spike_times_mat          = cl_var.Chirp_spike_times_mat(:,Overall_Quality_vec);
end
if p.Obj_clust_vec(3) == 1 || p.Obj_plot_vec(5) == 1 % FFF Noise
    cl_var.FFF_Noise_STA                  = cl_var.FFF_Noise_STA(Overall_Quality_vec,:,:);
    cl_var.FFF_Noise_Full_STA             = cl_var.FFF_Noise_Full_STA(Overall_Quality_vec,:);
    cl_var.FFF_Noise_STA_scaled           = cl_var.FFF_Noise_STA_scaled(Overall_Quality_vec,:,:);    % For normalisation code (10,06,2021)
    cl_var.FFF_Noise_Full_STA_scaled      = cl_var.FFF_Noise_Full_STA_scaled(Overall_Quality_vec,:); % For normalisation code (10,06,2021)
    %FFF_Noise_SD_logical_mat              = FFF_Noise_SD_logical_mat(Overall_Quality_vec,:); % PAR Mod 23,11,2021 - doesn't seem to be used after this
end
if p.Obj_clust_vec(4) == 1 || p.Obj_plot_vec(6) == 1 % Gratings
    cl_var.Gratings_400px_spike_times_mat = cl_var.Gratings_400px_spike_times_mat(:,Overall_Quality_vec);
end
if p.Obj_clust_vec(5) == 1 || p.Obj_plot_vec(7) == 1 % CNoise: Full RF Size
    cl_var.Full_RF_Size_vec               = cl_var.Full_RF_Size_vec(Overall_Quality_vec);
end
if p.Obj_clust_vec(6) == 1 || p.Obj_plot_vec(8) == 1 % CNoise: Full RF Ellipticity
    cl_var.Full_RF_Ellipticity_vec        = cl_var.Full_RF_Ellipticity_vec(Overall_Quality_vec);
end
if p.Obj_clust_vec(7) == 1 || p.Obj_plot_vec(9) == 1 % CNoise: Full RF Dominant Axis Angle
    cl_var.Full_RF_Dom_Ax_Ang_vec         = cl_var.Full_RF_Dom_Ax_Ang_vec(Overall_Quality_vec);
end
if p.Obj_clust_vec(8) == 1 || p.Obj_plot_vec(10) == 1 % FFF2
    cl_var.FFF2_spike_times_mat           = cl_var.FFF2_spike_times_mat(:,Overall_Quality_vec);
end
if p.Obj_clust_vec(9) == 1 || p.Obj_plot_vec(11) == 1 % Chirp2
    cl_var.Chirp2_spike_times_mat           = cl_var.Chirp2_spike_times_mat(:,Overall_Quality_vec);
end
if p.Obj_clust_vec(10) == 1 || p.Obj_plot_vec(12) == 1 % SSub
    cl_var.SSub_spike_times_mat           = cl_var.SSub_spike_times_mat(:,Overall_Quality_vec);
end
if p.Obj_clust_vec(11) == 1 || p.Obj_plot_vec(13) == 1 % CSteps
    cl_var.CSteps_spike_times_mat           = cl_var.CSteps_spike_times_mat(:,Overall_Quality_vec);
end
if p.Obj_clust_vec(12) == 1 || p.Obj_plot_vec(14) == 1 % Chirp3
    cl_var.Chirp3_spike_times_mat           = cl_var.Chirp3_spike_times_mat(:,Overall_Quality_vec);
end

% Final number of cells per experiment vector (after poor quality data removed)
cl_var.Num_Cell_Per_Exp_vec = NaN(p.Num_data_sets,1);
for i = 1:p.Num_data_sets
    if i==1
        cl_var.Num_Cell_Per_Exp_vec(i) = sum(Overall_Quality_vec(1:Num_Cell_Per_Exp_vec(1)));
    else
        cl_var.Num_Cell_Per_Exp_vec(i) = sum(Overall_Quality_vec(sum(Num_Cell_Per_Exp_vec(1:i-1)):sum(Num_Cell_Per_Exp_vec(1:i))));
    end
end

% Experiment identity vector
cl_var.Exp_Ident_vec = cl_var.Exp_Ident_vec(Overall_Quality_vec);

% Cell identity vector
cl_var.Cell_Ident_vec = cl_var.Cell_Ident_vec(Overall_Quality_vec);

% Find the number of cells responsive cells
cl_var.True_Num_Cells = sum(Overall_Quality_vec); % Num_Overall_Qual_Cells

%save('data_FFFNoise_FFF2_Cstep_Chirp3_19_11_2021_PreQC_3.mat','-v7.3'); % data 1-15 (No QC)
%save('data_FFFNoise_FFF2_Cstep_Chirp3_19_11_2021_PreQC_3_1.mat','-v7.3'); % data 1-16 (No QC)
%save('data_FFFNoise_FFF2_Cstep_Chirp3_19_11_2021_PreQC_3_2.mat','-v7.3'); % data 1-17 (loop but not FFFNoise code after) (No QC)
%save('data_FFFNoise_FFF2_Cstep_Chirp3_19_11_2021_PreQC_3_3.mat','-v7.3'); % Actually is QCed


%% Save Data

%save data_DataPrep4_Chick_22_06_2021_1;
%save data_DataPrep4_Chick_22_06_2021_2;
%save data_DataPrep4_Chick_22_06_2021_3;
%save data_DataPrep4_Chick_22_06_2021_4;
%save data_DataPrep4_Chick_22_06_2021_5;
%save data_DataPrep4_Chick_22_06_2021_6;
%save data_DataPrep4_Chick_22_06_2021_7;
%save data_DataPrep4_Chick_22_06_2021_8;
%save data_DataPrep4_Chick_22_06_2021_9;
%save data_DataPrep4_Chick_07_07_2021_1; % Accidentally saved over this 6th/7th/09/2021

%save data_DataPrep4_Chick_07_09_2021_1;
%save data_DataPrep4_Chick_07_09_2021_2;
%save data_DataPrep4_Chick_07_09_2021_3;
%save data_DataPrep4_Chick_07_09_2021_4;
%save data_DataPrep4_Chick_07_09_2021_5;
%save data_DataPrep4_Chick_07_09_2021_6;
%save data_DataPrep4_Chick_07_09_2021_7;
%save data_DataPrep4_Chick_07_09_2021_8;
%save data_DataPrep4_Chick_07_09_2021_9;

%save data_DataPrep5_Chick_24_09_2021_1;
%save data_DataPrep5_Chick_24_09_2021_2;
%save data_DataPrep5_Chick_24_09_2021_3;
%save data_DataPrep5_Chick_24_09_2021_3_AndFNoise;

%save data_DataPrep5_Chick_04_11_2021_1;

%save data_DataPrep5_Chick_15_11_2021_PreQC_1;
%save data_DataPrep5_Chick_15_11_2021_1;

%save data_DataPrep5_Chick_18_11_2021_1;

%save('data_FFFNoise_FFF2_Cstep_Chirp3_19_11_2021_QCed_1.mat','-v7.3');
%save('data_FFFNoise_FFF2_Cstep_Chirp3_19_11_2021_QCed_2.mat','-v7.3');

%save('data_FFFNoise_29_11_2021_QCed_1.mat','-v7.3');

%save('data_FFFNoise_FFF2_Cstep_Chirp3_19_11_2021_PreQC_3_3_Mod.mat','-v7.3');

%% Call Clustering Fn

%[] = CL_RF_Clustering_6_fn(cl_var,p);





