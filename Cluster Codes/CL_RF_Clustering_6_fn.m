function [] = CL_RF_Clustering_6_fn(cl_var,p)
% RF Clustering 6
% Function to ready data for clustering and call clustering function
% 23,09,2021 onwards: modified version of CL_RF_Clustering_5_fn, when
% added Chirp3 stimulus

%% Ready Data For Clustering

% Preallocate
All_Scores     = [];
Total_Num_Comp = 0;

%%% FFF
if p.Obj_clust_vec(1) == 1 || p.Obj_plot_vec(3) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    FFF_length_ksdensity_grid = 1e3; % now: 1e3, was: 100 choose a number divisible by 4 (200?)
    FFF_ksdensity_grid        = linspace(0,cl_var.FFF_stim_end_time,FFF_length_ksdensity_grid);
    FFF_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best
    FFF_spike_density_mat     = NaN(cl_var.True_Num_Cells,FFF_length_ksdensity_grid);%FFF_length_ksdensity_grid,cl_var.True_Num_Cells
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.FFF_spike_times_mat(:,i);cl_var.FFF_spike_times_mat(:,i);(2*cl_var.FFF_stim_end_time-cl_var.FFF_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            [FFF_spike_density_mat(i,:),~] = ksdensity(vec_loop,FFF_ksdensity_grid,'Bandwidth',FFF_ksdensity_bdwth); % [FFF_spike_density_mat(:,i),~] ksdensity(cl_var.FFF_spike_times_mat(:,i),FFF_ksdensity_grid,'Bandwidth',FFF_ksdensity_bdwth)
        else % If there were no spikes for this cell for this stim.
            FFF_spike_density_mat(i,:) = zeros(1,FFF_length_ksdensity_grid);
        end
    end
    % If want to limit to interval:
    % 'Support',[0,cl_var.FFF_stim_end_time],'BoundaryCorrection','log'/'reflection'(no good as makes go to zero at ends)
    
    % 1D plot of chosen cell FFF spike rate
    figure;
    plot(FFF_ksdensity_grid,FFF_spike_density_mat(2,:),'LineWidth',1.5); hold on; % FFF_spike_density_mat(:,2)
    plot(cl_var.FFF_spike_times_mat(:,2),zeros(size(cl_var.FFF_spike_times_mat,1),1),'bx');
    FFF_ks_v1 = vline(cl_var.FFF_trig_times_vec(2),'r');
    FFF_ks_v2 = vline(cl_var.FFF_trig_times_vec(3),'r');
    FFF_ks_v3 = vline(cl_var.FFF_trig_times_vec(4),'r');
    xlim([0 cl_var.FFF_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('FFF spike rate');
    set(FFF_ks_v1,'LineWidth',1.5);
    set(FFF_ks_v2,'LineWidth',1.5);
    set(FFF_ks_v3,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell FFF spike rates
    FFF_ks_fig = figure;
    imagesc(FFF_ksdensity_grid,[],FFF_spike_density_mat); hold on; % FFF_spike_density_mat'
    colormap(FFF_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    FFF_ks_v1 = vline(cl_var.FFF_trig_times_vec(2),'r');
    FFF_ks_v2 = vline(cl_var.FFF_trig_times_vec(3),'r');
    FFF_ks_v3 = vline(cl_var.FFF_trig_times_vec(4),'r');
    xlabel('time (sec)');
    ylabel('RF');
    title('FFF spike rates');
    set(FFF_ks_v1,'LineWidth',1.5);
    set(FFF_ks_v2,'LineWidth',1.5);
    set(FFF_ks_v3,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(1) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            %FFF_segment_time = 0.25*cl_var.FFF_stim_end_time; % Could have used this...
            
            FFF_R_spike_density_mat  = FFF_spike_density_mat(:,1:0.25*FFF_length_ksdensity_grid);
            FFF_G_spike_density_mat  = FFF_spike_density_mat(:,0.25*FFF_length_ksdensity_grid+1:0.5*FFF_length_ksdensity_grid);
            FFF_B_spike_density_mat  = FFF_spike_density_mat(:,0.5*FFF_length_ksdensity_grid+1:0.75*FFF_length_ksdensity_grid);
            FFF_UV_spike_density_mat = FFF_spike_density_mat(:,0.75*FFF_length_ksdensity_grid+1:FFF_length_ksdensity_grid);
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [FFF_coeff,FFF_score,FFF_latent,FFF_tsquared,FFF_explained,FFF_mu] = pca(FFF_spike_density_mat); %  FFF_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                FFF_explained_cumsum  = cumsum(FFF_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                FFF_PC_Num = find(FFF_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                FFF_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,FFF_score(:,1:FFF_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + FFF_PC_Num;
            
            disp(sprintf('Num FFF PCs = %i',FFF_PC_Num));
            
            FFF_pca_fig = figure;
            imagesc(FFF_ksdensity_grid,[],abs(FFF_coeff(:,1:FFF_PC_Num)'));%FFF_PC_Num, 10
            colormap(FFF_pca_fig,gray(256));
            colorbar;
            FFF_ks_v1 = vline(cl_var.FFF_trig_times_vec(2),'r');
            FFF_ks_v2 = vline(cl_var.FFF_trig_times_vec(3),'r');
            FFF_ks_v3 = vline(cl_var.FFF_trig_times_vec(4),'r');
            title('FFF - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(FFF_ks_v1,'LineWidth',1.5);
            set(FFF_ks_v2,'LineWidth',1.5);
            set(FFF_ks_v3,'LineWidth',1.5);
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [FFF_R_coeff,FFF_R_score,FFF_R_latent,FFF_R_tsquared,FFF_R_explained,FFF_R_mu]       = pca(FFF_R_spike_density_mat);
            [FFF_G_coeff,FFF_G_score,FFF_G_latent,FFF_G_tsquared,FFF_G_explained,FFF_G_mu]       = pca(FFF_G_spike_density_mat);
            [FFF_B_coeff,FFF_B_score,FFF_B_latent,FFF_B_tsquared,FFF_B_explained,FFF_B_mu]       = pca(FFF_B_spike_density_mat);
            [FFF_UV_coeff,FFF_UV_score,FFF_UV_latent,FFF_UV_tsquared,FFF_UV_explained,FFF_UV_mu] = pca(FFF_UV_spike_density_mat);
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                FFF_R_explained_cumsum  = cumsum(FFF_R_explained);
                FFF_G_explained_cumsum  = cumsum(FFF_G_explained);
                FFF_B_explained_cumsum  = cumsum(FFF_B_explained);
                FFF_UV_explained_cumsum = cumsum(FFF_UV_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                FFF_R_PC_Num  = find(FFF_R_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF_G_PC_Num  = find(FFF_G_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF_B_PC_Num  = find(FFF_B_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF_UV_PC_Num = find(FFF_UV_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                FFF_R_PC_Num  = p.PCA_comp_num;
                FFF_G_PC_Num  = p.PCA_comp_num;
                FFF_B_PC_Num  = p.PCA_comp_num;
                FFF_UV_PC_Num = p.PCA_comp_num;
                
            end
            
            FFF_PC_Num = FFF_R_PC_Num + FFF_G_PC_Num + FFF_B_PC_Num + FFF_UV_PC_Num;
            
            All_Scores     = [All_Scores,FFF_R_score(:,1:FFF_R_PC_Num),FFF_G_score(:,1:FFF_G_PC_Num),...
                FFF_B_score(:,1:FFF_B_PC_Num),FFF_UV_score(:,1:FFF_UV_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + FFF_PC_Num;
            
            disp(sprintf('Num FFF R PCs = %i',FFF_R_PC_Num));
            disp(sprintf('Num FFF G PCs = %i',FFF_G_PC_Num));
            disp(sprintf('Num FFF B PCs = %i',FFF_B_PC_Num));
            disp(sprintf('Num FFF UV PCs = %i',FFF_UV_PC_Num));
            disp(sprintf('Total Num FFF PCs = %i',FFF_PC_Num));
            
            FFF_pca_fig = figure;
            subplot(2,2,1);
            imagesc(FFF_ksdensity_grid(1:0.25*FFF_length_ksdensity_grid),[],abs(FFF_R_coeff(:,1:FFF_R_PC_Num)'));
            colormap(FFF_pca_fig,gray(256));
            colorbar;
            title('FFF R - PCA coeff');ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,2,2);
            imagesc(FFF_ksdensity_grid(0.25*FFF_length_ksdensity_grid+1:0.5*FFF_length_ksdensity_grid),[],abs(FFF_G_coeff(:,1:FFF_G_PC_Num)'));
            colormap(FFF_pca_fig,gray(256));
            colorbar;
            title('FFF G - PCA coeff');
            set(gca,'FontSize',12);
            subplot(2,2,3);
            imagesc(FFF_ksdensity_grid(0.5*FFF_length_ksdensity_grid+1:0.75*FFF_length_ksdensity_grid),[],abs(FFF_B_coeff(:,1:FFF_B_PC_Num)'));
            colormap(FFF_pca_fig,gray(256));
            colorbar;
            title('FFF B - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,2,4);
            imagesc(FFF_ksdensity_grid(0.75*FFF_length_ksdensity_grid+1:FFF_length_ksdensity_grid),[],abs(FFF_UV_coeff(:,1:FFF_UV_PC_Num)'));
            colormap(FFF_pca_fig,gray(256));
            colorbar;
            title('FFF UV - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% Chirp
if p.Obj_clust_vec(2) == 1 || p.Obj_plot_vec(4) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    Chirp_length_ksdensity_grid = 1e3; % now: 1e3, was: 100
    Chirp_ksdensity_grid        = linspace(0,cl_var.Chirp_stim_end_time,Chirp_length_ksdensity_grid);
    Chirp_ksdensity_bdwth       = 5*1e-2; % Tom said between 1e-2 and 1e-1 is best
    Chirp_spike_density_mat     = NaN(cl_var.True_Num_Cells,Chirp_length_ksdensity_grid);
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.Chirp_spike_times_mat(:,i);cl_var.Chirp_spike_times_mat(:,i);(2*cl_var.Chirp_stim_end_time-cl_var.Chirp_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            %[Chirp_spike_density_mat(i,:),~] = ksdensity(vec_loop,Chirp_ksdensity_grid,'Bandwidth',Chirp_ksdensity_bdwth); % ksdensity(cl_var.Chirp_spike_times_mat(:,i),Chirp_ksdensity_grid,'Bandwidth',Chirp_ksdensity_bdwth);
            [Chirp_spike_density_temp,~] = ksdensity(vec_loop,Chirp_ksdensity_grid,'Bandwidth',Chirp_ksdensity_bdwth); % Mod 14,06,2021
            %Chirp_spike_density_mat(i,:) = Chirp_spike_density_temp/max(Chirp_spike_density_temp);                     % Mod 14,06,2021
            Chirp_spike_density_mat(i,:) = Chirp_spike_density_temp*sum(~isnan(cl_var.Chirp_spike_times_mat(:,i)));                     % Mod 14,06,2021
        else % If there were no spikes for this cell for this stim.
            Chirp_spike_density_mat(i,:) = zeros(1,Chirp_length_ksdensity_grid);
        end
    end
    
    % Not necc as have 0 cl_var.Chirp_stim_end_time
    % % Spike density time start and end
    % Chirp_t_start = Chirp_ksdensity_grid(1);
    % Chirp_t_end   = Chirp_ksdensity_grid(end);
    
    % 1D plot of chosen cell Chirp spike rate
    figure;
    plot(Chirp_ksdensity_grid,Chirp_spike_density_mat(2,:),'LineWidth',1.5); hold on;
    plot(cl_var.Chirp_spike_times_mat(:,2),zeros(size(cl_var.Chirp_spike_times_mat,1),1),'bx');
    Chirp_ks_v1 = vline(cl_var.Chirp_trig_times_vec(2),'r');
    xlim([0 cl_var.Chirp_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('Chirp spike rate');
    set(Chirp_ks_v1,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell Chirp spike rates
    Chirp_ks_fig = figure;
    imagesc(Chirp_ksdensity_grid,[],Chirp_spike_density_mat); hold on;
    colormap(Chirp_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    Chirp_ks_v1 = vline(cl_var.Chirp_trig_times_vec(2),'r');
    xlabel('time (sec)');
    ylabel('RF');
    title('Chirp spike rates');
    set(Chirp_ks_v1,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(2) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            Num_Chirp_1_bins = sum(Chirp_ksdensity_grid<cl_var.Chirp_trig_times_vec(2));
            
            Chirp_1_spike_density_mat  = Chirp_spike_density_mat(:,1:Num_Chirp_1_bins);
            Chirp_2_spike_density_mat  = Chirp_spike_density_mat(:,Num_Chirp_1_bins+1:end);
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [Chirp_coeff,Chirp_score,Chirp_latent,Chirp_tsquared,Chirp_explained,Chirp_mu] = pca(Chirp_spike_density_mat); %  Chirp_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Chirp_explained_cumsum  = cumsum(Chirp_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Chirp_PC_Num = find(Chirp_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Chirp_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,Chirp_score(:,1:Chirp_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Chirp_PC_Num;
            
            disp(sprintf('Num Chirp PCs = %i',Chirp_PC_Num));
            
            Chirp_pca_fig = figure;
            imagesc(Chirp_ksdensity_grid,[],abs(Chirp_coeff(:,1:Chirp_PC_Num)'));%Chirp_PC_Num, 10
            colormap(Chirp_pca_fig,gray(256));
            colorbar;
            Chirp_ks_v1 = vline(cl_var.Chirp_trig_times_vec(2),'r');
            title('Chirp - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(Chirp_ks_v1,'LineWidth',1.5);
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [Chirp_1_coeff,Chirp_1_score,Chirp_1_latent,Chirp_1_tsquared,Chirp_1_explained,Chirp_1_mu] = pca(Chirp_1_spike_density_mat);
            [Chirp_2_coeff,Chirp_2_score,Chirp_2_latent,Chirp_2_tsquared,Chirp_2_explained,Chirp_2_mu] = pca(Chirp_2_spike_density_mat);
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Chirp_1_explained_cumsum  = cumsum(Chirp_1_explained);
                Chirp_2_explained_cumsum  = cumsum(Chirp_2_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Chirp_1_PC_Num  = find(Chirp_1_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Chirp_2_PC_Num  = find(Chirp_2_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Chirp_1_PC_Num  = p.PCA_comp_num;
                Chirp_2_PC_Num  = p.PCA_comp_num;
                
            end
            
            Chirp_PC_Num = Chirp_1_PC_Num + Chirp_2_PC_Num;
            
            All_Scores     = [All_Scores,Chirp_1_score(:,1:Chirp_1_PC_Num),Chirp_2_score(:,1:Chirp_2_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Chirp_PC_Num;
            
            disp(sprintf('Num Chirp 1 PCs = %i',Chirp_1_PC_Num));
            disp(sprintf('Num Chirp 2 PCs = %i',Chirp_2_PC_Num));
            disp(sprintf('Total Num Chirp PCs = %i',Chirp_PC_Num));
            
            Chirp_pca_fig = figure;
            subplot(1,2,1);
            imagesc(Chirp_ksdensity_grid(1:Num_Chirp_1_bins),[],abs(Chirp_1_coeff(:,1:Chirp_1_PC_Num)'));
            colormap(Chirp_pca_fig,gray(256));
            colorbar;
            title('Chirp 1 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(1,2,2);
            imagesc(Chirp_ksdensity_grid(Num_Chirp_1_bins+1:end),[],abs(Chirp_2_coeff(:,1:Chirp_2_PC_Num)'));
            colormap(Chirp_pca_fig,gray(256));
            title('Chirp 2 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% FFF_Noise
if p.Obj_clust_vec(3) == 1 || p.Obj_plot_vec(5) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    % N/A here as already have kernels
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(3) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
%             FFF_Noise_STA_R  = cl_var.FFF_Noise_STA(:,:,1); % FFF_R_spike_density_mat
%             FFF_Noise_STA_G  = cl_var.FFF_Noise_STA(:,:,2);
%             FFF_Noise_STA_B  = cl_var.FFF_Noise_STA(:,:,3);
%             FFF_Noise_STA_UV = cl_var.FFF_Noise_STA(:,:,4);
            
            % Try Scaled Versions
            FFF_Noise_STA_R  = cl_var.FFF_Noise_STA_scaled(:,:,1); % 14,06,2021 Mod
            FFF_Noise_STA_G  = cl_var.FFF_Noise_STA_scaled(:,:,2); % 14,06,2021 Mod
            FFF_Noise_STA_B  = cl_var.FFF_Noise_STA_scaled(:,:,3); % 14,06,2021 Mod
            FFF_Noise_STA_UV = cl_var.FFF_Noise_STA_scaled(:,:,4); % 14,06,2021 Mod
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [FFF_Noise_coeff,FFF_Noise_score,FFF_Noise_latent,FFF_Noise_tsquared,FFF_Noise_explained,FFF_Noise_mu] = pca(cl_var.FFF_Noise_Full_STA);
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                FFF_Noise_explained_cumsum  = cumsum(FFF_Noise_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                FFF_Noise_PC_Num = find(FFF_Noise_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                FFF_Noise_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,FFF_Noise_score(:,1:FFF_Noise_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + FFF_Noise_PC_Num;
            
            disp(sprintf('Num FFF Noise PCs = %i',FFF_Noise_PC_Num));
            
            FFF_Noise_pca_fig = figure;
            imagesc(cl_var.FFF_Noise_Full_t_vec+0.5*cl_var.FFF_Noise_Full_t_vec_int,[],abs(FFF_Noise_coeff(:,1:FFF_Noise_PC_Num)'));
            colormap(FFF_Noise_pca_fig,gray(256));
            colorbar;
            FFF_Noise_STA_v1 = vline(0.25*(cl_var.FFF_Noise_Full_t_vec_end+cl_var.FFF_Noise_Full_t_vec_int),'r');
            FFF_Noise_STA_v2 = vline(0.5*(cl_var.FFF_Noise_Full_t_vec_end+cl_var.FFF_Noise_Full_t_vec_int),'r');
            FFF_Noise_STA_v3 = vline(0.75*(cl_var.FFF_Noise_Full_t_vec_end+cl_var.FFF_Noise_Full_t_vec_int),'r');
            title('FFF Noise - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(FFF_Noise_STA_v1,'LineWidth',1.5);
            set(FFF_Noise_STA_v2,'LineWidth',1.5);
            set(FFF_Noise_STA_v3,'LineWidth',1.5);
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [FFF_Noise_R_coeff,FFF_Noise_R_score,FFF_Noise_R_latent,FFF_Noise_R_tsquared,FFF_Noise_R_explained,FFF_Noise_R_mu]       = pca(FFF_Noise_STA_R);
            [FFF_Noise_G_coeff,FFF_Noise_G_score,FFF_Noise_G_latent,FFF_Noise_G_tsquared,FFF_Noise_G_explained,FFF_Noise_G_mu]       = pca(FFF_Noise_STA_G);
            [FFF_Noise_B_coeff,FFF_Noise_B_score,FFF_Noise_B_latent,FFF_Noise_B_tsquared,FFF_Noise_B_explained,FFF_Noise_B_mu]       = pca(FFF_Noise_STA_B);
            [FFF_Noise_UV_coeff,FFF_Noise_UV_score,FFF_Noise_UV_latent,FFF_Noise_UV_tsquared,FFF_Noise_UV_explained,FFF_Noise_UV_mu] = pca(FFF_Noise_STA_UV);
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                FFF_Noise_R_explained_cumsum  = cumsum(FFF_Noise_R_explained);
                FFF_Noise_G_explained_cumsum  = cumsum(FFF_Noise_G_explained);
                FFF_Noise_B_explained_cumsum  = cumsum(FFF_Noise_B_explained);
                FFF_Noise_UV_explained_cumsum = cumsum(FFF_Noise_UV_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                FFF_Noise_R_PC_Num  = find(FFF_Noise_R_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF_Noise_G_PC_Num  = find(FFF_Noise_G_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF_Noise_B_PC_Num  = find(FFF_Noise_B_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF_Noise_UV_PC_Num = find(FFF_Noise_UV_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                FFF_Noise_R_PC_Num  = p.PCA_comp_num;
                FFF_Noise_G_PC_Num  = p.PCA_comp_num;
                FFF_Noise_B_PC_Num  = p.PCA_comp_num;
                FFF_Noise_UV_PC_Num = p.PCA_comp_num;
                
            end
            
            FFF_Noise_PC_Num = FFF_Noise_R_PC_Num + FFF_Noise_G_PC_Num + FFF_Noise_B_PC_Num + FFF_Noise_UV_PC_Num;
            
            All_Scores     = [All_Scores,FFF_Noise_R_score(:,1:FFF_Noise_R_PC_Num),FFF_Noise_G_score(:,1:FFF_Noise_G_PC_Num),...
                FFF_Noise_B_score(:,1:FFF_Noise_B_PC_Num),FFF_Noise_UV_score(:,1:FFF_Noise_UV_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + FFF_Noise_PC_Num;
            
            disp(sprintf('Num FFF Noise R PCs = %i',FFF_Noise_R_PC_Num));
            disp(sprintf('Num FFF Noise G PCs = %i',FFF_Noise_G_PC_Num));
            disp(sprintf('Num FFF Noise B PCs = %i',FFF_Noise_B_PC_Num));
            disp(sprintf('Num FFF Noise UV PCs = %i',FFF_Noise_UV_PC_Num));
            disp(sprintf('Total Num FFF Noise PCs = %i',FFF_Noise_PC_Num));
            
            FFF_Noise_pca_fig = figure;
            subplot(2,2,1);
            imagesc(p.stim_timesample_vec,[],abs(FFF_Noise_R_coeff(:,1:FFF_Noise_R_PC_Num)'));
            colormap(FFF_Noise_pca_fig,gray(256));
            colorbar;
            title('FFF Noise R - PCA coeff');ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,2,2);
            imagesc(p.stim_timesample_vec,[],abs(FFF_Noise_G_coeff(:,1:FFF_Noise_G_PC_Num)'));
            colormap(FFF_Noise_pca_fig,gray(256));
            colorbar;
            title('FFF Noise G - PCA coeff');
            set(gca,'FontSize',12);
            subplot(2,2,3);
            imagesc(p.stim_timesample_vec,[],abs(FFF_Noise_B_coeff(:,1:FFF_Noise_B_PC_Num)'));
            colormap(FFF_Noise_pca_fig,gray(256));
            colorbar;
            title('FFF Noise B - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,2,4);
            imagesc(p.stim_timesample_vec,[],abs(FFF_Noise_UV_coeff(:,1:FFF_Noise_UV_PC_Num)'));
            colormap(FFF_Noise_pca_fig,gray(256));
            colorbar;
            title('FFF Noise UV - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% Gratings_400px
if p.Obj_clust_vec(4) == 1 || p.Obj_plot_vec(6) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    Gratings_400px_length_ksdensity_grid = 1024; % now: 1e3, was: 104 % 104 choose a number divisible by 8
    Gratings_400px_ksdensity_grid        = linspace(0,cl_var.Gratings_400px_stim_end_time,Gratings_400px_length_ksdensity_grid);
    Gratings_400px_ksdensity_bdwth       = 1e-2; % Tom said 1e-2 is best (of options presented)
    Gratings_400px_spike_density_mat     = NaN(cl_var.True_Num_Cells,Gratings_400px_length_ksdensity_grid);
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.Gratings_400px_spike_times_mat(:,i);cl_var.Gratings_400px_spike_times_mat(:,i);(2*cl_var.Gratings_400px_stim_end_time-cl_var.Gratings_400px_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            [Gratings_400px_spike_density_mat(i,:),~] = ksdensity(vec_loop,Gratings_400px_ksdensity_grid,'Bandwidth',Gratings_400px_ksdensity_bdwth); % ksdensity(cl_var.Gratings_400px_spike_times_mat(:,i),Gratings_400px_ksdensity_grid,'Bandwidth',Gratings_400px_ksdensity_bdwth);
        else % If there were no spikes for this cell for this stim.
            Gratings_400px_spike_density_mat(i,:) = zeros(1,Gratings_400px_length_ksdensity_grid);
        end
    end
    
    % 1D plot of chosen cell Gratings_400px spike rate
    figure;
    plot(Gratings_400px_ksdensity_grid,Gratings_400px_spike_density_mat(2,:),'LineWidth',1.5); hold on;
    plot(cl_var.Gratings_400px_spike_times_mat(:,2),zeros(size(cl_var.Gratings_400px_spike_times_mat,1),1),'bx');
    Gratings_400px_ks_v1 = vline(cl_var.Gratings_400px_trig_times_vec(2),'r');
    Gratings_400px_ks_v2 = vline(cl_var.Gratings_400px_trig_times_vec(3),'r');
    Gratings_400px_ks_v3 = vline(cl_var.Gratings_400px_trig_times_vec(4),'r');
    Gratings_400px_ks_v4 = vline(cl_var.Gratings_400px_trig_times_vec(5),'r');
    Gratings_400px_ks_v5 = vline(cl_var.Gratings_400px_trig_times_vec(6),'r');
    Gratings_400px_ks_v6 = vline(cl_var.Gratings_400px_trig_times_vec(7),'r');
    Gratings_400px_ks_v7 = vline(cl_var.Gratings_400px_trig_times_vec(8),'r');
    xlim([0 cl_var.Gratings_400px_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('Gratings 400px spike rate');
    set(Gratings_400px_ks_v1,'LineWidth',1.5);
    set(Gratings_400px_ks_v2,'LineWidth',1.5);
    set(Gratings_400px_ks_v3,'LineWidth',1.5);
    set(Gratings_400px_ks_v4,'LineWidth',1.5);
    set(Gratings_400px_ks_v5,'LineWidth',1.5);
    set(Gratings_400px_ks_v6,'LineWidth',1.5);
    set(Gratings_400px_ks_v7,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell Gratings_400px spike rates
    Gratings_400px_ks_fig = figure;
    imagesc(Gratings_400px_ksdensity_grid,[],Gratings_400px_spike_density_mat); hold on;
    colormap(Gratings_400px_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    Gratings_400px_ks_v1 = vline(cl_var.Gratings_400px_trig_times_vec(2),'r');
    Gratings_400px_ks_v2 = vline(cl_var.Gratings_400px_trig_times_vec(3),'r');
    Gratings_400px_ks_v3 = vline(cl_var.Gratings_400px_trig_times_vec(4),'r');
    Gratings_400px_ks_v4 = vline(cl_var.Gratings_400px_trig_times_vec(5),'r');
    Gratings_400px_ks_v5 = vline(cl_var.Gratings_400px_trig_times_vec(6),'r');
    Gratings_400px_ks_v6 = vline(cl_var.Gratings_400px_trig_times_vec(7),'r');
    Gratings_400px_ks_v7 = vline(cl_var.Gratings_400px_trig_times_vec(8),'r');
    xlabel('time (sec)');
    ylabel('RF');
    title('Gratings 400px spike rates');
    set(Gratings_400px_ks_v1,'LineWidth',1.5);
    set(Gratings_400px_ks_v2,'LineWidth',1.5);
    set(Gratings_400px_ks_v3,'LineWidth',1.5);
    set(Gratings_400px_ks_v4,'LineWidth',1.5);
    set(Gratings_400px_ks_v5,'LineWidth',1.5);
    set(Gratings_400px_ks_v6,'LineWidth',1.5);
    set(Gratings_400px_ks_v7,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    %%% 4.5. Summary Statistic Analysis
    
    % SVG Analysis
    
    %alpha   = [pi;3*pi;5*pi;7*pi;9*pi;11*pi;13*pi;15*pi]/8;
    alpha   = [0;pi;pi/4;5*pi/4;pi/2;3*pi/2;3*pi/4;7*pi/4];
    phi     = exp(complex(0,1)*alpha);
    % Check
%     figure;
%     for i=1:8
%         polarplot(phi(i),'o'); hold on;
%     end
%     legend;
    
    % True direction selectivity
    
    K_vec   = NaN(cl_var.True_Num_Cells,1);
    DSi_vec = NaN(cl_var.True_Num_Cells,1);
    
    for i = 1:cl_var.True_Num_Cells
        
        A_loop = [Gratings_400px_1_spike_density_mat(i,:)',...
            Gratings_400px_2_spike_density_mat(i,:)',...
            Gratings_400px_3_spike_density_mat(i,:)',...
            Gratings_400px_4_spike_density_mat(i,:)',...
            Gratings_400px_5_spike_density_mat(i,:)',...
            Gratings_400px_6_spike_density_mat(i,:)',...
            Gratings_400px_7_spike_density_mat(i,:)',...
            Gratings_400px_8_spike_density_mat(i,:)'];
        
        [U_loop,S_loop,V_loop] = svd(A_loop);
        
        V1_loop    = V_loop(:,1);       % V1 is negative for cell 14
        K_vec(i)   = dot(phi',V1_loop); % Dot product with complex numbers u and v is sum conj(u)*v, so take conj of phi to keep it as u*v.
        DSi_vec(i) = abs(K_vec(i));
        
    end
    
    figure;
    for i = 1 : cl_var.True_Num_Cells
        polarplot([0,K_vec(i)],'b-'); hold on;
    end
    polarplot(K_vec,'ro');
    title('direction selectivity');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % Bootstrapping check
    
    Gratings_400px_Num_BS_iter  = 1000; % Number of bootstrapping iterations
    Gratings_400px_BS_Amp       = cl_var.Gratings_400px_trig_times_vec(2); % BS noise amplitude (stimulus (direction) duration)
    
    Gratings_400px_BS_Noise_vec = 2*Gratings_400px_BS_Amp*(rand(Gratings_400px_Num_BS_iter,1) - 0.5);
    
    K_BS_mat   = NaN(cl_var.True_Num_Cells,Gratings_400px_Num_BS_iter); % Num cells x Num BS iter
    DSi_BS_mat = NaN(cl_var.True_Num_Cells,Gratings_400px_Num_BS_iter); % Num cells x Num BS iter
    
    for i = 1:cl_var.True_Num_Cells
        
        if any(~isnan(cl_var.Gratings_400px_spike_times_mat(:,i)))
            
            Gratings_400px_spike_times_vec_loop = cl_var.Gratings_400px_spike_times_mat(:,i);
            Gratings_400px_spike_times_vec_loop(isnan(Gratings_400px_spike_times_vec_loop)) = [];
            
            for j = 1:Gratings_400px_Num_BS_iter
                
                Gratings_400px_spike_times_vec_BS_loop    = mod(Gratings_400px_spike_times_vec_loop + Gratings_400px_BS_Noise_vec(j),cl_var.Gratings_400px_stim_end_time);
                vec_loop                                  = [-Gratings_400px_spike_times_vec_BS_loop;Gratings_400px_spike_times_vec_BS_loop;(2*cl_var.Gratings_400px_stim_end_time-Gratings_400px_spike_times_vec_BS_loop)];
                [Gratings_400px_spike_density_vec_loop,~] = ksdensity(vec_loop,Gratings_400px_ksdensity_grid,'Bandwidth',Gratings_400px_ksdensity_bdwth);
                
                A_loop = [Gratings_400px_spike_density_vec_loop(1:0.125*Gratings_400px_length_ksdensity_grid)',...
                          Gratings_400px_spike_density_vec_loop(0.125*Gratings_400px_length_ksdensity_grid+1:0.25*Gratings_400px_length_ksdensity_grid)',...
                          Gratings_400px_spike_density_vec_loop(0.25*Gratings_400px_length_ksdensity_grid+1:0.375*Gratings_400px_length_ksdensity_grid)',...
                          Gratings_400px_spike_density_vec_loop(0.375*Gratings_400px_length_ksdensity_grid+1:0.5*Gratings_400px_length_ksdensity_grid)',...
                          Gratings_400px_spike_density_vec_loop(0.5*Gratings_400px_length_ksdensity_grid+1:0.625*Gratings_400px_length_ksdensity_grid)',...
                          Gratings_400px_spike_density_vec_loop(0.625*Gratings_400px_length_ksdensity_grid+1:0.75*Gratings_400px_length_ksdensity_grid)',...
                          Gratings_400px_spike_density_vec_loop(0.75*Gratings_400px_length_ksdensity_grid+1:0.875*Gratings_400px_length_ksdensity_grid)',...
                          Gratings_400px_spike_density_vec_loop(0.875*Gratings_400px_length_ksdensity_grid+1:Gratings_400px_length_ksdensity_grid)'];
                
                      [U_loop,S_loop,V_loop] = svd(A_loop);
                      
                      V1_loop         = V_loop(:,1);
                      K_BS_mat(i,j)   = dot(phi',V1_loop); % Dot product with complex numbers u and v is sum conj(u)*v, so take conj of phi to keep it as u*v.
                      DSi_BS_mat(i,j) = abs(K_BS_mat(i,j));
                      
%                       figure;
%                       subplot(1,2,1);
%                       imagesc(Gratings_400px_ksdensity_grid(1:0.125*Gratings_400px_length_ksdensity_grid),alpha([1,3,5,7,2,4,6,8]),A_loop(:,[1,3,5,7,2,4,6,8])');
%                       set(gca,'YDir','normal');
%                       xlabel('time (s)');
%                       ylabel('direction (rad)');
%                       title('mean response matrix'); % not normalised but each direction equally weighted
%                       colorbar;
%                       set(gca,'FontSize',12);
%                       subplot(1,2,2);
%                       A_approx_loop = S_loop(1,1)*U_loop(:,1)*V_loop(:,1)';
%                       imagesc(Gratings_400px_ksdensity_grid(1:0.125*Gratings_400px_length_ksdensity_grid),alpha([1,3,5,7,2,4,6,8]),A_approx_loop(:,[1,3,5,7,2,4,6,8])');
%                       set(gca,'YDir','normal');
%                       xlabel('time (s)');
%                       %ylabel('direction (rad)');
%                       title('SVD reconstructed mean response matrix'); % not normalised but each direction equally weighted
%                       colorbar;
%                       set(gca,'FontSize',12);
%                       set(gcf,'color','w');
                      
    
            end
            
        else % If there were no spikes for this cell for this stim.
            
            %Gratings_400px_spike_density_vec_loop =
            %zeros(1,Gratings_400px_length_ksdensity_grid); % Not required.
            %K_BS_mat(i,j)   = NaN; % Not required as NaN already.
            %DSi_BS_mat(i,j) = NaN;% Not required as NaN already.
            
        end
        
    end
    
    figure;
    histogram(DSi_BS_mat(14,:)); hold on;
    Gratings_400px_ks_v1 = vline(DSi_vec(14),'red');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    set(Gratings_400px_ks_v1,'LineWidth',1.5);
    xlabel('DSi');
    ylabel('num. values.');
    title('direction selectivity boostrapping');
    
    % Fast Fourier Transform
    
%     Gratings_400px_SD_thresh_fac = 2;
%     
%     Gratings_400px_1_spike_density_mat = Gratings_400px_spike_density_mat(:,1:0.125*Gratings_400px_length_ksdensity_grid);
%     Gratings_400px_2_spike_density_mat = Gratings_400px_spike_density_mat(:,0.125*Gratings_400px_length_ksdensity_grid+1:0.25*Gratings_400px_length_ksdensity_grid);
%     Gratings_400px_3_spike_density_mat = Gratings_400px_spike_density_mat(:,0.25*Gratings_400px_length_ksdensity_grid+1:0.375*Gratings_400px_length_ksdensity_grid);
%     Gratings_400px_4_spike_density_mat = Gratings_400px_spike_density_mat(:,0.375*Gratings_400px_length_ksdensity_grid+1:0.5*Gratings_400px_length_ksdensity_grid);
%     Gratings_400px_5_spike_density_mat = Gratings_400px_spike_density_mat(:,0.5*Gratings_400px_length_ksdensity_grid+1:0.625*Gratings_400px_length_ksdensity_grid);
%     Gratings_400px_6_spike_density_mat = Gratings_400px_spike_density_mat(:,0.625*Gratings_400px_length_ksdensity_grid+1:0.75*Gratings_400px_length_ksdensity_grid);
%     Gratings_400px_7_spike_density_mat = Gratings_400px_spike_density_mat(:,0.75*Gratings_400px_length_ksdensity_grid+1:0.875*Gratings_400px_length_ksdensity_grid);
%     Gratings_400px_8_spike_density_mat = Gratings_400px_spike_density_mat(:,0.875*Gratings_400px_length_ksdensity_grid+1:Gratings_400px_length_ksdensity_grid);
%     
%     Gratings_400px_T  = cl_var.Gratings_400px_stim_end_time/Gratings_400px_length_ksdensity_grid; % Sampling period % cl_var.Gratings_400px_stim_end_time/8
%     Gratings_400px_Fs = 1/Gratings_400px_T;                                                       % Sampling frequency
%     Gratings_400px_L  = Gratings_400px_length_ksdensity_grid/8;                                   % Length of signal
%     
%     Gratings_400px_freq_vec = 0:(Gratings_400px_Fs/Gratings_400px_n):(Gratings_400px_Fs/2-Gratings_400px_Fs/Gratings_400px_n);
%     
%     Gratings_400px_n = 2^nextpow2(Gratings_400px_L+1); % Add one to get 256 rather than same value of 128.
%     
%     Gratings_400px_1_fft = fft(Gratings_400px_1_spike_density_mat,Gratings_400px_n,2);
%     Gratings_400px_2_fft = fft(Gratings_400px_2_spike_density_mat,Gratings_400px_n,2);
%     Gratings_400px_3_fft = fft(Gratings_400px_3_spike_density_mat,Gratings_400px_n,2);
%     Gratings_400px_4_fft = fft(Gratings_400px_4_spike_density_mat,Gratings_400px_n,2);
%     Gratings_400px_5_fft = fft(Gratings_400px_5_spike_density_mat,Gratings_400px_n,2);
%     Gratings_400px_6_fft = fft(Gratings_400px_6_spike_density_mat,Gratings_400px_n,2);
%     Gratings_400px_7_fft = fft(Gratings_400px_7_spike_density_mat,Gratings_400px_n,2);
%     Gratings_400px_8_fft = fft(Gratings_400px_8_spike_density_mat,Gratings_400px_n,2);
%     % If X is a matrix, then fft(X) treats the columns of X as vectors and returns the Fourier transform of each column.
%     % Y = fft(X,n) returns the n-point DFT. If no value is specified, Y is the same size as X.
%     % Y = fft(X,n,dim) returns the Fourier transform along the dimension dim. For example, if X is a matrix, then fft(X,n,2) returns the n-point Fourier transform of each row.
%     
%     Gratings_400px_1_P2            = abs(Gratings_400px_1_fft/Gratings_400px_L);
%     Gratings_400px_1_P1            = Gratings_400px_1_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_1_P1(:,2:end-1) = 2*Gratings_400px_1_P1(:,2:end-1);
%     Gratings_400px_1_P1_mean       = mean(Gratings_400px_1_P1,'all');
%     Gratings_400px_1_P1_std        = std(Gratings_400px_1_P1,0,'all');
%     
%     Gratings_400px_2_P2            = abs(Gratings_400px_2_fft/Gratings_400px_L);
%     Gratings_400px_2_P1            = Gratings_400px_2_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_2_P1(:,2:end-1) = 2*Gratings_400px_2_P1(:,2:end-1);
%     
%     Gratings_400px_3_P2            = abs(Gratings_400px_3_fft/Gratings_400px_L);
%     Gratings_400px_3_P1            = Gratings_400px_3_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_3_P1(:,2:end-1) = 2*Gratings_400px_3_P1(:,2:end-1);
%     
%     Gratings_400px_4_P2            = abs(Gratings_400px_4_fft/Gratings_400px_L);
%     Gratings_400px_4_P1            = Gratings_400px_4_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_4_P1(:,2:end-1) = 2*Gratings_400px_4_P1(:,2:end-1);
%     
%     Gratings_400px_5_P2            = abs(Gratings_400px_5_fft/Gratings_400px_L);
%     Gratings_400px_5_P1            = Gratings_400px_5_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_5_P1(:,2:end-1) = 2*Gratings_400px_5_P1(:,2:end-1);
%     
%     Gratings_400px_6_P2            = abs(Gratings_400px_6_fft/Gratings_400px_L);
%     Gratings_400px_6_P1            = Gratings_400px_6_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_6_P1(:,2:end-1) = 2*Gratings_400px_6_P1(:,2:end-1);
%     
%     Gratings_400px_7_P2            = abs(Gratings_400px_7_fft/Gratings_400px_L);
%     Gratings_400px_7_P1            = Gratings_400px_7_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_7_P1(:,2:end-1) = 2*Gratings_400px_7_P1(:,2:end-1);
%     
%     Gratings_400px_8_P2            = abs(Gratings_400px_8_fft/Gratings_400px_L);
%     Gratings_400px_8_P1            = Gratings_400px_8_P2(:,1:Gratings_400px_n/2+1);
%     Gratings_400px_8_P1(:,2:end-1) = 2*Gratings_400px_8_P1(:,2:end-1);
%     
%     figure;
%     cell_choice_temp = 14;
%     subplot(4,2,1);
%     plot(Gratings_400px_freq_vec,Gratings_400px_1_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     ylabel('P1(fft)');
%     set(gca,'FontSize',12);
%     subplot(4,2,2);
%     plot(Gratings_400px_freq_vec,Gratings_400px_2_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     set(gca,'FontSize',12);
%     subplot(4,2,3);
%     plot(Gratings_400px_freq_vec,Gratings_400px_3_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     ylabel('P1(fft)');
%     set(gca,'FontSize',12);
%     subplot(4,2,4);
%     plot(Gratings_400px_freq_vec,Gratings_400px_4_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     set(gca,'FontSize',12);
%     subplot(4,2,5);
%     plot(Gratings_400px_freq_vec,Gratings_400px_5_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     ylabel('P1(fft)');
%     set(gca,'FontSize',12);
%     subplot(4,2,6);
%     plot(Gratings_400px_freq_vec,Gratings_400px_6_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     set(gca,'FontSize',12);
%     subplot(4,2,7);
%     plot(Gratings_400px_freq_vec,Gratings_400px_7_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     xlabel('freq. (1/sec)');
%     ylabel('P1(fft)');
%     set(gca,'FontSize',12);
%     subplot(4,2,8);
%     plot(Gratings_400px_freq_vec,Gratings_400px_8_P1(cell_choice_temp,1:Gratings_400px_n/2),'LineWidth',1.5); hold on;
%     Gratings_400px_h1 = hline(Gratings_400px_1_P1_mean);
%     Gratings_400px_h2 = hline(Gratings_400px_1_P1_mean + Gratings_400px_SD_thresh_fac*Gratings_400px_1_P1_std);
%     set(Gratings_400px_h1,'LineWidth',1.5,'LineStyle','-','Color','r');
%     set(Gratings_400px_h2,'LineWidth',1.5,'LineStyle','-','Color','g');
%     xlim([0 Gratings_400px_freq_vec(end)]);
%     xlabel('freq. (1/sec)');
%     set(gca,'FontSize',12);
%     set(gcf,'color','w');
%     
%     figure;
%     subplot(4,2,1);
%     imagesc(Gratings_400px_1_P1(:,1:Gratings_400px_n/2));
%     ylabel('cell');
%     set(gca,'FontSize',12);
%     subplot(4,2,2);
%     imagesc(Gratings_400px_2_P1(:,1:Gratings_400px_n/2));
%     set(gca,'FontSize',12);
%     subplot(4,2,3);
%     imagesc(Gratings_400px_3_P1(:,1:Gratings_400px_n/2));
%     ylabel('cell');
%     set(gca,'FontSize',12);
%     subplot(4,2,4);
%     imagesc(Gratings_400px_4_P1(:,1:Gratings_400px_n/2));
%     set(gca,'FontSize',12);
%     subplot(4,2,5);
%     imagesc(Gratings_400px_5_P1(:,1:Gratings_400px_n/2));
%     ylabel('cell');
%     set(gca,'FontSize',12);
%     subplot(4,2,6);
%     imagesc(Gratings_400px_6_P1(:,1:Gratings_400px_n/2));
%     set(gca,'FontSize',12);
%     subplot(4,2,7);
%     imagesc(Gratings_400px_7_P1(:,1:Gratings_400px_n/2));
%     ylabel('cell');
%     xlabel('freq. (1/sec)');
%     set(gca,'FontSize',12);
%     subplot(4,2,8);
%     imagesc(Gratings_400px_8_P1(:,1:Gratings_400px_n/2));
%     xlabel('freq. (1/sec)');
%     set(gca,'FontSize',12);
%     set(gcf,'color','w');
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(4) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            %Gratings_400px_segment_time = 0.25*cl_var.Gratings_400px_stim_end_time; % Could have used this...
            
            Gratings_400px_1_spike_density_mat = Gratings_400px_spike_density_mat(:,1:0.125*Gratings_400px_length_ksdensity_grid);
            Gratings_400px_2_spike_density_mat = Gratings_400px_spike_density_mat(:,0.125*Gratings_400px_length_ksdensity_grid+1:0.25*Gratings_400px_length_ksdensity_grid);
            Gratings_400px_3_spike_density_mat = Gratings_400px_spike_density_mat(:,0.25*Gratings_400px_length_ksdensity_grid+1:0.375*Gratings_400px_length_ksdensity_grid);
            Gratings_400px_4_spike_density_mat = Gratings_400px_spike_density_mat(:,0.375*Gratings_400px_length_ksdensity_grid+1:0.5*Gratings_400px_length_ksdensity_grid);
            Gratings_400px_5_spike_density_mat = Gratings_400px_spike_density_mat(:,0.5*Gratings_400px_length_ksdensity_grid+1:0.625*Gratings_400px_length_ksdensity_grid);
            Gratings_400px_6_spike_density_mat = Gratings_400px_spike_density_mat(:,0.625*Gratings_400px_length_ksdensity_grid+1:0.75*Gratings_400px_length_ksdensity_grid);
            Gratings_400px_7_spike_density_mat = Gratings_400px_spike_density_mat(:,0.75*Gratings_400px_length_ksdensity_grid+1:0.875*Gratings_400px_length_ksdensity_grid);
            Gratings_400px_8_spike_density_mat = Gratings_400px_spike_density_mat(:,0.875*Gratings_400px_length_ksdensity_grid+1:Gratings_400px_length_ksdensity_grid);
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [Gratings_400px_coeff,Gratings_400px_score,Gratings_400px_latent,Gratings_400px_tsquared,Gratings_400px_explained,Gratings_400px_mu] = pca(Gratings_400px_spike_density_mat); %  Gratings_400px_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Gratings_400px_explained_cumsum = cumsum(Gratings_400px_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Gratings_400px_PC_Num = find(Gratings_400px_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Gratings_400px_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,Gratings_400px_score(:,1:Gratings_400px_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Gratings_400px_PC_Num;
            
            disp(sprintf('Num Gratings_400px PCs = %i',Gratings_400px_PC_Num));
            
            Gratings_400px_pca_fig = figure;
            imagesc(Gratings_400px_ksdensity_grid,[],abs(Gratings_400px_coeff(:,1:Gratings_400px_PC_Num)'));%Gratings_400px_PC_Num, 10
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            Gratings_400px_ks_v1 = vline(cl_var.Gratings_400px_trig_times_vec(2),'r');
            Gratings_400px_ks_v2 = vline(cl_var.Gratings_400px_trig_times_vec(3),'r');
            Gratings_400px_ks_v3 = vline(cl_var.Gratings_400px_trig_times_vec(4),'r');
            Gratings_400px_ks_v4 = vline(cl_var.Gratings_400px_trig_times_vec(5),'r');
            Gratings_400px_ks_v5 = vline(cl_var.Gratings_400px_trig_times_vec(6),'r');
            Gratings_400px_ks_v6 = vline(cl_var.Gratings_400px_trig_times_vec(7),'r');
            Gratings_400px_ks_v7 = vline(cl_var.Gratings_400px_trig_times_vec(8),'r');
            title('Gratings 400px - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(Gratings_400px_ks_v1,'LineWidth',1.5);
            set(Gratings_400px_ks_v2,'LineWidth',1.5);
            set(Gratings_400px_ks_v3,'LineWidth',1.5);
            set(Gratings_400px_ks_v4,'LineWidth',1.5);
            set(Gratings_400px_ks_v5,'LineWidth',1.5);
            set(Gratings_400px_ks_v6,'LineWidth',1.5);
            set(Gratings_400px_ks_v7,'LineWidth',1.5);
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [Gratings_400px_1_coeff,Gratings_400px_1_score,Gratings_400px_1_latent,Gratings_400px_1_tsquared,Gratings_400px_1_explained,Gratings_400px_1_mu] = pca(Gratings_400px_1_spike_density_mat);
            [Gratings_400px_2_coeff,Gratings_400px_2_score,Gratings_400px_2_latent,Gratings_400px_2_tsquared,Gratings_400px_2_explained,Gratings_400px_2_mu] = pca(Gratings_400px_2_spike_density_mat);
            [Gratings_400px_3_coeff,Gratings_400px_3_score,Gratings_400px_3_latent,Gratings_400px_3_tsquared,Gratings_400px_3_explained,Gratings_400px_3_mu] = pca(Gratings_400px_3_spike_density_mat);
            [Gratings_400px_4_coeff,Gratings_400px_4_score,Gratings_400px_4_latent,Gratings_400px_4_tsquared,Gratings_400px_4_explained,Gratings_400px_4_mu] = pca(Gratings_400px_4_spike_density_mat);
            [Gratings_400px_5_coeff,Gratings_400px_5_score,Gratings_400px_5_latent,Gratings_400px_5_tsquared,Gratings_400px_5_explained,Gratings_400px_5_mu] = pca(Gratings_400px_5_spike_density_mat);
            [Gratings_400px_6_coeff,Gratings_400px_6_score,Gratings_400px_6_latent,Gratings_400px_6_tsquared,Gratings_400px_6_explained,Gratings_400px_6_mu] = pca(Gratings_400px_6_spike_density_mat);
            [Gratings_400px_7_coeff,Gratings_400px_7_score,Gratings_400px_7_latent,Gratings_400px_7_tsquared,Gratings_400px_7_explained,Gratings_400px_7_mu] = pca(Gratings_400px_7_spike_density_mat);
            [Gratings_400px_8_coeff,Gratings_400px_8_score,Gratings_400px_8_latent,Gratings_400px_8_tsquared,Gratings_400px_8_explained,Gratings_400px_8_mu] = pca(Gratings_400px_8_spike_density_mat);
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Gratings_400px_1_explained_cumsum = cumsum(Gratings_400px_1_explained);
                Gratings_400px_2_explained_cumsum = cumsum(Gratings_400px_2_explained);
                Gratings_400px_3_explained_cumsum = cumsum(Gratings_400px_3_explained);
                Gratings_400px_4_explained_cumsum = cumsum(Gratings_400px_4_explained);
                Gratings_400px_5_explained_cumsum = cumsum(Gratings_400px_5_explained);
                Gratings_400px_6_explained_cumsum = cumsum(Gratings_400px_6_explained);
                Gratings_400px_7_explained_cumsum = cumsum(Gratings_400px_7_explained);
                Gratings_400px_8_explained_cumsum = cumsum(Gratings_400px_8_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Gratings_400px_1_PC_Num = find(Gratings_400px_1_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Gratings_400px_2_PC_Num = find(Gratings_400px_2_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Gratings_400px_3_PC_Num = find(Gratings_400px_3_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Gratings_400px_4_PC_Num = find(Gratings_400px_4_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Gratings_400px_5_PC_Num = find(Gratings_400px_5_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Gratings_400px_6_PC_Num = find(Gratings_400px_6_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Gratings_400px_7_PC_Num = find(Gratings_400px_7_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Gratings_400px_8_PC_Num = find(Gratings_400px_8_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Gratings_400px_1_PC_Num = p.PCA_comp_num;
                Gratings_400px_2_PC_Num = p.PCA_comp_num;
                Gratings_400px_3_PC_Num = p.PCA_comp_num;
                Gratings_400px_4_PC_Num = p.PCA_comp_num;
                Gratings_400px_5_PC_Num = p.PCA_comp_num;
                Gratings_400px_6_PC_Num = p.PCA_comp_num;
                Gratings_400px_7_PC_Num = p.PCA_comp_num;
                Gratings_400px_8_PC_Num = p.PCA_comp_num;
                
            end
            
            Gratings_400px_PC_Num = Gratings_400px_1_PC_Num + Gratings_400px_2_PC_Num + Gratings_400px_3_PC_Num + Gratings_400px_4_PC_Num +...
                Gratings_400px_5_PC_Num + Gratings_400px_6_PC_Num + Gratings_400px_7_PC_Num + Gratings_400px_8_PC_Num;
            
            All_Scores     = [All_Scores,Gratings_400px_1_score(:,1:Gratings_400px_1_PC_Num),Gratings_400px_2_score(:,1:Gratings_400px_2_PC_Num),...
                Gratings_400px_3_score(:,1:Gratings_400px_3_PC_Num),Gratings_400px_4_score(:,1:Gratings_400px_4_PC_Num),...
                Gratings_400px_5_score(:,1:Gratings_400px_5_PC_Num),Gratings_400px_6_score(:,1:Gratings_400px_6_PC_Num),...
                Gratings_400px_7_score(:,1:Gratings_400px_7_PC_Num),Gratings_400px_8_score(:,1:Gratings_400px_8_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Gratings_400px_PC_Num;
            
            disp(sprintf('Num Gratings 400px 1 PCs = %i',Gratings_400px_1_PC_Num));
            disp(sprintf('Num Gratings 400px 2 PCs = %i',Gratings_400px_2_PC_Num));
            disp(sprintf('Num Gratings 400px 3 PCs = %i',Gratings_400px_3_PC_Num));
            disp(sprintf('Num Gratings 400px 4 PCs = %i',Gratings_400px_4_PC_Num));
            disp(sprintf('Num Gratings 400px 5 PCs = %i',Gratings_400px_5_PC_Num));
            disp(sprintf('Num Gratings 400px 6 PCs = %i',Gratings_400px_6_PC_Num));
            disp(sprintf('Num Gratings 400px 7 PCs = %i',Gratings_400px_7_PC_Num));
            disp(sprintf('Num Gratings 400px 8 PCs = %i',Gratings_400px_8_PC_Num));
            disp(sprintf('Total Num Gratings 400px PCs = %i',Gratings_400px_PC_Num));
            
            Gratings_400px_pca_fig = figure;
            subplot(2,4,1);
            imagesc(Gratings_400px_ksdensity_grid(1:0.125*Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_1_coeff(:,1:Gratings_400px_1_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 1 - PCA coeff'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,4,2);
            imagesc(Gratings_400px_ksdensity_grid(0.125*Gratings_400px_length_ksdensity_grid+1:0.25*Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_2_coeff(:,1:Gratings_400px_2_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 2 - PCA coeff');
            set(gca,'FontSize',12);
            subplot(2,4,3);
            imagesc(Gratings_400px_ksdensity_grid(0.25*Gratings_400px_length_ksdensity_grid+1:0.375*Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_3_coeff(:,1:Gratings_400px_3_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 3 - PCA coeff');
            set(gca,'FontSize',12);
            subplot(2,4,4);
            imagesc(Gratings_400px_ksdensity_grid(0.375*Gratings_400px_length_ksdensity_grid+1:0.5*Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_4_coeff(:,1:Gratings_400px_4_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 4 - PCA coeff');
            set(gca,'FontSize',12);
            subplot(2,4,5);
            imagesc(Gratings_400px_ksdensity_grid(0.5*Gratings_400px_length_ksdensity_grid+1:0.625*Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_5_coeff(:,1:Gratings_400px_5_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 5 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,4,6);
            imagesc(Gratings_400px_ksdensity_grid(0.625*Gratings_400px_length_ksdensity_grid+1:0.75*Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_6_coeff(:,1:Gratings_400px_6_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 6 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            subplot(2,4,7);
            imagesc(Gratings_400px_ksdensity_grid(0.75*Gratings_400px_length_ksdensity_grid+1:0.875*Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_7_coeff(:,1:Gratings_400px_7_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 7 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            subplot(2,4,8);
            imagesc(Gratings_400px_ksdensity_grid(0.875*Gratings_400px_length_ksdensity_grid+1:Gratings_400px_length_ksdensity_grid),[],abs(Gratings_400px_8_coeff(:,1:Gratings_400px_8_PC_Num)'));
            colormap(Gratings_400px_pca_fig,gray(256));
            colorbar;
            title('Gratings 400px 8 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% CNoise: Full RF Size
if p.Obj_clust_vec(5) == 1% || p.Obj_plot_vec(7) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
    
    % Scale all scores later
%     Full_RF_Size_vec_min    = min(cl_var.Full_RF_Size_vec);
%     Full_RF_Size_vec_max    = max(cl_var.Full_RF_Size_vec);
%     Full_RF_Size_vec_scaled = (cl_var.Full_RF_Size_vec-Full_RF_Size_vec_min)/(Full_RF_Size_vec_max - Full_RF_Size_vec_min);
%     All_Scores              = [All_Scores,Full_RF_Size_vec_scaled];
    All_Scores              = [All_Scores,cl_var.Full_RF_Size_vec];

end



%%% CNoise: Full RF Ellipticity
if p.Obj_clust_vec(6) == 1% || p.Obj_plot_vec(8) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
    
    % Scale all scores later
%     Full_RF_Ellip_vec_min    = min(cl_var.Full_RF_Ellipticity_vec);
%     Full_RF_Ellip_vec_max    = max(cl_var.Full_RF_Ellipticity_vec);
%     Full_RF_Ellip_vec_scaled = (cl_var.Full_RF_Ellipticity_vec-Full_RF_Ellip_vec_min)/(Full_RF_Ellip_vec_max - Full_RF_Ellip_vec_min);
%     All_Scores               = [All_Scores,Full_RF_Ellip_vec_scaled];
    All_Scores               = [All_Scores,cl_var.Full_RF_Ellipticity_vec];
    
end



%%% CNoise: Full RF Dominant Axis Angle
if p.Obj_clust_vec(7) == 1% || p.Obj_plot_vec(9) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
    
    All_Scores              = [All_Scores,cl_var.Full_RF_Dom_Ax_Ang_vec];

end



%%% FFF2
if p.Obj_clust_vec(8) == 1 || p.Obj_plot_vec(10) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    FFF2_length_ksdensity_grid = 1e3; % now: 1e3, was: 100 choose a number divisible by 4 (200?)
    FFF2_ksdensity_grid        = linspace(0,cl_var.FFF2_stim_end_time,FFF2_length_ksdensity_grid);
    FFF2_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best
    FFF2_spike_density_mat     = NaN(cl_var.True_Num_Cells,FFF2_length_ksdensity_grid);%FFF2_length_ksdensity_grid,cl_var.True_Num_Cells
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.FFF2_spike_times_mat(:,i);cl_var.FFF2_spike_times_mat(:,i);(2*cl_var.FFF2_stim_end_time-cl_var.FFF2_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            [FFF2_spike_density_mat(i,:),~] = ksdensity(vec_loop,FFF2_ksdensity_grid,'Bandwidth',FFF2_ksdensity_bdwth); % [FFF2_spike_density_mat(:,i),~] ksdensity(cl_var.FFF2_spike_times_mat(:,i),FFF2_ksdensity_grid,'Bandwidth',FFF2_ksdensity_bdwth)
        else % If there were no spikes for this cell for this stim.
            FFF2_spike_density_mat(i,:) = zeros(1,FFF2_length_ksdensity_grid);
        end
    end
    % If want to limit to interval:
    % 'Support',[0,cl_var.FFF2_stim_end_time],'BoundaryCorrection','log'/'reflection'(no good as makes go to zero at ends)
    
    % 1D plot of chosen cell FFF2 spike rate
    figure;
    plot(FFF2_ksdensity_grid,FFF2_spike_density_mat(2,:),'LineWidth',1.5); hold on; % FFF2_spike_density_mat(:,2)
    plot(cl_var.FFF2_spike_times_mat(:,2),zeros(size(cl_var.FFF2_spike_times_mat,1),1),'bx');
    FFF2_ks_v1 = vline(cl_var.FFF2_trig_times_vec(2),'r');
    FFF2_ks_v2 = vline(cl_var.FFF2_trig_times_vec(3),'r');
    FFF2_ks_v3 = vline(cl_var.FFF2_trig_times_vec(4),'r');
    FFF2_ks_v4 = vline(cl_var.FFF2_trig_times_vec(5),'r');
    FFF2_ks_v5 = vline(cl_var.FFF2_trig_times_vec(6),'r');
    xlim([0 cl_var.FFF2_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('FFF2 spike rate');
    set(FFF2_ks_v1,'LineWidth',1.5);
    set(FFF2_ks_v2,'LineWidth',1.5);
    set(FFF2_ks_v3,'LineWidth',1.5);
    set(FFF2_ks_v4,'LineWidth',1.5);
    set(FFF2_ks_v5,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell FFF2 spike rates
    FFF2_ks_fig = figure;
    imagesc(FFF2_ksdensity_grid,[],FFF2_spike_density_mat); hold on; % FFF2_spike_density_mat'
    colormap(FFF2_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    FFF2_ks_v1 = vline(cl_var.FFF2_trig_times_vec(2),'r');
    FFF2_ks_v2 = vline(cl_var.FFF2_trig_times_vec(3),'r');
    FFF2_ks_v3 = vline(cl_var.FFF2_trig_times_vec(4),'r');
    FFF2_ks_v4 = vline(cl_var.FFF2_trig_times_vec(5),'r');
    FFF2_ks_v5 = vline(cl_var.FFF2_trig_times_vec(6),'r');
    xlabel('time (sec)');
    ylabel('RF');
    title('FFF2 spike rates');
    set(FFF2_ks_v1,'LineWidth',1.5);
    set(FFF2_ks_v2,'LineWidth',1.5);
    set(FFF2_ks_v3,'LineWidth',1.5);
    set(FFF2_ks_v4,'LineWidth',1.5);
    set(FFF2_ks_v5,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(8) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            %FFF2_segment_time = 0.25*cl_var.FFF2_stim_end_time; % Could have used this...
            
            FFF2_630_spike_density_mat = FFF2_spike_density_mat(:,1:(1/6)*FFF2_length_ksdensity_grid);
            FFF2_560_spike_density_mat = FFF2_spike_density_mat(:,(1/6)*FFF2_length_ksdensity_grid+1:(2/6)*FFF2_length_ksdensity_grid);
            FFF2_505_spike_density_mat = FFF2_spike_density_mat(:,(2/6)*FFF2_length_ksdensity_grid+1:(3/6)*FFF2_length_ksdensity_grid);
            FFF2_480_spike_density_mat = FFF2_spike_density_mat(:,(3/6)*FFF2_length_ksdensity_grid+1:(4/6)*FFF2_length_ksdensity_grid);
            FFF2_430_spike_density_mat = FFF2_spike_density_mat(:,(4/6)*FFF2_length_ksdensity_grid+1:(5/6)*FFF2_length_ksdensity_grid);
            FFF2_360_spike_density_mat = FFF2_spike_density_mat(:,(5/6)*FFF2_length_ksdensity_grid+1:FFF2_length_ksdensity_grid);
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [FFF2_coeff,FFF2_score,FFF2_latent,FFF2_tsquared,FFF2_explained,FFF2_mu] = pca(FFF2_spike_density_mat); %  FFF2_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                FFF2_explained_cumsum  = cumsum(FFF2_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                FFF2_PC_Num = find(FFF2_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                FFF2_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,FFF2_score(:,1:FFF2_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + FFF2_PC_Num;
            
            disp(sprintf('Num FFF2 PCs = %i',FFF2_PC_Num));
            
            FFF2_pca_fig = figure;
            imagesc(FFF2_ksdensity_grid,[],abs(FFF2_coeff(:,1:FFF2_PC_Num)'));%FFF2_PC_Num, 10
            colormap(FFF2_pca_fig,gray(256));
            colorbar;
            FFF2_ks_v1 = vline(cl_var.FFF2_trig_times_vec(2),'r');
            FFF2_ks_v2 = vline(cl_var.FFF2_trig_times_vec(3),'r');
            FFF2_ks_v3 = vline(cl_var.FFF2_trig_times_vec(4),'r');
            FFF2_ks_v4 = vline(cl_var.FFF2_trig_times_vec(5),'r');
            FFF2_ks_v5 = vline(cl_var.FFF2_trig_times_vec(6),'r');
            title('FFF2 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(FFF2_ks_v1,'LineWidth',1.5);
            set(FFF2_ks_v2,'LineWidth',1.5);
            set(FFF2_ks_v3,'LineWidth',1.5);
            set(FFF2_ks_v4,'LineWidth',1.5);
            set(FFF2_ks_v5,'LineWidth',1.5);
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [FFF2_630_coeff,FFF2_630_score,FFF2_630_latent,FFF2_630_tsquared,FFF2_630_explained,FFF2_630_mu] = pca(FFF2_630_spike_density_mat);
            [FFF2_560_coeff,FFF2_560_score,FFF2_560_latent,FFF2_560_tsquared,FFF2_560_explained,FFF2_560_mu] = pca(FFF2_560_spike_density_mat);
            [FFF2_505_coeff,FFF2_505_score,FFF2_505_latent,FFF2_505_tsquared,FFF2_505_explained,FFF2_505_mu] = pca(FFF2_505_spike_density_mat);
            [FFF2_480_coeff,FFF2_480_score,FFF2_480_latent,FFF2_480_tsquared,FFF2_480_explained,FFF2_480_mu] = pca(FFF2_480_spike_density_mat);
            [FFF2_430_coeff,FFF2_430_score,FFF2_430_latent,FFF2_430_tsquared,FFF2_430_explained,FFF2_430_mu] = pca(FFF2_430_spike_density_mat);
            [FFF2_360_coeff,FFF2_360_score,FFF2_360_latent,FFF2_360_tsquared,FFF2_360_explained,FFF2_360_mu] = pca(FFF2_360_spike_density_mat);
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                FFF2_630_explained_cumsum = cumsum(FFF2_630_explained);
                FFF2_560_explained_cumsum = cumsum(FFF2_560_explained);
                FFF2_505_explained_cumsum = cumsum(FFF2_505_explained);
                FFF2_480_explained_cumsum = cumsum(FFF2_480_explained);
                FFF2_430_explained_cumsum = cumsum(FFF2_430_explained);
                FFF2_360_explained_cumsum = cumsum(FFF2_360_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                FFF2_630_PC_Num = find(FFF2_630_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF2_560_PC_Num = find(FFF2_560_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF2_505_PC_Num = find(FFF2_505_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF2_480_PC_Num = find(FFF2_480_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF2_430_PC_Num = find(FFF2_430_explained_cumsum>p.PCA_ExpVar_thresh,1);
                FFF2_360_PC_Num = find(FFF2_360_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                FFF2_630_PC_Num = p.PCA_comp_num;
                FFF2_560_PC_Num = p.PCA_comp_num;
                FFF2_505_PC_Num = p.PCA_comp_num;
                FFF2_480_PC_Num = p.PCA_comp_num;
                FFF2_430_PC_Num = p.PCA_comp_num;
                FFF2_360_PC_Num = p.PCA_comp_num;
                
            end
            
            FFF2_PC_Num = FFF2_630_PC_Num + FFF2_560_PC_Num + FFF2_505_PC_Num + FFF2_480_PC_Num + FFF2_430_PC_Num + FFF2_360_PC_Num;
            
            All_Scores     = [All_Scores,FFF2_630_score(:,1:FFF2_630_PC_Num),FFF2_560_score(:,1:FFF2_560_PC_Num),...
                FFF2_505_score(:,1:FFF2_505_PC_Num),FFF2_480_score(:,1:FFF2_480_PC_Num),FFF2_430_score(:,1:FFF2_430_PC_Num),FFF2_360_score(:,1:FFF2_360_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + FFF2_PC_Num;
            
            disp(sprintf('Num FFF2 R PCs = %i',FFF2_630_PC_Num));
            disp(sprintf('Num FFF2 G PCs = %i',FFF2_560_PC_Num));
            disp(sprintf('Num FFF2 B PCs = %i',FFF2_505_PC_Num));
            disp(sprintf('Num FFF2 UV PCs = %i',FFF2_480_PC_Num));
            disp(sprintf('Num FFF2 B PCs = %i',FFF2_430_PC_Num));
            disp(sprintf('Num FFF2 UV PCs = %i',FFF2_360_PC_Num));
            disp(sprintf('Total Num FFF2 PCs = %i',FFF2_PC_Num));
            
            FFF2_pca_fig = figure;
            subplot(2,3,1);
            imagesc(FFF2_ksdensity_grid(1:(1/6)*FFF2_length_ksdensity_grid),[],abs(FFF2_630_coeff(:,1:FFF2_630_PC_Num)'));
            colormap(FFF2_pca_fig,gray(256));
            colorbar;
            title('FFF2 630 - PCA coeff');ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,3,2);
            imagesc(FFF2_ksdensity_grid((1/6)*FFF2_length_ksdensity_grid+1:(2/6)*FFF2_length_ksdensity_grid),[],abs(FFF2_560_coeff(:,1:FFF2_560_PC_Num)'));
            colormap(FFF2_pca_fig,gray(256));
            colorbar;
            title('FFF2 560 - PCA coeff');
            set(gca,'FontSize',12);
            subplot(2,3,3);
            imagesc(FFF2_ksdensity_grid((2/6)*FFF2_length_ksdensity_grid+1:(3/6)*FFF2_length_ksdensity_grid),[],abs(FFF2_505_coeff(:,1:FFF2_505_PC_Num)'));
            colormap(FFF2_pca_fig,gray(256));
            colorbar;
            title('FFF2 505 - PCA coeff');
            set(gca,'FontSize',12);
            subplot(2,3,4);
            imagesc(FFF2_ksdensity_grid((3/6)*FFF2_length_ksdensity_grid+1:(4/6)*FFF2_length_ksdensity_grid),[],abs(FFF2_480_coeff(:,1:FFF2_480_PC_Num)'));
            colormap(FFF2_pca_fig,gray(256));
            colorbar;
            title('FFF2 480 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(2,3,5);
            imagesc(FFF2_ksdensity_grid((4/6)*FFF2_length_ksdensity_grid+1:(5/6)*FFF2_length_ksdensity_grid),[],abs(FFF2_430_coeff(:,1:FFF2_430_PC_Num)'));
            colormap(FFF2_pca_fig,gray(256));
            colorbar;
            title('FFF2 430 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            subplot(2,3,6);
            imagesc(FFF2_ksdensity_grid((5/6)*FFF2_length_ksdensity_grid+1:FFF2_length_ksdensity_grid),[],abs(FFF2_360_coeff(:,1:FFF2_360_PC_Num)'));
            colormap(FFF2_pca_fig,gray(256));
            colorbar;
            title('FFF2 360 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% Chirp2
if p.Obj_clust_vec(9) == 1 || p.Obj_plot_vec(11) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    Chirp2_length_ksdensity_grid = 1e3; % now: 1e3, was: 100
    Chirp2_ksdensity_grid        = linspace(0,cl_var.Chirp2_stim_end_time,Chirp2_length_ksdensity_grid);
    Chirp2_ksdensity_bdwth       = 5*1e-2; % Tom said between 1e-2 and 1e-1 is best for original Chirp
    Chirp2_spike_density_mat     = NaN(cl_var.True_Num_Cells,Chirp2_length_ksdensity_grid);
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.Chirp2_spike_times_mat(:,i);cl_var.Chirp2_spike_times_mat(:,i);(2*cl_var.Chirp2_stim_end_time-cl_var.Chirp2_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            %[Chirp2_spike_density_mat(i,:),~] = ksdensity(vec_loop,Chirp2_ksdensity_grid,'Bandwidth',Chirp2_ksdensity_bdwth); % ksdensity(cl_var.Chirp2_spike_times_mat(:,i),Chirp2_ksdensity_grid,'Bandwidth',Chirp2_ksdensity_bdwth);
            [Chirp2_spike_density_temp,~] = ksdensity(vec_loop,Chirp2_ksdensity_grid,'Bandwidth',Chirp2_ksdensity_bdwth); % Mod 14,06,2021
            %Chirp2_spike_density_mat(i,:) = Chirp2_spike_density_temp/max(Chirp2_spike_density_temp);                     % Mod 14,06,2021
            Chirp2_spike_density_mat(i,:) = Chirp2_spike_density_temp*sum(~isnan(cl_var.Chirp2_spike_times_mat(:,i)));                     % Mod 14,06,2021
        else % If there were no spikes for this cell for this stim.
            Chirp2_spike_density_mat(i,:) = zeros(1,Chirp2_length_ksdensity_grid);
        end
    end
    
    % Not necc as have 0 cl_var.Chirp2_stim_end_time
    % % Spike density time start and end
    % Chirp2_t_start = Chirp2_ksdensity_grid(1);
    % Chirp2_t_end   = Chirp2_ksdensity_grid(end);
    
    % 1D plot of chosen cell Chirp2 spike rate
    figure;
    plot(Chirp2_ksdensity_grid,Chirp2_spike_density_mat(2,:),'LineWidth',1.5); hold on;
    plot(cl_var.Chirp2_spike_times_mat(:,2),zeros(size(cl_var.Chirp2_spike_times_mat,1),1),'bx');
    Chirp2_ks_v1 = vline(cl_var.Chirp2_trig_times_vec(2),'r');
    xlim([0 cl_var.Chirp2_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('Chirp2 spike rate');
    set(Chirp2_ks_v1,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell Chirp2 spike rates
    Chirp2_ks_fig = figure;
    imagesc(Chirp2_ksdensity_grid,[],Chirp2_spike_density_mat); hold on;
    colormap(Chirp2_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    Chirp2_ks_v1 = vline(cl_var.Chirp2_trig_times_vec(2),'r');
    xlabel('time (sec)');
    ylabel('RF');
    title('Chirp2 spike rates');
    set(Chirp2_ks_v1,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(9) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            Num_Chirp2_1_bins = sum(Chirp2_ksdensity_grid<cl_var.Chirp2_trig_times_vec(2));
            
            Chirp2_1_spike_density_mat  = Chirp2_spike_density_mat(:,1:Num_Chirp2_1_bins);
            Chirp2_2_spike_density_mat  = Chirp2_spike_density_mat(:,Num_Chirp2_1_bins+1:end);
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [Chirp2_coeff,Chirp2_score,Chirp2_latent,Chirp2_tsquared,Chirp2_explained,Chirp2_mu] = pca(Chirp2_spike_density_mat); %  Chirp2_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Chirp2_explained_cumsum  = cumsum(Chirp2_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Chirp2_PC_Num = find(Chirp2_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Chirp2_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,Chirp2_score(:,1:Chirp2_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Chirp2_PC_Num;
            
            disp(sprintf('Num Chirp2 PCs = %i',Chirp2_PC_Num));
            
            Chirp2_pca_fig = figure;
            imagesc(Chirp2_ksdensity_grid,[],abs(Chirp2_coeff(:,1:Chirp2_PC_Num)'));%Chirp2_PC_Num, 10
            colormap(Chirp2_pca_fig,gray(256));
            colorbar;
            Chirp2_ks_v1 = vline(cl_var.Chirp2_trig_times_vec(2),'r');
            title('Chirp2 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(Chirp2_ks_v1,'LineWidth',1.5);
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [Chirp2_1_coeff,Chirp2_1_score,Chirp2_1_latent,Chirp2_1_tsquared,Chirp2_1_explained,Chirp2_1_mu] = pca(Chirp2_1_spike_density_mat);
            [Chirp2_2_coeff,Chirp2_2_score,Chirp2_2_latent,Chirp2_2_tsquared,Chirp2_2_explained,Chirp2_2_mu] = pca(Chirp2_2_spike_density_mat);
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Chirp2_1_explained_cumsum  = cumsum(Chirp2_1_explained);
                Chirp2_2_explained_cumsum  = cumsum(Chirp2_2_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Chirp2_1_PC_Num  = find(Chirp2_1_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Chirp2_2_PC_Num  = find(Chirp2_2_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Chirp2_1_PC_Num  = p.PCA_comp_num;
                Chirp2_2_PC_Num  = p.PCA_comp_num;
                
            end
            
            Chirp2_PC_Num = Chirp2_1_PC_Num + Chirp2_2_PC_Num;
            
            All_Scores     = [All_Scores,Chirp2_1_score(:,1:Chirp2_1_PC_Num),Chirp2_2_score(:,1:Chirp2_2_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Chirp2_PC_Num;
            
            disp(sprintf('Num Chirp2 1 PCs = %i',Chirp2_1_PC_Num));
            disp(sprintf('Num Chirp2 2 PCs = %i',Chirp2_2_PC_Num));
            disp(sprintf('Total Num Chirp2 PCs = %i',Chirp2_PC_Num));
            
            Chirp2_pca_fig = figure;
            subplot(1,2,1);
            imagesc(Chirp2_ksdensity_grid(1:Num_Chirp2_1_bins),[],abs(Chirp2_1_coeff(:,1:Chirp2_1_PC_Num)'));
            colormap(Chirp2_pca_fig,gray(256));
            colorbar;
            title('Chirp2 1 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(1,2,2);
            imagesc(Chirp2_ksdensity_grid(Num_Chirp2_1_bins+1:end),[],abs(Chirp2_2_coeff(:,1:Chirp2_2_PC_Num)'));
            colormap(Chirp2_pca_fig,gray(256));
            title('Chirp2 2 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% SSub
if p.Obj_clust_vec(10) == 1 || p.Obj_plot_vec(12) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    SSub_length_ksdensity_grid = 1.8*1e3; % now: 1e3, was: 100 choose a number divisible by 18 (200?)
    SSub_ksdensity_grid        = linspace(0,cl_var.SSub_stim_end_time,SSub_length_ksdensity_grid);
    SSub_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best for FFF
    SSub_spike_density_mat     = NaN(cl_var.True_Num_Cells,SSub_length_ksdensity_grid);%SSub_length_ksdensity_grid,cl_var.True_Num_Cells
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.SSub_spike_times_mat(:,i);cl_var.SSub_spike_times_mat(:,i);(2*cl_var.SSub_stim_end_time-cl_var.SSub_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            [SSub_spike_density_mat(i,:),~] = ksdensity(vec_loop,SSub_ksdensity_grid,'Bandwidth',SSub_ksdensity_bdwth); % [SSub_spike_density_mat(:,i),~] ksdensity(cl_var.SSub_spike_times_mat(:,i),SSub_ksdensity_grid,'Bandwidth',SSub_ksdensity_bdwth)
        else % If there were no spikes for this cell for this stim.
            SSub_spike_density_mat(i,:) = zeros(1,SSub_length_ksdensity_grid);
        end
    end
    % If want to limit to interval:
    % 'Support',[0,cl_var.SSub_stim_end_time],'BoundaryCorrection','log'/'reflection'(no good as makes go to zero at ends)
    
    % 1D plot of chosen cell SSub spike rate
    figure;
    plot(SSub_ksdensity_grid,SSub_spike_density_mat(2,:),'LineWidth',1.5); hold on; % SSub_spike_density_mat(:,2)
    plot(cl_var.SSub_spike_times_mat(:,2),zeros(size(cl_var.SSub_spike_times_mat,1),1),'bx');
    for i = 1:p.SSub_NumTrigPerRep-1
        SSub_ks_loop = vline(cl_var.SSub_trig_times_vec(i+1),'r');
        set(SSub_ks_loop,'LineWidth',1.5);
    end
    xlim([0 cl_var.SSub_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('SSub spike rate');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell SSub spike rates
    SSub_ks_fig = figure;
    imagesc(SSub_ksdensity_grid,[],SSub_spike_density_mat); hold on; % SSub_spike_density_mat'
    colormap(SSub_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    for i = 1:p.SSub_NumTrigPerRep-1
        SSub_ks_loop = vline(cl_var.SSub_trig_times_vec(i+1),'r');
        set(SSub_ks_loop,'LineWidth',1.5);
    end
    xlabel('time (sec)');
    ylabel('RF');
    title('SSub spike rates');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(10) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            %SSub_segment_time = 0.25*cl_var.SSub_stim_end_time; % Could have used this...
            
            SSub_spike_density_arr = NaN(p.SSub_NumTrigPerRep,cl_var.True_Num_Cells,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep);
            
            for i = 1:p.SSub_NumTrigPerRep
                SSub_spike_density_arr(i,:,:) = SSub_spike_density_mat(:,((i-1)/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid+1:(i/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid);
            end
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [SSub_coeff,SSub_score,SSub_latent,SSub_tsquared,SSub_explained,SSub_mu] = pca(SSub_spike_density_mat); %  SSub_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                SSub_explained_cumsum  = cumsum(SSub_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                SSub_PC_Num = find(SSub_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                SSub_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,SSub_score(:,1:SSub_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + SSub_PC_Num;
            
            disp(sprintf('Num SSub PCs = %i',SSub_PC_Num));
            
            SSub_pca_fig = figure;
            imagesc(SSub_ksdensity_grid,[],abs(SSub_coeff(:,1:SSub_PC_Num)'));%SSub_PC_Num, 10
            colormap(SSub_pca_fig,gray(256));
            colorbar;
            for i = 1:p.SSub_NumTrigPerRep-1
                SSub_ks_loop = vline(cl_var.SSub_trig_times_vec(i+1),'r');
                set(SSub_ks_loop,'LineWidth',1.5);
            end
            title('SSub - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [SSub_coeff_temp,SSub_score_temp,SSub_latent_temp,SSub_tsquared_temp,SSub_explained_temp,SSub_mu_temp] = pca(squeeze(SSub_spike_density_arr(1,:,:)));
            [SSub_coeff_rows,SSub_coeff_cols] = size(SSub_coeff_temp);
            [SSub_score_rows,SSub_score_cols] = size(SSub_score_temp);
            SSub_latent_length                = length(SSub_latent_temp);
            SSub_tsquared_length              = length(SSub_tsquared_temp);
            SSub_explained_length             = length(SSub_explained_temp);
            SSub_mu_length                    = length(SSub_mu_temp);
            
%             SSub_coeff_arr     = NaN(p.SSub_NumTrigPerRep,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep); % pxp = varxvar, p = num time pts (was cl_var.True_Num_Cells-1 at end)
%             SSub_score_arr     = NaN(p.SSub_NumTrigPerRep,cl_var.True_Num_Cells,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep);   % num cells x num var (comp) (was cl_var.True_Num_Cells-1 at end)
%             SSub_latent_arr    = NaN(p.SSub_NumTrigPerRep,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep);                         % num var (was cl_var.True_Num_Cells-1 at end)
%             SSub_tsquared_arr  = NaN(p.SSub_NumTrigPerRep,cl_var.True_Num_Cells);                                                   % num cells
%             SSub_explained_arr = NaN(p.SSub_NumTrigPerRep,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep);                         % num comp (var) (was cl_var.True_Num_Cells-1 at end)
%             SSub_mu_arr        = NaN(p.SSub_NumTrigPerRep,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep);                         % num comp (var)

            SSub_coeff_arr     = NaN(p.SSub_NumTrigPerRep,SSub_coeff_rows,SSub_coeff_cols); % pxp = varxvar, p = num time pts
            SSub_score_arr     = NaN(p.SSub_NumTrigPerRep,SSub_score_rows,SSub_score_cols); % num cells x num var (comp)
            SSub_latent_arr    = NaN(p.SSub_NumTrigPerRep,SSub_latent_length);              % num var
            SSub_tsquared_arr  = NaN(p.SSub_NumTrigPerRep,SSub_tsquared_length);            % num cells
            SSub_explained_arr = NaN(p.SSub_NumTrigPerRep,SSub_explained_length);           % num comp (var) (was cl_var.True_Num_Cells-1 at end)
            SSub_mu_arr        = NaN(p.SSub_NumTrigPerRep,SSub_mu_length);                  % num comp (var)
            
            
            for i = 1:p.SSub_NumTrigPerRep
                [SSub_coeff_arr(i,:,:),SSub_score_arr(i,:,:),SSub_latent_arr(i,:),SSub_tsquared_arr(i,:),SSub_explained_arr(i,:),SSub_mu_arr(i,:)] = pca(squeeze(SSub_spike_density_arr(i,:,:)));
%                 [SSub_coeff_loop,SSub_score_arr_loop,SSub_latent_arr_loop,SSub_tsquared_arr_loop,SSub_explained_arr_loop,SSub_mu_arr_loop] = pca(squeeze(SSub_spike_density_arr(i,:,:)));
%                 SSub_coeff_arr(i,:,:)   = SSub_coeff_loop;
%                 SSub_score_arr(i,:,:)   = SSub_score_arr_loop;
%                 SSub_latent_arr(i,:)    = SSub_latent_arr_loop;
%                 SSub_tsquared_arr(i,:)  = SSub_tsquared_arr_loop;
%                 SSub_explained_arr(i,:) = SSub_explained_arr_loop;
%                 SSub_mu_arr(i,:)        = SSub_mu_arr_loop;
            end
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                %SSub_explained_cumsum_arr = NaN(p.SSub_NumTrigPerRep,SSub_length_ksdensity_grid/p.SSub_NumTrigPerRep); % (was cl_var.True_Num_Cells-1 at end)
                SSub_explained_cumsum_arr = NaN(p.SSub_NumTrigPerRep,SSub_explained_length);
                for i = 1:p.SSub_NumTrigPerRep
                    SSub_explained_cumsum_arr(i,:) = cumsum(SSub_explained_arr(i,:));
                end
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                
                SSub_PC_Num_vec = NaN(p.SSub_NumTrigPerRep,1);
                for i = 1:p.SSub_NumTrigPerRep
                    SSub_PC_Num_vec(i) = find(SSub_explained_cumsum_arr(i,:)>p.PCA_ExpVar_thresh,1);
                end
                
            else % p.PCA_thresh_type == 2 % number of components
                
                SSub_PC_Num_vec = p.PCA_comp_num*ones(p.SSub_NumTrigPerRep,1);
                
            end
            
            SSub_PC_Num = sum(SSub_PC_Num_vec);
            
            for i = 1:p.SSub_NumTrigPerRep
                if SSub_PC_Num_vec(i) > 1
                    All_Scores = [All_Scores,squeeze(SSub_score_arr(i,:,1:SSub_PC_Num_vec(i)))];
                else % SSub_PC_Num_vec(i) == 1 Need to take transpose in thise case
                    All_Scores = [All_Scores,squeeze(SSub_score_arr(i,:,1:SSub_PC_Num_vec(i)))'];
                end
            end
            
            Total_Num_Comp = Total_Num_Comp + SSub_PC_Num;
            
            for i = 1:p.SSub_NumTrigPerRep
                disp(sprintf('Num SSub PCs %i = %i',i,SSub_PC_Num_vec(i)));
            end
            
            
            SSub_PC_Plot_Dim = ceil(sqrt(p.SSub_NumTrigPerRep));
            SSub_PC_Plot_Row = ceil(p.SSub_NumTrigPerRep/SSub_PC_Plot_Dim);
            SSub_PC_Plot_Col = SSub_PC_Plot_Dim;
            
            SSub_pca_fig = figure;
            for i = 1:p.SSub_NumTrigPerRep
                subplot(SSub_PC_Plot_Row,SSub_PC_Plot_Col,i);
                if SSub_PC_Num_vec(i)>1
                            %SSub_spike_density_arr(i,:,:) = SSub_spike_density_mat(:,((i-1)/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid+1:(i/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid); SSub_coeff_arr(i,:,:)
                    imagesc(SSub_ksdensity_grid(((i-1)/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid+1:(i/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid),[],abs(squeeze(SSub_coeff_arr(i,:,1:SSub_PC_Num_vec(i)))'));
                else
                    imagesc(SSub_ksdensity_grid(((i-1)/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid+1:(i/p.SSub_NumTrigPerRep)*SSub_length_ksdensity_grid),[],abs(squeeze(SSub_coeff_arr(i,:,1:SSub_PC_Num_vec(i)))));
                end
                colormap(SSub_pca_fig,gray(256));
                colorbar;
                %title('SSub 630 - PCA coeff');
                if (i/SSub_PC_Plot_Col > SSub_PC_Plot_Row-1) || ((i/SSub_PC_Plot_Col > SSub_PC_Plot_Row-2)&&(i/SSub_PC_Plot_Col < SSub_PC_Plot_Row-1)&&(mod(p.SSub_NumTrigPerRep,SSub_PC_Plot_Col)~=0)&&(mod(i,SSub_PC_Plot_Col)>mod(p.SSub_NumTrigPerRep,SSub_PC_Plot_Col)))
                    xlabel('time (sec)');
                end
                if mod(i,SSub_PC_Plot_Col) == 1
                    ylabel('component');
                end
                set(gca,'FontSize',12);
            end
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% CSteps
if p.Obj_clust_vec(11) == 1 || p.Obj_plot_vec(13) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    CSteps_length_ksdensity_grid = 1e3; % now: 1e3, was: 100 choose a number divisible by 10 (200?)
    CSteps_ksdensity_grid        = linspace(0,cl_var.CSteps_stim_end_time,CSteps_length_ksdensity_grid);
    CSteps_ksdensity_bdwth       = 5*1e-2; % Tom said 5*1e-2 is best for FFF
    CSteps_spike_density_mat     = NaN(cl_var.True_Num_Cells,CSteps_length_ksdensity_grid);%CSteps_length_ksdensity_grid,cl_var.True_Num_Cells
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.CSteps_spike_times_mat(:,i);cl_var.CSteps_spike_times_mat(:,i);(2*cl_var.CSteps_stim_end_time-cl_var.CSteps_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            [CSteps_spike_density_mat(i,:),~] = ksdensity(vec_loop,CSteps_ksdensity_grid,'Bandwidth',CSteps_ksdensity_bdwth); % [CSteps_spike_density_mat(:,i),~] ksdensity(cl_var.CSteps_spike_times_mat(:,i),CSteps_ksdensity_grid,'Bandwidth',CSteps_ksdensity_bdwth)
        else % If there were no spikes for this cell for this stim.
            CSteps_spike_density_mat(i,:) = zeros(1,CSteps_length_ksdensity_grid);
        end
    end
    % If want to limit to interval:
    % 'Support',[0,cl_var.CSteps_stim_end_time],'BoundaryCorrection','log'/'reflection'(no good as makes go to zero at ends)
    
    % 1D plot of chosen cell CSteps spike rate
    figure;
    plot(CSteps_ksdensity_grid,CSteps_spike_density_mat(2,:),'LineWidth',1.5); hold on; % CSteps_spike_density_mat(:,2)
    plot(cl_var.CSteps_spike_times_mat(:,2),zeros(size(cl_var.CSteps_spike_times_mat,1),1),'bx');
    for i = 1:p.CSteps_NumTrigPerRep-1
        CSteps_ks_loop = vline(cl_var.CSteps_trig_times_vec(i+1),'r');
        set(CSteps_ks_loop,'LineWidth',1.5);
    end
    xlim([0 cl_var.CSteps_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('CSteps spike rate');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell CSteps spike rates
    CSteps_ks_fig = figure;
    imagesc(CSteps_ksdensity_grid,[],CSteps_spike_density_mat); hold on; % CSteps_spike_density_mat'
    colormap(CSteps_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    for i = 1:p.CSteps_NumTrigPerRep-1
        CSteps_ks_loop = vline(cl_var.CSteps_trig_times_vec(i+1),'r');
        set(CSteps_ks_loop,'LineWidth',1.5);
    end
    xlabel('time (sec)');
    ylabel('RF');
    title('CSteps spike rates');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(11) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            %CSteps_segment_time = 0.25*cl_var.CSteps_stim_end_time; % Could have used this...
            
            CSteps_spike_density_arr = NaN(p.CSteps_NumTrigPerRep,cl_var.True_Num_Cells,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep);
            
            for i = 1:p.CSteps_NumTrigPerRep
                CSteps_spike_density_arr(i,:,:) = CSteps_spike_density_mat(:,((i-1)/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid+1:(i/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid);
            end
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [CSteps_coeff,CSteps_score,CSteps_latent,CSteps_tsquared,CSteps_explained,CSteps_mu] = pca(CSteps_spike_density_mat); %  CSteps_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                CSteps_explained_cumsum  = cumsum(CSteps_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                CSteps_PC_Num = find(CSteps_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                CSteps_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,CSteps_score(:,1:CSteps_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + CSteps_PC_Num;
            
            disp(sprintf('Num CSteps PCs = %i',CSteps_PC_Num));
            
            CSteps_pca_fig = figure;
            imagesc(CSteps_ksdensity_grid,[],abs(CSteps_coeff(:,1:CSteps_PC_Num)'));%CSteps_PC_Num, 10
            colormap(CSteps_pca_fig,gray(256));
            colorbar;
            for i = 1:p.CSteps_NumTrigPerRep-1
                CSteps_ks_loop = vline(cl_var.CSteps_trig_times_vec(i+1),'r');
                set(CSteps_ks_loop,'LineWidth',1.5);
            end
            title('CSteps - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
                
            [CSteps_coeff_temp,CSteps_score_temp,CSteps_latent_temp,CSteps_tsquared_temp,CSteps_explained_temp,CSteps_mu_temp] = pca(squeeze(CSteps_spike_density_arr(1,:,:)));
            [CSteps_coeff_rows,CSteps_coeff_cols] = size(CSteps_coeff_temp);
            [CSteps_score_rows,CSteps_score_cols] = size(CSteps_score_temp);
            CSteps_latent_length                  = length(CSteps_latent_temp);
            CSteps_tsquared_length                = length(CSteps_tsquared_temp);
            CSteps_explained_length               = length(CSteps_explained_temp);
            CSteps_mu_length                      = length(CSteps_mu_temp);
            
%             CSteps_coeff_arr     = NaN(p.CSteps_NumTrigPerRep,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep); % pxp = varxvar, p = num time pts
%             CSteps_score_arr     = NaN(p.CSteps_NumTrigPerRep,cl_var.True_Num_Cells,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep); % num cells x num var (comp)
%             CSteps_latent_arr    = NaN(p.CSteps_NumTrigPerRep,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep); % num var
%             CSteps_tsquared_arr  = NaN(p.CSteps_NumTrigPerRep,cl_var.True_Num_Cells);                               % num cells
%             CSteps_explained_arr = NaN(p.CSteps_NumTrigPerRep,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep); % num comp (var)
%             CSteps_mu_arr        = NaN(p.CSteps_NumTrigPerRep,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep); % num comp (var)
            
              CSteps_coeff_arr     = NaN(p.CSteps_NumTrigPerRep,CSteps_coeff_rows,CSteps_coeff_cols); % pxp = varxvar, p = num time pts
              CSteps_score_arr     = NaN(p.CSteps_NumTrigPerRep,CSteps_score_rows,CSteps_score_cols); % num cells x num var (comp)
              CSteps_latent_arr    = NaN(p.CSteps_NumTrigPerRep,CSteps_latent_length);                % num var
              CSteps_tsquared_arr  = NaN(p.CSteps_NumTrigPerRep,CSteps_tsquared_length);              % num cells
              CSteps_explained_arr = NaN(p.CSteps_NumTrigPerRep,CSteps_explained_length);             % num comp (var) (was cl_var.True_Num_Cells-1 at end)
              CSteps_mu_arr        = NaN(p.CSteps_NumTrigPerRep,CSteps_mu_length);                    % num comp (var)

            for i = 1:p.CSteps_NumTrigPerRep
                [CSteps_coeff_arr(i,:,:),CSteps_score_arr(i,:,:),CSteps_latent_arr(i,:),CSteps_tsquared_arr(i,:),CSteps_explained_arr(i,:),CSteps_mu_arr(i,:)] = pca(squeeze(CSteps_spike_density_arr(i,:,:)));
%                 [CSteps_coeff_loop,CSteps_score_arr_loop,CSteps_latent_arr_loop,CSteps_tsquared_arr_loop,CSteps_explained_arr_loop,CSteps_mu_arr_loop] = pca(squeeze(CSteps_spike_density_arr(i,:,:)));
%                 CSteps_coeff_arr(i,:,:)   = CSteps_coeff_loop;
%                 CSteps_score_arr(i,:,:)   = CSteps_score_arr_loop;
%                 CSteps_latent_arr(i,:)    = CSteps_latent_arr_loop;
%                 CSteps_tsquared_arr(i,:)  = CSteps_tsquared_arr_loop;
%                 CSteps_explained_arr(i,:) = CSteps_explained_arr_loop;
%                 CSteps_mu_arr(i,:)        = CSteps_mu_arr_loop;
            end
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                %CSteps_explained_cumsum_arr = NaN(p.CSteps_NumTrigPerRep,CSteps_length_ksdensity_grid/p.CSteps_NumTrigPerRep);%cl_var.True_Num_Cells-1
                CSteps_explained_cumsum_arr = NaN(p.CSteps_NumTrigPerRep,CSteps_explained_length);
                for i = 1:p.CSteps_NumTrigPerRep
                    CSteps_explained_cumsum_arr(i,:) = cumsum(CSteps_explained_arr(i,:));
                end
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                
                CSteps_PC_Num_vec = NaN(p.CSteps_NumTrigPerRep,1);
                for i = 1:p.CSteps_NumTrigPerRep
                    CSteps_PC_Num_vec(i) = find(CSteps_explained_cumsum_arr(i,:)>p.PCA_ExpVar_thresh,1);
                end
                
            else % p.PCA_thresh_type == 2 % number of components
                
                CSteps_PC_Num_vec = p.PCA_comp_num*ones(p.CSteps_NumTrigPerRep,1);
                
            end
            
            CSteps_PC_Num = sum(CSteps_PC_Num_vec);
            
            for i = 1:p.CSteps_NumTrigPerRep
                if CSteps_PC_Num_vec(i) > 1
                    All_Scores = [All_Scores,squeeze(CSteps_score_arr(i,:,1:CSteps_PC_Num_vec(i)))];
                else % CSteps_PC_Num_vec(i) == 1 Need to take transpose in thise case
                    All_Scores = [All_Scores,squeeze(CSteps_score_arr(i,:,1:CSteps_PC_Num_vec(i)))'];
                end
            end
            
            Total_Num_Comp = Total_Num_Comp + CSteps_PC_Num;
            
            for i = 1:p.CSteps_NumTrigPerRep
                disp(sprintf('Num CSteps PCs %i = %i',i,CSteps_PC_Num_vec(i)));
            end
            
            
            CSteps_PC_Plot_Dim = ceil(sqrt(p.CSteps_NumTrigPerRep));
            CSteps_PC_Plot_Row = ceil(p.CSteps_NumTrigPerRep/CSteps_PC_Plot_Dim);
            CSteps_PC_Plot_Col = CSteps_PC_Plot_Dim;
            
            CSteps_pca_fig = figure;
            for i = 1:p.CSteps_NumTrigPerRep
                subplot(CSteps_PC_Plot_Row,CSteps_PC_Plot_Col,i);
                if CSteps_PC_Num_vec(i)>1
                       %CSteps_spike_density_arr(i,:,:) = CSteps_spike_density_mat(:,((i-1)/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid+1:(i/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid); CSteps_coeff_arr(i,:,:)
                       imagesc(CSteps_ksdensity_grid(((i-1)/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid+1:(i/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid),[],abs(squeeze(CSteps_coeff_arr(i,:,1:CSteps_PC_Num_vec(i)))'));
                else
                    imagesc(CSteps_ksdensity_grid(((i-1)/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid+1:(i/p.CSteps_NumTrigPerRep)*CSteps_length_ksdensity_grid),[],abs(squeeze(CSteps_coeff_arr(i,:,1:CSteps_PC_Num_vec(i)))));
                end
                colormap(CSteps_pca_fig,gray(256));
                colorbar;
                %title('CSteps 630 - PCA coeff');
                if (i/CSteps_PC_Plot_Col > CSteps_PC_Plot_Row-1) || ((i/CSteps_PC_Plot_Col > CSteps_PC_Plot_Row-2)&&(i/CSteps_PC_Plot_Col < CSteps_PC_Plot_Row-1)&&(mod(p.CSteps_NumTrigPerRep,CSteps_PC_Plot_Col)~=0)&&(mod(i,CSteps_PC_Plot_Col)>mod(p.CSteps_NumTrigPerRep,CSteps_PC_Plot_Col)))
                    xlabel('time (sec)');
                end
                if mod(i,CSteps_PC_Plot_Col) == 1
                    ylabel('component');
                end
                set(gca,'FontSize',12);
            end
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end


%%% Chirp3
if p.Obj_clust_vec(12) == 1 || p.Obj_plot_vec(14) == 1
    
    %%% 4. Apply Kernel Density Smoothing to Each Cell
    Chirp3_length_ksdensity_grid = 1e3; % now: 1e3, was: 100
    Chirp3_ksdensity_grid        = linspace(0,cl_var.Chirp3_stim_end_time,Chirp3_length_ksdensity_grid);
    Chirp3_ksdensity_bdwth       = 5*1e-2; % Tom said between 1e-2 and 1e-1 is best for original Chirp
    Chirp3_spike_density_mat     = NaN(cl_var.True_Num_Cells,Chirp3_length_ksdensity_grid);
    
    for i = 1:cl_var.True_Num_Cells
        vec_loop = [-cl_var.Chirp3_spike_times_mat(:,i);cl_var.Chirp3_spike_times_mat(:,i);(2*cl_var.Chirp3_stim_end_time-cl_var.Chirp3_spike_times_mat(:,i))];
        if any(~isnan(vec_loop)) % If there were spikes for this cell for this stim.
            %[Chirp3_spike_density_mat(i,:),~] = ksdensity(vec_loop,Chirp3_ksdensity_grid,'Bandwidth',Chirp3_ksdensity_bdwth); % ksdensity(cl_var.Chirp3_spike_times_mat(:,i),Chirp3_ksdensity_grid,'Bandwidth',Chirp3_ksdensity_bdwth);
            [Chirp3_spike_density_temp,~] = ksdensity(vec_loop,Chirp3_ksdensity_grid,'Bandwidth',Chirp3_ksdensity_bdwth); % Mod 14,06,2021
            %Chirp3_spike_density_mat(i,:) = Chirp3_spike_density_temp/max(Chirp3_spike_density_temp);                     % Mod 14,06,2021
            Chirp3_spike_density_mat(i,:) = Chirp3_spike_density_temp*sum(~isnan(cl_var.Chirp3_spike_times_mat(:,i)));                     % Mod 14,06,2021
        else % If there were no spikes for this cell for this stim.
            Chirp3_spike_density_mat(i,:) = zeros(1,Chirp3_length_ksdensity_grid);
        end
    end
    
    % Not necc as have 0 cl_var.Chirp3_stim_end_time
    % % Spike density time start and end
    % Chirp3_t_start = Chirp3_ksdensity_grid(1);
    % Chirp3_t_end   = Chirp3_ksdensity_grid(end);
    
    % 1D plot of chosen cell Chirp3 spike rate
    figure;
    plot(Chirp3_ksdensity_grid,Chirp3_spike_density_mat(2,:),'LineWidth',1.5); hold on;
    plot(cl_var.Chirp3_spike_times_mat(:,2),zeros(size(cl_var.Chirp3_spike_times_mat,1),1),'bx');
    Chirp3_ks_v1 = vline(cl_var.Chirp3_trig_times_vec(2),'r');
    Chirp3_ks_v2 = vline(cl_var.Chirp3_trig_times_vec(3),'r');
    xlim([0 cl_var.Chirp3_stim_end_time]);
    xlabel('time (sec)');
    ylabel('spike rate');
    title('Chirp3 spike rate');
    set(Chirp3_ks_v1,'LineWidth',1.5);
    set(Chirp3_ks_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    % 2D plot of all cell Chirp3 spike rates
    Chirp3_ks_fig = figure;
    imagesc(Chirp3_ksdensity_grid,[],Chirp3_spike_density_mat); hold on;
    colormap(Chirp3_ks_fig,gray(256)); % gray(256), parula(256)
    colorbar;
    Chirp3_ks_v1 = vline(cl_var.Chirp3_trig_times_vec(2),'r');
    Chirp3_ks_v2 = vline(cl_var.Chirp3_trig_times_vec(3),'r');
    xlabel('time (sec)');
    ylabel('RF');
    title('Chirp3 spike rates');
    set(Chirp3_ks_v1,'LineWidth',1.5);
    set(Chirp3_ks_v2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
    
    
    %%% 5. Perform PCA
    
    if p.Obj_clust_vec(12) == 1 % Only do this if clustering on it, not if just plotting (11,05,2021)
        
        if p.Dim_red_meth == 2 || p.Dim_red_meth == 4 % PCA segmented or sparse PCA segmented
            
            Num_Chirp3_1_bins = sum(Chirp3_ksdensity_grid<cl_var.Chirp3_trig_times_vec(2));
            Num_Chirp3_2_bins = sum(Chirp3_ksdensity_grid<cl_var.Chirp3_trig_times_vec(3))-Num_Chirp3_1_bins;
            
            Chirp3_1_spike_density_mat  = Chirp3_spike_density_mat(:,1:Num_Chirp3_1_bins);
            Chirp3_2_spike_density_mat  = Chirp3_spike_density_mat(:,Num_Chirp3_1_bins+1:Num_Chirp3_1_bins+Num_Chirp3_2_bins);
            Chirp3_3_spike_density_mat  = Chirp3_spike_density_mat(:,Num_Chirp3_1_bins+Num_Chirp3_2_bins+1:end);
            
        end
        
        if p.Dim_red_meth == 1 % PCA global
            
            [Chirp3_coeff,Chirp3_score,Chirp3_latent,Chirp3_tsquared,Chirp3_explained,Chirp3_mu] = pca(Chirp3_spike_density_mat); %  Chirp3_spike_density_mat'
            % Rows of X correspond to observations and columns correspond to variables.
            % needs data matrix in form: Rows = RFs, columns = times.
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Chirp3_explained_cumsum  = cumsum(Chirp3_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Chirp3_PC_Num = find(Chirp3_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Chirp3_PC_Num = p.PCA_comp_num;
                
            end
            
            All_Scores     = [All_Scores,Chirp3_score(:,1:Chirp3_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Chirp3_PC_Num;
            
            disp(sprintf('Num Chirp3 PCs = %i',Chirp3_PC_Num));
            
            Chirp3_pca_fig = figure;
            imagesc(Chirp3_ksdensity_grid,[],abs(Chirp3_coeff(:,1:Chirp3_PC_Num)'));%Chirp3_PC_Num, 10
            colormap(Chirp3_pca_fig,gray(256));
            colorbar;
            Chirp3_ks_v1 = vline(cl_var.Chirp3_trig_times_vec(2),'r');
            Chirp3_ks_v2 = vline(cl_var.Chirp3_trig_times_vec(3),'r');
            title('Chirp3 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(Chirp3_ks_v1,'LineWidth',1.5);
            set(Chirp3_ks_v2,'LineWidth',1.5);
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 2 % PCA segmented
            
            [Chirp3_1_coeff,Chirp3_1_score,Chirp3_1_latent,Chirp3_1_tsquared,Chirp3_1_explained,Chirp3_1_mu] = pca(Chirp3_1_spike_density_mat);
            [Chirp3_2_coeff,Chirp3_2_score,Chirp3_2_latent,Chirp3_2_tsquared,Chirp3_2_explained,Chirp3_2_mu] = pca(Chirp3_2_spike_density_mat);
            [Chirp3_3_coeff,Chirp3_3_score,Chirp3_3_latent,Chirp3_3_tsquared,Chirp3_3_explained,Chirp3_3_mu] = pca(Chirp3_3_spike_density_mat);
            
            if p.PCA_thresh_type == 1     % percentage variance explained
                
                % Calculate cumulative explained variances
                Chirp3_1_explained_cumsum  = cumsum(Chirp3_1_explained);
                Chirp3_2_explained_cumsum  = cumsum(Chirp3_2_explained);
                Chirp3_3_explained_cumsum  = cumsum(Chirp3_3_explained);
                
                % Find minimum number of components required to explain >PCA_ExpVar_Thresh variance in each case.
                Chirp3_1_PC_Num  = find(Chirp3_1_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Chirp3_2_PC_Num  = find(Chirp3_2_explained_cumsum>p.PCA_ExpVar_thresh,1);
                Chirp3_3_PC_Num  = find(Chirp3_3_explained_cumsum>p.PCA_ExpVar_thresh,1);
                
            else % p.PCA_thresh_type == 2 % number of components
                
                Chirp3_1_PC_Num  = p.PCA_comp_num;
                Chirp3_2_PC_Num  = p.PCA_comp_num;
                Chirp3_3_PC_Num  = p.PCA_comp_num;
                
            end
            
            Chirp3_PC_Num = Chirp3_1_PC_Num + Chirp3_2_PC_Num + Chirp3_3_PC_Num;
            
            All_Scores     = [All_Scores,Chirp3_1_score(:,1:Chirp3_1_PC_Num),Chirp3_2_score(:,1:Chirp3_2_PC_Num),Chirp3_3_score(:,1:Chirp3_3_PC_Num)];
            Total_Num_Comp = Total_Num_Comp + Chirp3_PC_Num;
            
            disp(sprintf('Num Chirp3 1 PCs = %i',Chirp3_1_PC_Num));
            disp(sprintf('Num Chirp3 2 PCs = %i',Chirp3_2_PC_Num));
            disp(sprintf('Num Chirp3 3 PCs = %i',Chirp3_3_PC_Num));
            disp(sprintf('Total Num Chirp3 PCs = %i',Chirp3_PC_Num));
            
            Chirp3_pca_fig = figure;
            subplot(1,3,1);
            imagesc(Chirp3_ksdensity_grid(1:Num_Chirp3_1_bins),[],abs(Chirp3_1_coeff(:,1:Chirp3_1_PC_Num)'));
            colormap(Chirp3_pca_fig,gray(256));
            colorbar;
            title('Chirp3 1 - PCA coeff'); xlabel('time (sec)'); ylabel('component');
            set(gca,'FontSize',12);
            subplot(1,3,2);
            imagesc(Chirp3_ksdensity_grid(Num_Chirp3_1_bins+1:Num_Chirp3_1_bins+Num_Chirp3_2_bins),[],abs(Chirp3_2_coeff(:,1:Chirp3_2_PC_Num)'));
            colormap(Chirp3_pca_fig,gray(256));
            title('Chirp3 2 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            subplot(1,3,3);
            imagesc(Chirp3_ksdensity_grid(Num_Chirp3_1_bins+Num_Chirp3_2_bins+1:end),[],abs(Chirp3_3_coeff(:,1:Chirp3_3_PC_Num)'));
            colormap(Chirp3_pca_fig,gray(256));
            title('Chirp3 3 - PCA coeff'); xlabel('time (sec)');
            set(gca,'FontSize',12);
            set(gcf,'color','w');
            
        elseif p.Dim_red_meth == 3 % sparse PCA global
            
        else % p.Dim_red_meth == 4 % sparse PCA segmented
            
        end
        
    end
    
end



%%% Scale scores
if p.Scale_Scores == 1
    for i = 1:size(All_Scores,2)
        All_Scores_min_loop = min(All_Scores(:,i));
        All_Scores_max_loop = max(All_Scores(:,i));
        All_Scores(:,i) = (All_Scores(:,i) - All_Scores_min_loop)/(All_Scores_max_loop - All_Scores_min_loop);
    end
end



%%% PC Components
% Used to find column of all scores corresponding to CNoise: Full RF Dominant Axis Angle
% so that can change NaNs to -1s.
% And similarly for CNoise: Full RF Ellipticity (which is NaN where RF Size = 0)
% 1. FFF_PC_Num
% 2. Chirp_PC_Num
% 3. FFF_Noise_PC_Num
% 4. Gratings_400px_PC_Num
% 5. 1 for CNoise: Full RF Size
% 6. 1 for CNoise: Full RF Ellipticity 
% 7. 1 for CNoise: Full RF Dominant Axis Angle
% 8. FFF2_PC_Num
% 9. Chirp2_PC_Num
% 10. SSub_PC_Num
% 11. CSteps_PC_Num
% 12. Chirp3_PC_Num
p.PC_Num_Partial = 0;
p.PC_Num_Partial_2 = 0;
if p.Obj_clust_vec(7) == 1
    if p.Obj_clust_vec(1) == 1 % FFF
        p.PC_Num_Partial   = p.PC_Num_Partial   + FFF_PC_Num;
        p.PC_Num_Partial_2 = p.PC_Num_Partial_2 + FFF_PC_Num;
    end
    if p.Obj_clust_vec(2) == 1 % Chirp
        p.PC_Num_Partial   = p.PC_Num_Partial   + Chirp_PC_Num;
        p.PC_Num_Partial_2 = p.PC_Num_Partial_2 + Chirp_PC_Num;
    end
    if p.Obj_clust_vec(3) == 1 % FFF_Noise
        p.PC_Num_Partial   = p.PC_Num_Partial   + FFF_Noise_PC_Num;
        p.PC_Num_Partial_2 = p.PC_Num_Partial_2 + FFF_Noise_PC_Num;
    end
    if p.Obj_clust_vec(4) == 1 % Gratings_400px
        p.PC_Num_Partial   = p.PC_Num_Partial   + Gratings_400px_PC_Num;
        p.PC_Num_Partial_2 = p.PC_Num_Partial_2 + Gratings_400px_PC_Num;
    end
    if p.Obj_clust_vec(5) == 1 % CNoise: Full RF Size
        p.PC_Num_Partial   = p.PC_Num_Partial   + 1;
        p.PC_Num_Partial_2 = p.PC_Num_Partial_2 + 1;
    end
    if p.Obj_clust_vec(6) == 1 % CNoise: Full RF Ellipticity
        p.PC_Num_Partial   = p.PC_Num_Partial   + 1;
        p.PC_Num_Partial_2 = p.PC_Num_Partial_2 + 1;
    end
    if p.Obj_clust_vec(7) == 1 % CNoise: Full RF Dominant Axis Angle
        p.PC_Num_Partial = p.PC_Num_Partial + 1;
    end
end




%% Cluster

if p.Clust_Meth == 1     % Gaussian Mixture Models
    [gm,aic,bic,converged,allConverge,clust_vec,Num_clust] = GMM_Clust_fn_3(All_Scores,p); % PAR Mod 29,11,2021 (was GMM_Clust_fn_2)
else % p.Clust_Meth == 2 % Hierarchical
    [clust_vec,Num_clust,cophenetic_corr_coeff] = Hier_Clust_fn_2(All_Scores,p);
end

%save('data_DP5_RFC6_19_11_2021_2.mat','-v7.3');
%save('data_DP5_RFC6_29_11_2021_1.mat','-v7.3');


%% Eliminate clusters with fewer than a critical number of ROIs

Num_Cluster_Vec = zeros(Num_clust,1); % Number of members in each cluster

for i = 1:Num_clust
    Num_Cluster_Vec(i) = sum(clust_vec==i);
end

TrueClusIndex_Vec = find(Num_Cluster_Vec>=p.MinClusROI_thresh); % Indicies of clusters with greater than/equal to the critical number of members
Num_clust_True    = length(TrueClusIndex_Vec);                  % The number of clusters with greater than/equal to the critical number of members


%% Calculate Scaling Factors for Plots

%%% FFF
if p.Obj_plot_vec(3) == 1
    FFF_Clus_mean_mat = zeros(Num_clust_True,FFF_length_ksdensity_grid);
    FFF_Clus_std_mat  = zeros(Num_clust_True,FFF_length_ksdensity_grid);
    for i = 1:Num_clust_True
        FFF_Mat_Clusi          = FFF_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        FFF_Clus_mean_mat(i,:) = mean(FFF_Mat_Clusi,1);
        FFF_Clus_std_mat(i,:)  = std(FFF_Mat_Clusi,0,1);
    end
    % Use these to scale the FFF traces
    FFF_trace_LwrBd = min(FFF_Clus_mean_mat-FFF_Clus_std_mat,[],'all');
    FFF_trace_UprBd = max(FFF_Clus_mean_mat+FFF_Clus_std_mat,[],'all');
    % Use these to scale the FFF heat maps
    FFF_min_resp = min(FFF_spike_density_mat,[],'all');
    FFF_max_resp = max(FFF_spike_density_mat,[],'all');
end

%%% Chirp
if p.Obj_plot_vec(4) == 1
    Chirp_Clus_mean_mat = zeros(Num_clust_True,Chirp_length_ksdensity_grid);
    Chirp_Clus_std_mat  = zeros(Num_clust_True,Chirp_length_ksdensity_grid);
    for i = 1:Num_clust_True
        Chirp_Mat_Clusi          = Chirp_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        Chirp_Clus_mean_mat(i,:) = mean(Chirp_Mat_Clusi,1);
        Chirp_Clus_std_mat(i,:)  = std(Chirp_Mat_Clusi,0,1);
    end
    % Use these to scale the Chirp traces
    Chirp_trace_LwrBd = min(Chirp_Clus_mean_mat-Chirp_Clus_std_mat,[],'all');
    Chirp_trace_UprBd = max(Chirp_Clus_mean_mat+Chirp_Clus_std_mat,[],'all');
    % Use these to scale the Chirp heat maps
    Chirp_min_resp = min(Chirp_spike_density_mat,[],'all');
    Chirp_max_resp = max(Chirp_spike_density_mat,[],'all');
end

%%% FFF Noise
if p.Obj_plot_vec(5) == 1
    FFF_Noise_Clus_mean_mat = zeros(Num_clust_True,p.FFF_Noise_Spectral_Dim*p.Num_STE_bins);
   %FFF_Noise_Clus_std_mat  = zeros(Num_clust_True,p.FFF_Noise_Spectral_Dim*p.Num_STE_bins);
    for i = 1:Num_clust_True
        %FFF_Noise_Mat_Clusi          = cl_var.FFF_Noise_Full_STA(clust_vec==TrueClusIndex_Vec(i),:);
        FFF_Noise_Mat_Clusi          = cl_var.FFF_Noise_Full_STA_scaled(clust_vec==TrueClusIndex_Vec(i),:); % Mod 14,06,2021
        FFF_Noise_Clus_mean_mat(i,:) = nanmean(FFF_Noise_Mat_Clusi,1);  % mean
       %FFF_Noise_Clus_std_mat(i,:)  = nanstd(FFF_Noise_Mat_Clusi,0,1); % std
    end
    % Use these to scale the FFF_Noise traces
   %FFF_Noise_trace_LwrBd = min([FFF_Noise_Clus_mean_mat-FFF_Noise_Clus_std_mat,FFF_Noise_Clus_mean_mat-FFF_Noise_Clus_std_mat],[],'all');
   %FFF_Noise_trace_UprBd = max([FFF_Noise_Clus_mean_mat+FFF_Noise_Clus_std_mat,FFF_Noise_Clus_mean_mat+FFF_Noise_Clus_std_mat],[],'all');
    FFF_Noise_trace_LwrBd = min(FFF_Noise_Clus_mean_mat,[],'all');
    FFF_Noise_trace_UprBd = max(FFF_Noise_Clus_mean_mat,[],'all');
    % Use these to scale the FFF_Noise heat maps
    %FFF_Noise_min_resp = min(cl_var.FFF_Noise_Full_STA,[],'all');
    %FFF_Noise_max_resp = max(cl_var.FFF_Noise_Full_STA,[],'all');
    FFF_Noise_min_resp = min(cl_var.FFF_Noise_Full_STA_scaled,[],'all');% Mod 14,06,2021
    FFF_Noise_max_resp = max(cl_var.FFF_Noise_Full_STA_scaled,[],'all');% Mod 14,06,2021
end

%%% Gratings 400px
if p.Obj_plot_vec(6) == 1
    Gratings_400px_Clus_mean_mat = zeros(Num_clust_True,Gratings_400px_length_ksdensity_grid);
    Gratings_400px_Clus_std_mat  = zeros(Num_clust_True,Gratings_400px_length_ksdensity_grid);
    for i = 1:Num_clust_True
        Gratings_400px_Mat_Clusi          = Gratings_400px_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        Gratings_400px_Clus_mean_mat(i,:) = mean(Gratings_400px_Mat_Clusi,1);
        Gratings_400px_Clus_std_mat(i,:)  = std(Gratings_400px_Mat_Clusi,0,1);
    end
    % Use these to scale the Gratings 400px traces
    Gratings_400px_trace_LwrBd = min(Gratings_400px_Clus_mean_mat-Gratings_400px_Clus_std_mat,[],'all');
    Gratings_400px_trace_UprBd = max(Gratings_400px_Clus_mean_mat+Gratings_400px_Clus_std_mat,[],'all');
    % Use these to scale the Gratings 400px heat maps
    Gratings_400px_min_resp = min(Gratings_400px_spike_density_mat,[],'all');
    Gratings_400px_max_resp = max(Gratings_400px_spike_density_mat,[],'all');
end

%%% CNoise: Full RF Size
if p.Obj_plot_vec(7) == 1
    RF_Size_hist_Clus_MaxValue_vec = zeros(Num_clust_True,1);
    min_Full_RF_Size_hist_BinLim   = min(cl_var.Full_RF_Size_hist_BinLimits_mat(:,1));
    max_Full_RF_Size_hist_BinLim   = max(cl_var.Full_RF_Size_hist_BinLimits_mat(:,2));
    min_Full_RF_Size_hist_BinWidth = min(cl_var.Full_RF_Size_hist_BinWidth_vec);
    RF_Size_hist_NumBins_Final     = ceil((max_Full_RF_Size_hist_BinLim-min_Full_RF_Size_hist_BinLim)/min_Full_RF_Size_hist_BinWidth);
    RF_Size_hist_BinEdges_Final    = linspace(min_Full_RF_Size_hist_BinLim,max_Full_RF_Size_hist_BinLim,RF_Size_hist_NumBins_Final+1);
    for i = 1:Num_clust_True
        Full_RF_Size_vec_Clusi = cl_var.Full_RF_Size_vec(clust_vec==TrueClusIndex_Vec(i));
        %hist_loop = histogram(Full_RF_Size_vec_Clusi,cl_var.RF_Size_hist_BinEdges);
        [values_loop,~] = histcounts(Full_RF_Size_vec_Clusi,RF_Size_hist_BinEdges_Final); % cl_var.RF_Size_hist_BinEdges
        RF_Size_hist_Clus_MaxValue_vec(i) = max(values_loop); % max(hist_loop.Values)
    end
    RF_Size_hist_MaxValue = max(RF_Size_hist_Clus_MaxValue_vec);
end

%%% CNoise: Full RF Ellipticity
if p.Obj_plot_vec(8) == 1
    RF_Ellip_hist_Clus_MaxValue_vec = zeros(Num_clust_True,1);
    min_Full_RF_Ellip_hist_BinLim   = min(cl_var.Full_RF_Ellipticity_hist_BinLimits_mat(:,1));
    max_Full_RF_Ellip_hist_BinLim   = max(cl_var.Full_RF_Ellipticity_hist_BinLimits_mat(:,2));
    min_Full_RF_Ellip_hist_BinWidth = min(cl_var.Full_RF_Ellipticity_hist_BinWidth_vec);
    RF_Ellip_hist_NumBins_Final     = ceil((max_Full_RF_Ellip_hist_BinLim-min_Full_RF_Ellip_hist_BinLim)/min_Full_RF_Ellip_hist_BinWidth);
    RF_Ellip_hist_BinEdges_Final    = linspace(min_Full_RF_Ellip_hist_BinLim,max_Full_RF_Ellip_hist_BinLim,RF_Ellip_hist_NumBins_Final+1);
    for i = 1:Num_clust_True
        Full_RF_Ellip_vec_Clusi = cl_var.Full_RF_Ellipticity_vec(clust_vec==TrueClusIndex_Vec(i));
        %hist_loop = histogram(Full_RF_Ellip_vec_Clusi,cl_var.RF_Ellip_hist_BinEdges);
        [values_loop,~] = histcounts(Full_RF_Ellip_vec_Clusi,RF_Ellip_hist_BinEdges_Final); % cl_var.RF_Ellip_hist_BinEdges
        RF_Ellip_hist_Clus_MaxValue_vec(i) = max(values_loop); % max(hist_loop.Values
    end
    RF_Ellip_hist_MaxValue = max(RF_Ellip_hist_Clus_MaxValue_vec);
end

%%% CNoise: Full RF Dominant Axis Angle --> This isn't needed.
% if p.Obj_plot_vec(9) == 1
%     RF_Dom_Ax_Ang_hist_Clus_MaxValue_vec = zeros(Num_clust_True,1);
%     min_Full_RF_Dom_Ax_Ang_hist_BinLim   = min(cl_var.Full_RF_Dom_Ax_Ang_hist_BinLimits_mat(:,1));
%     max_Full_RF_Dom_Ax_Ang_hist_BinLim   = max(cl_var.Full_RF_Dom_Ax_Ang_hist_BinLimits_mat(:,2));
%     min_Full_RF_Dom_Ax_Ang_hist_BinWidth = min(cl_var.Full_RF_Dom_Ax_Ang_hist_BinWidth_vec);
%     RF_Dom_Ax_Ang_hist_NumBins_Final     = ceil((max_Full_RF_Dom_Ax_Ang_hist_BinLim-min_Full_RF_Dom_Ax_Ang_hist_BinLim)/min_Full_RF_Dom_Ax_Ang_hist_BinWidth);
%     RF_Dom_Ax_Ang_hist_BinEdges_Final    = linspace(min_Full_RF_Dom_Ax_Ang_hist_BinLim,max_Full_RF_Dom_Ax_Ang_hist_BinLim,RF_Dom_Ax_Ang_hist_NumBins_Final+1);
%     for i = 1:Num_clust_True
%         Full_RF_Dom_Ax_Ang_vec_Clusi = cl_var.Full_RF_Dom_Ax_Ang_vec(clust_vec==TrueClusIndex_Vec(i));
%         [values_loop,~] = histcounts(Full_RF_Dom_Ax_Ang_vec_Clusi,RF_Dom_Ax_Ang_hist_BinEdges_Final);
%         RF_Dom_Ax_Ang_hist_Clus_MaxValue_vec(i) = max(values_loop);
%     end
%     RF_Dom_Ax_Ang_hist_MaxValue = max(RF_Dom_Ax_Ang_hist_Clus_MaxValue_vec);
% end


%%% FFF2
if p.Obj_plot_vec(10) == 1
    FFF2_Clus_mean_mat = zeros(Num_clust_True,FFF2_length_ksdensity_grid);
    FFF2_Clus_std_mat  = zeros(Num_clust_True,FFF2_length_ksdensity_grid);
    for i = 1:Num_clust_True
        FFF2_Mat_Clusi          = FFF2_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        FFF2_Clus_mean_mat(i,:) = mean(FFF2_Mat_Clusi,1);
        FFF2_Clus_std_mat(i,:)  = std(FFF2_Mat_Clusi,0,1);
    end
    % Use these to scale the FFF2 traces
    FFF2_trace_LwrBd = min(FFF2_Clus_mean_mat-FFF2_Clus_std_mat,[],'all');
    FFF2_trace_UprBd = max(FFF2_Clus_mean_mat+FFF2_Clus_std_mat,[],'all');
    % Use these to scale the FFF2 heat maps
    FFF2_min_resp = min(FFF2_spike_density_mat,[],'all');
    FFF2_max_resp = max(FFF2_spike_density_mat,[],'all');
end

%%% Chirp2
if p.Obj_plot_vec(11) == 1
    Chirp2_Clus_mean_mat = zeros(Num_clust_True,Chirp2_length_ksdensity_grid);
    Chirp2_Clus_std_mat  = zeros(Num_clust_True,Chirp2_length_ksdensity_grid);
    for i = 1:Num_clust_True
        Chirp2_Mat_Clusi          = Chirp2_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        Chirp2_Clus_mean_mat(i,:) = mean(Chirp2_Mat_Clusi,1);
        Chirp2_Clus_std_mat(i,:)  = std(Chirp2_Mat_Clusi,0,1);
    end
    % Use these to scale the Chirp2 traces
    Chirp2_trace_LwrBd = min(Chirp2_Clus_mean_mat-Chirp2_Clus_std_mat,[],'all');
    Chirp2_trace_UprBd = max(Chirp2_Clus_mean_mat+Chirp2_Clus_std_mat,[],'all');
    % Use these to scale the Chirp2 heat maps
    Chirp2_min_resp = min(Chirp2_spike_density_mat,[],'all');
    Chirp2_max_resp = max(Chirp2_spike_density_mat,[],'all');
end

%%% SSub
if p.Obj_plot_vec(12) == 1
    SSub_Clus_mean_mat = zeros(Num_clust_True,SSub_length_ksdensity_grid);
    SSub_Clus_std_mat  = zeros(Num_clust_True,SSub_length_ksdensity_grid);
    for i = 1:Num_clust_True
        SSub_Mat_Clusi          = SSub_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        SSub_Clus_mean_mat(i,:) = mean(SSub_Mat_Clusi,1);
        SSub_Clus_std_mat(i,:)  = std(SSub_Mat_Clusi,0,1);
    end
    % Use these to scale the SSub traces
    SSub_trace_LwrBd = min(SSub_Clus_mean_mat-SSub_Clus_std_mat,[],'all');
    SSub_trace_UprBd = max(SSub_Clus_mean_mat+SSub_Clus_std_mat,[],'all');
    % Use these to scale the SSub heat maps
    SSub_min_resp = min(SSub_spike_density_mat,[],'all');
    SSub_max_resp = max(SSub_spike_density_mat,[],'all');
end


%%% CSteps
if p.Obj_plot_vec(13) == 1
    CSteps_Clus_mean_mat = zeros(Num_clust_True,CSteps_length_ksdensity_grid);
    CSteps_Clus_std_mat  = zeros(Num_clust_True,CSteps_length_ksdensity_grid);
    for i = 1:Num_clust_True
        CSteps_Mat_Clusi          = CSteps_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        CSteps_Clus_mean_mat(i,:) = mean(CSteps_Mat_Clusi,1);
        CSteps_Clus_std_mat(i,:)  = std(CSteps_Mat_Clusi,0,1);
    end
    % Use these to scale the CSteps traces
    CSteps_trace_LwrBd = min(CSteps_Clus_mean_mat-CSteps_Clus_std_mat,[],'all');
    CSteps_trace_UprBd = max(CSteps_Clus_mean_mat+CSteps_Clus_std_mat,[],'all');
    % Use these to scale the CSteps heat maps
    CSteps_min_resp = min(CSteps_spike_density_mat,[],'all');
    CSteps_max_resp = max(CSteps_spike_density_mat,[],'all');
end


%%% Chirp3
if p.Obj_plot_vec(14) == 1
    Chirp3_Clus_mean_mat = zeros(Num_clust_True,Chirp3_length_ksdensity_grid);
    Chirp3_Clus_std_mat  = zeros(Num_clust_True,Chirp3_length_ksdensity_grid);
    for i = 1:Num_clust_True
        Chirp3_Mat_Clusi          = Chirp3_spike_density_mat(clust_vec==TrueClusIndex_Vec(i),:);
        Chirp3_Clus_mean_mat(i,:) = mean(Chirp3_Mat_Clusi,1);
        Chirp3_Clus_std_mat(i,:)  = std(Chirp3_Mat_Clusi,0,1);
    end
    % Use these to scale the Chirp3 traces
    Chirp3_trace_LwrBd = min(Chirp3_Clus_mean_mat-Chirp3_Clus_std_mat,[],'all');
    Chirp3_trace_UprBd = max(Chirp3_Clus_mean_mat+Chirp3_Clus_std_mat,[],'all');
    % Use these to scale the Chirp3 heat maps
    Chirp3_min_resp = min(Chirp3_spike_density_mat,[],'all');
    Chirp3_max_resp = max(Chirp3_spike_density_mat,[],'all');
end


%% Reorder Clusters

%clust_vec		        Cluster labels for each cell
%Num_Cluster_Vec		Number of members in each cluster
%TrueClusIndex_Vec		Indicies of clusters with greater than/equal to the critical number of members
%Num_clust_True 		The number of clusters with greater than/equal to the critical number of members
%Cell_per_Cluster_mat = NaN(Num_clust_True,p.Num_data_sets); Number of Cells from each Data Set in Each Cluster


% Reorder clusters?
% 1. no;
% 2. yes.
Reorder_Choice = 2; % 2 usually

% Choose what to reorder by
if Reorder_Choice == 2
    % 1. Cluster size
    Reorder_Type = 1;
end

if Reorder_Choice == 2
    if Reorder_Type == 1
        Num_Cluster_Vec_True = Num_Cluster_Vec(TrueClusIndex_Vec);
        [Num_Cluster_Vec_True_ordered,Sorting_Vec] = sort(Num_Cluster_Vec_True,'descend');
        % Sorting Vec needs to be relabelled if we filtered out small
        % clusters as the cluster and index numbers in TrueClusIndex_Vec are
        % not the same in this case.
        Sorting_Vec_2 = zeros(Num_clust_True,1);
        for i = 1:Num_clust_True
            Sorting_Vec_2(i) = TrueClusIndex_Vec(Sorting_Vec(i));
        end
    end
    % Create Vector of Reordered Number of Cells in Each Cluster if don't
    % order by this originally
end

% Create Final Cluster Index Vec
if Reorder_Choice == 1
    FinalClusIndex_Vec = TrueClusIndex_Vec;
else %  Reorder_Choice == 2
    FinalClusIndex_Vec = Sorting_Vec_2;
end


%% Calculate Number of Cells from each Data Set in Each Cluster
% Moved here from above after added reordering command 22,09,2021

Cell_per_Cluster_mat = NaN(Num_clust_True,p.Num_data_sets);

for i = 1:Num_clust_True
    for j = 1:p.Num_data_sets
        %Exp_Ident_vec_Clusi       = cl_var.Exp_Ident_vec(clust_vec==TrueClusIndex_Vec(i));
        Exp_Ident_vec_Clusi       = cl_var.Exp_Ident_vec(clust_vec==FinalClusIndex_Vec(i));
        Cell_per_Cluster_mat(i,j) = sum(Exp_Ident_vec_Clusi==j);
    end
end


%%% Data set bar chart
if p.Obj_plot_vec(14) == 1
%     Cell_per_Cluster_hist_Clus_MaxValue_vec = zeros(Num_clust_True,1);
%     for i = 1:Num_clust_True
%         Cell_per_Cluster_vec_Clusi = cl_var.Full_RF_Ellipticity_vec(clust_vec==TrueClusIndex_Vec(i));
%         [values_loop,~] = histcounts(Cell_per_Cluster_vec_Clusi,cl_var.RF_Ellip_hist_BinEdges);
%         Cell_per_Cluster_hist_Clus_MaxValue_vec(i) = max(values_loop);
%     end
%     Cell_per_Cluster_hist_MaxValue = max(Cell_per_Cluster_hist_Clus_MaxValue_vec);
    Cell_per_Cluster_bar_MaxValue = max(Cell_per_Cluster_mat,[],'all');
end

%% Plot Clustering Results

% Distribute clusters between separate figures when there are too many
% clusters for a single figure.
First_Entry = 1;                            % 1
Last_Entry  = 10;               % Num_clust_True
Fig_Rows    = Last_Entry - First_Entry + 1; % Num_clust_True

%%%
% Plotting Style Commands
v_line_width  = 1;   % vertical line width
v_line_transp = 0.5; % vertical line transparency
%%%

% Choose subplot or subaxis
% 1. subplot;
% 2. subaxis.
Plot_Axis_Choice = 2;

% Remove x-tick labels for all but last row?
% 1. no;
% 2. yes.
xtick_Choice = 2;

% Normalise data set bar plots
% 1. no;
% 2. yes.
Norm_Data_Bar_Choice = 2;

if Norm_Data_Bar_Choice == 2
    Cell_per_Cluster_mat_norm = NaN(Num_clust_True,p.Num_data_sets);
    for i = 1:p.Num_data_sets
        Cell_per_Cluster_mat_norm(:,i) = Cell_per_Cluster_mat(:,i)/cl_var.Num_Cell_Per_Exp_vec(i);
    end
end

if Plot_Axis_Choice == 2
    if xtick_Choice == 1
        Vert_Spacing        = 0.5*1e-1;  % 0.5*1e-1 (0.5*1e-1 works for 9 rows x 10 col)
    else
        Vert_Spacing        = 0.25*1e-1; % 0.5*1e-1 (0.25*1e-1 works for 9 rows x 10 col)
    end
    Horiz_Spacing       = 0.25*1e-1; % 0.5*1e-1 (0.25*1e-1 works for 9 rows x 10 col)
    Horiz_Spacing_Polar = 1e-2;      % 1e-2     (0 works for 9 rows x 10 col, need more spacing for fewer rows unless remove y-tick labels)
    Left_Margin         = 2.5*1e-2;  % 1e-1     (2.5*1e-2 works for 9 rows x 10 col)
    Right_Margin        = 2.5*1e-2;  % 1e-1     (2.5*1e-2 works for 9 rows x 10 col)
    Top_Margin          = 2.5*1e-2;  % 1e-1     (2.5*1e-2 works for 9 rows x 10 col)
    Bottom_Margin       = 1e-1;      % 1e-1     (5*1e-2 works for 9 rows x 10 col) % PAR Mod 22,09,2021 was 5*1e-2 made bigger for new data labels % PAR Mod 04,11,2021 was 7.5*1e-2 made bigger for new data labels
end
% Defaults
% 'SpacingVertical',0.05,'SpacingHorizontal',0.05,...
% 'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0,...
% 'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1,

Num_Columns = p.Num_plots;

Col_index_vec = [0,cumsum(p.Obj_plot_vec.*p.Plot_per_obj_vec)];
Col_index_vec = Col_index_vec(1:end-1)+1;
    
figure;

for i = First_Entry:Last_Entry
    
    %%% Cell positions
    if p.Obj_plot_vec(1) == 1
        Full_RF_Pos_Mat_Clusi = cl_var.Full_RF_Pos_mat(clust_vec==FinalClusIndex_Vec(i),:);
        
        % Cell position scatter plot
        if Plot_Axis_Choice == 1
            ax1 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(1));
        else % Plot_Axis_Choice == 2
            ax1 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(1),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        scatter(Full_RF_Pos_Mat_Clusi(:,1),Full_RF_Pos_Mat_Clusi(:,2),'x','LineWidth',1.5);
        axis equal;
        set(gca,'Ydir','reverse');
        xlim([0.5 (p.Num_rows + 0.5)]);
        ylim([0.5 (p.Num_cols + 0.5)]);
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('x');
        end
        %ylabel('y');
        if i==First_Entry
            title(p.Plot_obj_name_vec{1}); % title('cell positions');
        end
        set(gca,'FontSize',12);
    end
    
    %%% FFF
    if p.Obj_plot_vec(3) == 1
        FFF_Mat_Clusi = FFF_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        
        % FFF trace
        if Plot_Axis_Choice == 1
            ax3 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(3));
        else % Plot_Axis_Choice == 2
            ax3 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(3),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(FFF_ksdensity_grid,FFF_Clus_mean_mat(i,:),FFF_Clus_std_mat(i,:),{'k-','markerfacecolor','k'}); hold on;
        xlim([0 cl_var.FFF_stim_end_time]);
        %ylim([FFF_trace_LwrBd FFF_trace_UprBd]);
        FFF_ks_v1 = vline(cl_var.FFF_trig_times_vec(2),'r'); FFF_ks_v1.Color = [FFF_ks_v1.Color v_line_transp];
        FFF_ks_v2 = vline(cl_var.FFF_trig_times_vec(3),'r'); FFF_ks_v2.Color = [FFF_ks_v2.Color v_line_transp];
        FFF_ks_v3 = vline(cl_var.FFF_trig_times_vec(4),'r'); FFF_ks_v3.Color = [FFF_ks_v3.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{3}); % title('FFF');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(FFF_ks_v1,'LineWidth',v_line_width);
        set(FFF_ks_v2,'LineWidth',v_line_width);
        set(FFF_ks_v3,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
        
        % FFF heatmap
        if Plot_Axis_Choice == 1
            ax4 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(3)+1);
        else % Plot_Axis_Choice == 2
            ax4 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(3)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        imagesc(FFF_ksdensity_grid,[],FFF_Mat_Clusi,[FFF_min_resp FFF_max_resp]); colormap(ax4,gray(256)); hold on;
        FFF_ks_v1 = vline(cl_var.FFF_trig_times_vec(2),'r'); FFF_ks_v1.Color = [FFF_ks_v1.Color v_line_transp];
        FFF_ks_v2 = vline(cl_var.FFF_trig_times_vec(3),'r'); FFF_ks_v2.Color = [FFF_ks_v2.Color v_line_transp];
        FFF_ks_v3 = vline(cl_var.FFF_trig_times_vec(4),'r'); FFF_ks_v3.Color = [FFF_ks_v3.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{3}); % title('FFF');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(FFF_ks_v1,'LineWidth',v_line_width);
        set(FFF_ks_v2,'LineWidth',v_line_width);
        set(FFF_ks_v3,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
    end
    
    %%% Chirp
    if p.Obj_plot_vec(4) == 1
        Chirp_Mat_Clusi = Chirp_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        
        % Chirp trace
        if Plot_Axis_Choice == 1
            ax5 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(4));
        else % Plot_Axis_Choice == 2
            ax5 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(4),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(Chirp_ksdensity_grid,Chirp_Clus_mean_mat(i,:),Chirp_Clus_std_mat(i,:),{'k-','markerfacecolor','k'}); hold on;
        xlim([0 cl_var.Chirp_stim_end_time]);
        %ylim([Chirp_trace_LwrBd Chirp_trace_UprBd]);
        Chirp_ks_v1 = vline(cl_var.Chirp_trig_times_vec(2),'r'); Chirp_ks_v1.Color = [Chirp_ks_v1.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{4}); % title('Chirp');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Chirp_ks_v1,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
        
        % Chirp heatmap
        if Plot_Axis_Choice == 1
            ax6 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(4)+1);
        else % Plot_Axis_Choice == 2
            ax6 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(4)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        imagesc(Chirp_ksdensity_grid,[],Chirp_Mat_Clusi,[Chirp_min_resp Chirp_max_resp]); colormap(ax6,gray(256)); hold on;
        Chirp_ks_v1 = vline(cl_var.Chirp_trig_times_vec(2),'r'); Chirp_ks_v1.Color = [Chirp_ks_v1.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{4}); % title('Chirp');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Chirp_ks_v1,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
    end
    
    %%% FFF Noise
    if p.Obj_plot_vec(5) == 1
        %FFF_Noise_Mat_Clusi = cl_var.FFF_Noise_Full_STA(clust_vec==FinalClusIndex_Vec(i),:);
        FFF_Noise_Mat_Clusi      = cl_var.FFF_Noise_Full_STA_scaled(clust_vec==FinalClusIndex_Vec(i),:); % Mod 14,06,2021
        FFF_Noise_mean_mat_Clusi = FFF_Noise_Clus_mean_mat(FinalClusIndex_Vec(i),:); % PAR Mod 24,11,2021
        
        % FFF_Noise trace
        if Plot_Axis_Choice == 1
            ax7 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(5));
        else % Plot_Axis_Choice == 2
            ax7 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(5),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        plot(p.stim_timesample_vec,FFF_Noise_mean_mat_Clusi(1:p.Num_STE_bins),'r','LineWidth',1.5); hold on;           % PAR Mod 24,11,2021 % plot(p.stim_timesample_vec,FFF_Noise_Clus_mean_mat(i,1:p.Num_STE_bins),'r','LineWidth',1.5); hold on; 
        plot(p.stim_timesample_vec,FFF_Noise_mean_mat_Clusi(p.Num_STE_bins+1:2*p.Num_STE_bins),'g','LineWidth',1.5);   % PAR Mod 24,11,2021 plot(p.stim_timesample_vec,FFF_Noise_Clus_mean_mat(i,p.Num_STE_bins+1:2*p.Num_STE_bins),'g','LineWidth',1.5);
        plot(p.stim_timesample_vec,FFF_Noise_mean_mat_Clusi(2*p.Num_STE_bins+1:3*p.Num_STE_bins),'b','LineWidth',1.5); % PAR Mod 24,11,2021 plot(p.stim_timesample_vec,FFF_Noise_Clus_mean_mat(i,2*p.Num_STE_bins+1:3*p.Num_STE_bins),'b','LineWidth',1.5);
        plot(p.stim_timesample_vec,FFF_Noise_mean_mat_Clusi(3*p.Num_STE_bins+1:4*p.Num_STE_bins),'m','LineWidth',1.5); % PAR Mod 24,11,2021 plot(p.stim_timesample_vec,FFF_Noise_Clus_mean_mat(i,3*p.Num_STE_bins+1:4*p.Num_STE_bins),'m','LineWidth',1.5);
        xlim([p.stim_timesample_vec(1) p.stim_timesample_vec(end)]); % xlim([p.stim_timesample_vec(1) 0]);
        %ylim([FFF_Noise_trace_LwrBd FFF_Noise_trace_UprBd]);
        if i==First_Entry
            title(p.Plot_obj_name_vec{5}); % title('FFF Noise');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(gca,'FontSize',12);
        
        % FFF_Noise heatmap
        if Plot_Axis_Choice == 1
            ax8 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(5)+1);
        else % Plot_Axis_Choice == 2
            ax8 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(5)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        imagesc(cl_var.FFF_Noise_Full_t_vec+0.5*cl_var.FFF_Noise_Full_t_vec_int,[],FFF_Noise_Mat_Clusi,[FFF_Noise_min_resp FFF_Noise_max_resp]); colormap(ax8,bluewhitered(256)); hold on;
        FFF_Noise_STA_v1 = vline(0.25*(cl_var.FFF_Noise_Full_t_vec_end+cl_var.FFF_Noise_Full_t_vec_int),'k'); FFF_Noise_STA_v1.Color = [FFF_Noise_STA_v1.Color v_line_transp];
        FFF_Noise_STA_v2 = vline(0.5*(cl_var.FFF_Noise_Full_t_vec_end+cl_var.FFF_Noise_Full_t_vec_int),'k');  FFF_Noise_STA_v2.Color = [FFF_Noise_STA_v2.Color v_line_transp];
        FFF_Noise_STA_v3 = vline(0.75*(cl_var.FFF_Noise_Full_t_vec_end+cl_var.FFF_Noise_Full_t_vec_int),'k'); FFF_Noise_STA_v3.Color = [FFF_Noise_STA_v3.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{5}); % title('FFF Noise');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(FFF_Noise_STA_v1,'LineWidth',v_line_width);
        set(FFF_Noise_STA_v2,'LineWidth',v_line_width);
        set(FFF_Noise_STA_v3,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
    end
    
    %%% Gratings_400px
    if p.Obj_plot_vec(6) == 1
        Gratings_400px_Mat_Clusi = Gratings_400px_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        
        % Gratings 400px trace
        if Plot_Axis_Choice == 1
            ax9 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(6));
        else % Plot_Axis_Choice == 2
            ax9 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(6),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(Gratings_400px_ksdensity_grid,Gratings_400px_Clus_mean_mat(i,:),Gratings_400px_Clus_std_mat(i,:),{'k-','markerfacecolor','k'}); hold on;
        xlim([0 cl_var.Gratings_400px_stim_end_time]);
        ylim([Gratings_400px_trace_LwrBd Gratings_400px_trace_UprBd]);
        Gratings_400px_ks_v1 = vline(cl_var.Gratings_400px_trig_times_vec(2),'r'); Gratings_400px_ks_v1.Color = [Gratings_400px_ks_v1.Color v_line_transp];
        Gratings_400px_ks_v2 = vline(cl_var.Gratings_400px_trig_times_vec(3),'r'); Gratings_400px_ks_v2.Color = [Gratings_400px_ks_v2.Color v_line_transp];
        Gratings_400px_ks_v3 = vline(cl_var.Gratings_400px_trig_times_vec(4),'r'); Gratings_400px_ks_v3.Color = [Gratings_400px_ks_v3.Color v_line_transp];
        Gratings_400px_ks_v4 = vline(cl_var.Gratings_400px_trig_times_vec(5),'r'); Gratings_400px_ks_v4.Color = [Gratings_400px_ks_v4.Color v_line_transp];
        Gratings_400px_ks_v5 = vline(cl_var.Gratings_400px_trig_times_vec(6),'r'); Gratings_400px_ks_v5.Color = [Gratings_400px_ks_v5.Color v_line_transp];
        Gratings_400px_ks_v6 = vline(cl_var.Gratings_400px_trig_times_vec(7),'r'); Gratings_400px_ks_v6.Color = [Gratings_400px_ks_v6.Color v_line_transp];
        Gratings_400px_ks_v7 = vline(cl_var.Gratings_400px_trig_times_vec(8),'r'); Gratings_400px_ks_v7.Color = [Gratings_400px_ks_v7.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{6}); % title('Gratings 400px');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Gratings_400px_ks_v1,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v2,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v3,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v4,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v5,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v6,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v7,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
        
        % Gratings 400px heatmap
        if Plot_Axis_Choice == 1
            ax10 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(6)+1);
        else % Plot_Axis_Choice == 2
            ax10 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(6)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        imagesc(Gratings_400px_ksdensity_grid,[],Gratings_400px_Mat_Clusi,[Gratings_400px_min_resp Gratings_400px_max_resp]); colormap(ax10,gray(256)); hold on;
        Gratings_400px_ks_v1 = vline(cl_var.Gratings_400px_trig_times_vec(2),'r'); Gratings_400px_ks_v1.Color = [Gratings_400px_ks_v1.Color v_line_transp];
        Gratings_400px_ks_v2 = vline(cl_var.Gratings_400px_trig_times_vec(3),'r'); Gratings_400px_ks_v2.Color = [Gratings_400px_ks_v2.Color v_line_transp];
        Gratings_400px_ks_v3 = vline(cl_var.Gratings_400px_trig_times_vec(4),'r'); Gratings_400px_ks_v3.Color = [Gratings_400px_ks_v3.Color v_line_transp];
        Gratings_400px_ks_v4 = vline(cl_var.Gratings_400px_trig_times_vec(5),'r'); Gratings_400px_ks_v4.Color = [Gratings_400px_ks_v4.Color v_line_transp];
        Gratings_400px_ks_v5 = vline(cl_var.Gratings_400px_trig_times_vec(6),'r'); Gratings_400px_ks_v5.Color = [Gratings_400px_ks_v5.Color v_line_transp];
        Gratings_400px_ks_v6 = vline(cl_var.Gratings_400px_trig_times_vec(7),'r'); Gratings_400px_ks_v6.Color = [Gratings_400px_ks_v6.Color v_line_transp];
        Gratings_400px_ks_v7 = vline(cl_var.Gratings_400px_trig_times_vec(8),'r'); Gratings_400px_ks_v7.Color = [Gratings_400px_ks_v7.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{6}); % title('Gratings 400px');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Gratings_400px_ks_v1,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v2,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v3,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v4,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v5,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v6,'LineWidth',v_line_width);
        set(Gratings_400px_ks_v7,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
    end
    
    %%% CNoise: Full RF Size
    if p.Obj_plot_vec(7) == 1
        Full_RF_Size_vec_Clusi = cl_var.Full_RF_Size_vec(clust_vec==FinalClusIndex_Vec(i));
        
        % CNoise: Full RF Size histogram
        if Plot_Axis_Choice == 1
            ax11 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(7));
        else % Plot_Axis_Choice == 2
            ax11 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(7),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        RF_Size_hist_loop = histogram(Full_RF_Size_vec_Clusi,RF_Size_hist_BinEdges_Final); % cl_var.RF_Size_hist_BinEdges
        %RF_Size_hist_loop = histogram(Full_RF_Size_vec_Clusi);
        xlim([RF_Size_hist_BinEdges_Final(1) RF_Size_hist_BinEdges_Final(end)]); % cl_var.RF_Size_hist_BinEdges(1) cl_var.RF_Size_hist_BinEdges(end
        if p.hist_y_lim == 1     % same limits across all panels
            ylim([0 RF_Size_hist_MaxValue]);
        else % p.hist_y_lim == 2 % different limits on each panel
            ylim([0 max(RF_Size_hist_loop.Values)]);
        end
        if i==First_Entry
            title(p.Plot_obj_name_vec{7}); % title('RF size');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('num. pixels');
        end
        set(gca,'FontSize',12);
    end
    
    %%% CNoise: Full RF Ellipticity
    if p.Obj_plot_vec(8) == 1
        Full_RF_Ellipticity_vec_Clusi = cl_var.Full_RF_Ellipticity_vec(clust_vec==FinalClusIndex_Vec(i));
        
        % CNoise: Full RF Size histogram
        if Plot_Axis_Choice == 1
            ax12 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(8));
        else % Plot_Axis_Choice == 2
            ax12 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(8),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        RF_Ellip_hist_loop = histogram(Full_RF_Ellipticity_vec_Clusi,RF_Ellip_hist_BinEdges_Final); % cl_var.RF_Ellip_hist_BinEdges
        xlim([RF_Ellip_hist_BinEdges_Final(1) RF_Ellip_hist_BinEdges_Final(end)]); % cl_var.RF_Ellip_hist_BinEdges(1) cl_var.RF_Ellip_hist_BinEdges(end)
        if p.hist_y_lim == 1     % same limits across all panels
            ylim([0 RF_Ellip_hist_MaxValue]);
        else % p.hist_y_lim == 2 % different limits on each panel
            if max(RF_Ellip_hist_loop.Values)~=0
                ylim([0 max(RF_Ellip_hist_loop.Values)]);
            else
                ylim([0 1]);
            end
        end
        if i==First_Entry
            title(p.Plot_obj_name_vec{8}); % title('RF ellipticity');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('RF ellipticity');
        end
        set(gca,'FontSize',12);
    end
    
    %%% CNoise: Full RF Dominant Axis Angle
    if p.Obj_plot_vec(9) == 1
        Full_RF_Dom_Ax_Ang_vec_Clusi = cl_var.Full_RF_Dom_Ax_Ang_vec(clust_vec==FinalClusIndex_Vec(i));
        
        % CNoise: Full RF Dominant Axis Angle histogram
        if Plot_Axis_Choice == 1
            ax13 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(9));
        else % Plot_Axis_Choice == 2
            ax13 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(9),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing_Polar,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        %RF_Dom_Ax_Ang_hist_loop = histogram(Full_RF_Dom_Ax_Ang_vec_Clusi,RF_Dom_Ax_Ang_hist_BinEdges_Final);
        RF_Dom_Ax_Ang_hist_loop = polarhistogram(Full_RF_Dom_Ax_Ang_vec_Clusi,(pi/180)*[-90 -67.5 -45 -22.5 0 22.5 45 67.5 90]);
        %xlim([RF_Dom_Ax_Ang_hist_BinEdges_Final(1) RF_Dom_Ax_Ang_hist_BinEdges_Final(end)]);
        %if p.hist_y_lim == 1     % same limits across all panels
        %    ylim([0 RF_Dom_Ax_Ang_hist_MaxValue]);
        %else % p.hist_y_lim == 2 % different limits on each panel
        %    ylim([0 max(RF_Dom_Ax_Ang_hist_loop.Values)]);
        %end
        %set(gca,'ThetaTick',[0 45 90 270 315]);
        %set(gca,'ThetaTickLabel',{'0';'45';'90';'-90';'-45'});
        set(gca,'ThetaTickLabel',[],'RTickLabel',[]);
        if i==First_Entry
            title(p.Plot_obj_name_vec{9}); % title('RF dom. axis angle');
            %title('x'); % Test
        end
        %if i==Last_Entry
        %    xlabel('dom. axis angle');
        %end
        set(gca,'FontSize',12);
    end
    
    
    %%% FFF2
    if p.Obj_plot_vec(10) == 1
        FFF2_Mat_Clusi = FFF2_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        FFF2_mean_mat_Clusi = FFF2_Clus_mean_mat(FinalClusIndex_Vec(i),:); % PAR Mod 15,11,2021
        FFF2_std_mat_Clusi  = FFF2_Clus_std_mat(FinalClusIndex_Vec(i),:);  % PAR Mod 15,11,2021
        
        
        % FFF2 trace
        if Plot_Axis_Choice == 1
            ax14 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(10));
        else % Plot_Axis_Choice == 2
            ax14 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(10),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(FFF2_ksdensity_grid,FFF2_mean_mat_Clusi,FFF2_std_mat_Clusi,{'k-','markerfacecolor','k'}); hold on; % PAR Mod 15,11,2021
        xlim([0 cl_var.FFF2_stim_end_time]);
        %ylim([FFF2_trace_LwrBd FFF2_trace_UprBd]);
        FFF2_ks_v1 = vline(cl_var.FFF2_trig_times_vec(2),'r'); FFF2_ks_v1.Color = [FFF2_ks_v1.Color v_line_transp];
        FFF2_ks_v2 = vline(cl_var.FFF2_trig_times_vec(3),'r'); FFF2_ks_v2.Color = [FFF2_ks_v2.Color v_line_transp];
        FFF2_ks_v3 = vline(cl_var.FFF2_trig_times_vec(4),'r'); FFF2_ks_v3.Color = [FFF2_ks_v3.Color v_line_transp];
        FFF2_ks_v4 = vline(cl_var.FFF2_trig_times_vec(5),'r'); FFF2_ks_v4.Color = [FFF2_ks_v4.Color v_line_transp];
        FFF2_ks_v5 = vline(cl_var.FFF2_trig_times_vec(6),'r'); FFF2_ks_v5.Color = [FFF2_ks_v5.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{10}); % title('FFF2');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(FFF2_ks_v1,'LineWidth',v_line_width);
        set(FFF2_ks_v2,'LineWidth',v_line_width);
        set(FFF2_ks_v3,'LineWidth',v_line_width);
        set(FFF2_ks_v4,'LineWidth',v_line_width);
        set(FFF2_ks_v5,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
        
        % FFF2 heatmap
        if Plot_Axis_Choice == 1
            ax15 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(10)+1);
        else % Plot_Axis_Choice == 2
            ax15 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(10)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        %imagesc(FFF2_ksdensity_grid,[],FFF2_Mat_Clusi,[FFF2_min_resp FFF2_max_resp]); colormap(ax15,gray(256)); hold on;
        imagesc(FFF2_ksdensity_grid,[],FFF2_Mat_Clusi); colormap(ax15,parula(256)); hold on; % colormap(ax15,gray(256)); hold on; % PAR Mod 15,11,2021
        FFF2_ks_v1 = vline(cl_var.FFF2_trig_times_vec(2),'r'); FFF2_ks_v1.Color = [FFF2_ks_v1.Color v_line_transp];
        FFF2_ks_v2 = vline(cl_var.FFF2_trig_times_vec(3),'r'); FFF2_ks_v2.Color = [FFF2_ks_v2.Color v_line_transp];
        FFF2_ks_v3 = vline(cl_var.FFF2_trig_times_vec(4),'r'); FFF2_ks_v3.Color = [FFF2_ks_v3.Color v_line_transp];
        FFF2_ks_v4 = vline(cl_var.FFF2_trig_times_vec(5),'r'); FFF2_ks_v4.Color = [FFF2_ks_v4.Color v_line_transp];
        FFF2_ks_v5 = vline(cl_var.FFF2_trig_times_vec(6),'r'); FFF2_ks_v5.Color = [FFF2_ks_v5.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{10}); % title('FFF2');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(FFF2_ks_v1,'LineWidth',v_line_width);
        set(FFF2_ks_v2,'LineWidth',v_line_width);
        set(FFF2_ks_v3,'LineWidth',v_line_width);
        set(FFF2_ks_v4,'LineWidth',v_line_width);
        set(FFF2_ks_v5,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
    end
    
    
    %%% Chirp2
    if p.Obj_plot_vec(11) == 1
        Chirp2_Mat_Clusi = Chirp2_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        
        % Chirp2 trace
        if Plot_Axis_Choice == 1
            ax16 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(11));
        else % Plot_Axis_Choice == 2
            ax16 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(11),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(Chirp2_ksdensity_grid,Chirp2_Clus_mean_mat(i,:),Chirp2_Clus_std_mat(i,:),{'k-','markerfacecolor','k'}); hold on;
        xlim([0 cl_var.Chirp2_stim_end_time]);
        %ylim([Chirp2_trace_LwrBd Chirp2_trace_UprBd]);
        Chirp2_ks_v1 = vline(cl_var.Chirp2_trig_times_vec(2),'r'); Chirp2_ks_v1.Color = [Chirp2_ks_v1.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{11}); % title('Chirp2');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Chirp2_ks_v1,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
        
        % Chirp2 heatmap
        if Plot_Axis_Choice == 1
            ax17 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(11)+1);
        else % Plot_Axis_Choice == 2
            ax17 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(11)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        imagesc(Chirp2_ksdensity_grid,[],Chirp2_Mat_Clusi,[Chirp2_min_resp Chirp2_max_resp]); colormap(ax17,gray(256)); hold on;
        Chirp2_ks_v1 = vline(cl_var.Chirp2_trig_times_vec(2),'r'); Chirp2_ks_v1.Color = [Chirp2_ks_v1.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{11}); % title('Chirp2');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Chirp2_ks_v1,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
    end
    
    %%% SSub
    if p.Obj_plot_vec(12) == 1
        SSub_Mat_Clusi = SSub_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        
        % SSub trace
        if Plot_Axis_Choice == 1
            ax18 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(12));
        else % Plot_Axis_Choice == 2
            ax18 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(12),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(SSub_ksdensity_grid,SSub_Clus_mean_mat(i,:),SSub_Clus_std_mat(i,:),{'k-','markerfacecolor','k'}); hold on;
        xlim([0 cl_var.SSub_stim_end_time]);
        %ylim([SSub_trace_LwrBd SSub_trace_UprBd]);
        for j = 1:p.SSub_NumTrigPerRep-1
            SSub_ks_loop = vline(cl_var.SSub_trig_times_vec(j+1),'r'); SSub_ks_loop.Color = [SSub_ks_loop.Color v_line_transp];
            set(SSub_ks_loop,'LineWidth',v_line_width);
        end
        if i==First_Entry
            title(p.Plot_obj_name_vec{12}); % title('SSub');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(gca,'FontSize',12);
        
        % SSub heatmap
        if Plot_Axis_Choice == 1
            ax19 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(12)+1);
        else % Plot_Axis_Choice == 2
            ax19 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(12)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        imagesc(SSub_ksdensity_grid,[],SSub_Mat_Clusi,[SSub_min_resp SSub_max_resp]); colormap(ax19,gray(256)); hold on;
        for j = 1:p.SSub_NumTrigPerRep-1
            SSub_ks_loop = vline(cl_var.SSub_trig_times_vec(j+1),'r'); SSub_ks_loop.Color = [SSub_ks_loop.Color v_line_transp];
            set(SSub_ks_loop,'LineWidth',v_line_width);
        end
        if i==First_Entry
            title(p.Plot_obj_name_vec{12}); % title('SSub');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(gca,'FontSize',12);
    end
    
    
    %%% CSteps
    if p.Obj_plot_vec(13) == 1
        CSteps_Mat_Clusi      = CSteps_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        CSteps_mean_mat_Clusi = CSteps_Clus_mean_mat(FinalClusIndex_Vec(i),:); % PAR Mod 24,11,2021
        CSteps_std_mat_Clusi  = CSteps_Clus_std_mat(FinalClusIndex_Vec(i),:);  % PAR Mod 24,11,2021
        
        % CSteps trace
        if Plot_Axis_Choice == 1
            ax20 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(13));
        else % Plot_Axis_Choice == 2
            ax20 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(13),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(CSteps_ksdensity_grid,CSteps_mean_mat_Clusi,CSteps_std_mat_Clusi,{'k-','markerfacecolor','k'}); hold on; % PAR Mod 24,11,2021 % shadedErrorBar(CSteps_ksdensity_grid,CSteps_Clus_mean_mat(i,:),CSteps_Clus_std_mat(i,:),{'k-','markerfacecolor','k'}); hold on;
        xlim([0 cl_var.CSteps_stim_end_time]);
        %ylim([CSteps_trace_LwrBd CSteps_trace_UprBd]);
        for j = 1:p.CSteps_NumTrigPerRep-1
            CSteps_ks_loop = vline(cl_var.CSteps_trig_times_vec(j+1),'r'); CSteps_ks_loop.Color = [CSteps_ks_loop.Color v_line_transp];
            set(CSteps_ks_loop,'LineWidth',v_line_width);
        end
        if i==First_Entry
            title(p.Plot_obj_name_vec{13}); % title('CSteps');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(gca,'FontSize',12);
        
        % CSteps heatmap
        if Plot_Axis_Choice == 1
            ax21 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(13)+1);
        else % Plot_Axis_Choice == 2
            ax21 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(13)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        %imagesc(CSteps_ksdensity_grid,[],CSteps_Mat_Clusi,[CSteps_min_resp CSteps_max_resp]); colormap(ax21,gray(256)); hold on;
        imagesc(CSteps_ksdensity_grid,[],CSteps_Mat_Clusi); colormap(ax21,parula(256)); hold on; % PAR Mod 24,11,2021
        for j = 1:p.CSteps_NumTrigPerRep-1
            CSteps_ks_loop = vline(cl_var.CSteps_trig_times_vec(j+1),'r'); CSteps_ks_loop.Color = [CSteps_ks_loop.Color v_line_transp];
            set(CSteps_ks_loop,'LineWidth',v_line_width);
        end
        if i==First_Entry
            title(p.Plot_obj_name_vec{13}); % title('CSteps');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(gca,'FontSize',12);
    end
    
    
    %%% Chirp3
    if p.Obj_plot_vec(14) == 1
        Chirp3_Mat_Clusi      = Chirp3_spike_density_mat(clust_vec==FinalClusIndex_Vec(i),:);
        Chirp3_mean_mat_Clusi = Chirp3_Clus_mean_mat(FinalClusIndex_Vec(i),:); % PAR Mod 24,11,2021
        Chirp3_std_mat_Clusi  = Chirp3_Clus_std_mat(FinalClusIndex_Vec(i),:);  % PAR Mod 24,11,2021
        
        % Chirp3 trace
        if Plot_Axis_Choice == 1
            ax22 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(14));
        else % Plot_Axis_Choice == 2
            ax22 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(14),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        shadedErrorBar(Chirp3_ksdensity_grid,Chirp3_mean_mat_Clusi,Chirp3_std_mat_Clusi,{'k-','markerfacecolor','k'}); hold on; % PAR Mod 24,11,2021 % shadedErrorBar(Chirp3_ksdensity_grid,Chirp3_Clus_mean_mat(i,:),Chirp3_Clus_std_mat(i,:),{'k-','markerfacecolor','k'});
        xlim([0 cl_var.Chirp3_stim_end_time]);
        %ylim([Chirp3_trace_LwrBd Chirp3_trace_UprBd]);
        Chirp3_ks_v1 = vline(cl_var.Chirp3_trig_times_vec(2),'r'); Chirp3_ks_v1.Color = [Chirp3_ks_v1.Color v_line_transp];
        Chirp3_ks_v2 = vline(cl_var.Chirp3_trig_times_vec(3),'r'); Chirp3_ks_v2.Color = [Chirp3_ks_v2.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{14}); % title('Chirp3');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Chirp3_ks_v1,'LineWidth',v_line_width);
        set(Chirp3_ks_v2,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
        
        % Chirp3 heatmap
        if Plot_Axis_Choice == 1
            ax23 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(14)+1);
        else % Plot_Axis_Choice == 2
            ax23 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(14)+1,...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        %imagesc(Chirp3_ksdensity_grid,[],Chirp3_Mat_Clusi,[Chirp3_min_resp Chirp3_max_resp]); colormap(ax23,gray(256)); hold on;
        imagesc(Chirp3_ksdensity_grid,[],Chirp3_Mat_Clusi,[Chirp3_min_resp Chirp3_max_resp]); colormap(ax23,parula(256)); hold on; % PAR Mod 24,11,2021
        Chirp3_ks_v1 = vline(cl_var.Chirp3_trig_times_vec(2),'r'); Chirp3_ks_v1.Color = [Chirp3_ks_v1.Color v_line_transp];
        Chirp3_ks_v2 = vline(cl_var.Chirp3_trig_times_vec(3),'r'); Chirp3_ks_v2.Color = [Chirp3_ks_v2.Color v_line_transp];
        if i==First_Entry
            title(p.Plot_obj_name_vec{14}); % title('Chirp3');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('time (s)');
        end
        set(Chirp3_ks_v1,'LineWidth',v_line_width);
        set(Chirp3_ks_v2,'LineWidth',v_line_width);
        set(gca,'FontSize',12);
    end
    
    
    %%% Data Set
    if p.Obj_plot_vec(15) == 1
        
        % Data Set bar chart
        if Plot_Axis_Choice == 1
            ax24 = subplot(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(15));
        else % Plot_Axis_Choice == 2
            ax24 = subaxis(Fig_Rows,Num_Columns,Num_Columns*(i-First_Entry)+Col_index_vec(15),...
                'SpacingVertical',Vert_Spacing,'SpacingHorizontal',Horiz_Spacing,...
                'MarginLeft',Left_Margin,'MarginRight',Right_Margin,...
                'MarginTop',Top_Margin,'MarginBottom',Bottom_Margin);
        end
        if Norm_Data_Bar_Choice == 1 % Don't normalise data set bar plots
            Cell_per_Cluster_bar_loop = bar(Gen_Data_Names_bar_vec,Cell_per_Cluster_mat(i,:));
        else % Norm_Data_Bar_Choice == 2 % Normalise data set bar plots
            Cell_per_Cluster_bar_loop = bar(Gen_Data_Names_bar_vec,Cell_per_Cluster_mat_norm(i,:));
        end
        %xlim([cl_var.RF_Ellip_hist_BinEdges(1) cl_var.RF_Ellip_hist_BinEdges(end)]);
        if Norm_Data_Bar_Choice == 1 % Don't normalise data set bar plots
            if p.hist_y_lim == 1     % same limits across all panels
                ylim([0 Cell_per_Cluster_bar_MaxValue]);
            else % p.hist_y_lim == 2 % different limits on each panel
                ylim([0 max(Cell_per_Cluster_mat(i,:))]);
            end
        else % Norm_Data_Bar_Choice == 2 % Normalise data set bar plots
            if p.hist_y_lim == 1     % same limits across all panels
                ylim([0 1]);
            else % p.hist_y_lim == 2 % different limits on each panel
                ylim([0 max(Cell_per_Cluster_mat_norm(i,:))]);
            end
        end
        if i==First_Entry
            title(p.Plot_obj_name_vec{15}); % title('data set');
        end
        if xtick_Choice == 2 && i < Last_Entry
            set(gca,'xticklabel',[]);
        end
        if i==Last_Entry
            xlabel('data set');
        end
        set(gca,'FontSize',12);
    end
    
end

set(gcf,'color','w');

%% Plot Cluster Dendrogram (New - 18,01,2022)

%%% If standardise first
% All_Scores_Dendro = zscore(All_Scores,1); % for a matrix X, if dim = 1, then zscore uses the means and standard deviations along the columns of X
All_Scores_Dendro = NaN(size(All_Scores,1),size(All_Scores,2));
for i = 1:size(All_Scores,2)
    All_Scores_Dendro(:,i) = (All_Scores(:,i) - nanmean(All_Scores(:,i)))/nanstd(All_Scores(:,i));
end
% Test
% figure;
% histogram(All_Scores_Dendro(:,50));
%%% If don't standardise
%All_Scores_Dendro = All_Scores;

MeanClustScores = NaN(Num_clust_True,size(All_Scores_Dendro,2)); % Rows cluster averages, columns score variables

for i = 1:Num_clust_True
    MeanClustScores(i,:) = nanmean(All_Scores_Dendro(clust_vec==FinalClusIndex_Vec(i),:),1);
end


Z = linkage(MeanClustScores,'average','correlation'); % linkage(MeanClustScores,'average','correlation')
D = pdist(MeanClustScores);
leafOrder = optimalleaforder(Z,D);
figure;
dendro = dendrogram(Z,100,'reorder',leafOrder,'Orientation','left','ColorThreshold','default');
set(dendro,'LineWidth',1.5);
set(gca,'FontSize',12);
set(gcf,'color','w');
% figure;
% dendrogram(Z);



%% Save Data

%save data_DP3_RFC4_14_06_2021_1;
%save data_DP3_RFC4_Chick_NoUV_27_05_2021_1;
%save data_DP3_RFC4_Chick_UV_27_05_2021_1;

%save data_DP4_RFC5_14_06_2021_1;
%save data_DP4_RFC5_22_06_2021_1;
%save data_DP4_RFC5_22_06_2021_3;
%save data_DP4_RFC5_22_06_2021_5;
%save data_DP4_RFC5_22_06_2021_7;
%save data_DP4_RFC5_22_06_2021_8;
%save data_DP4_RFC5_22_06_2021_9;

%save data_DP4_RFC5_07_07_2021_2;

%save data_DP4_RFC5_07_09_2021_1;
%save data_DP4_RFC5_07_09_2021_2;
%save data_DP4_RFC5_07_09_2021_3;
%save data_DP4_RFC5_07_09_2021_4;
%save data_DP4_RFC5_07_09_2021_5;
%save data_DP4_RFC5_07_09_2021_6;
%save data_DP4_RFC5_07_09_2021_7;
%save data_DP4_RFC5_07_09_2021_8;

%save data_DP5_RFC6_24_09_2021_1;
%save data_DP5_RFC6_24_09_2021_2;
%save data_DP5_RFC6_24_09_2021_2_2;
%save data_DP5_RFC6_24_09_2021_3_AndFNoise;
%save data_DP5_RFC6_24_09_2021_3_2_AndFNoise;

%save data_DP5_RFC6_04_11_2021_1;
%save data_DP5_RFC6_04_11_2021_2;

%save data_DP5_RFC6_15_11_2021_1;
%save data_DP5_RFC6_15_11_2021_2;

%save('data_DP5_RFC6_19_11_2021_1.mat','-v7.3');
%save('data_DP5_RFC6_19_11_2021_2_2.mat.mat','-v7.3');
%save('data_DP5_RFC6_29_11_2021_2.mat','-v7.3');
%save('data_DP5_RFC6_29_11_2021_3.mat','-v7.3');
%save('data_DP5_RFC6_29_11_2021_4.mat','-v7.3');
%save('data_DP5_RFC6_29_11_2021_5.mat','-v7.3');

%save('data_DP5_RFC6_17_01_2022_1.mat','-v7.3');

