function [gm,aic,bic,converged,allConverge,clusterOpt,k_Opt] = GMM_Clust_fn_3(All_Scores,p)
% 29,11,2021 Mod: Correcting ln 128 for choice of gmopt for p.Info_crit == 3.

% Modify dom axis angle where it exists to remove NaN for clustering
if p.Obj_clust_vec(7) == 1
%     DomAxAng_Col_index = sum(p.Obj_clust_vec(1:7));
%     All_Scores(isnan(All_Scores(:,DomAxAng_Col_index)),DomAxAng_Col_index) = -1;
    All_Scores(isnan(All_Scores(:,p.PC_Num_Partial)),p.PC_Num_Partial) = -1; % PAR Mod 08,07,2021
end
% Modify ellipticity where it exists to remove NaN for clustering
if p.Obj_clust_vec(6) == 1 % PAR Mod 08,07,2021 (All new)
    All_Scores(isnan(All_Scores(:,p.PC_Num_Partial_2)),p.PC_Num_Partial_2) = -1;
end


%%% Settings
% Top 3 now defined in CL
% nK       = numel(p.k_vec);                % number of clusters examined
% nSigma   = numel(p.Sigma);                % 2 covariance matrix types examined
% nSC      = numel(p.SharedCovariance);     % 2 covarience options
options  = statset('MaxIter',p.Max_Iter); % increase beyond the std. num. of iterations (100). Matlab help: 10000.                
%Num_Comp = p.Total_PC_Num;               % number of principal components to use. (unused)

%%% Preallocation
gm        = cell(p.nK,p.nSigma,p.nSC);     % GMM output
aic       = zeros(p.nK,p.nSigma,p.nSC);    % AIC output
bic       = zeros(p.nK,p.nSigma,p.nSC);    % BIC output
converged = false(p.nK,p.nSigma,p.nSC);    % convergence checks record

if p.Parpool == 1 % parallel processing - on
    
    parpool(p.Num_Cores);
    for m = 1:p.nSC
        for j = 1:p.nSigma
            parfor i = 1:p.nK
                
                rng(1); % for reproducability (always get the same random numbers)
                
                gm{i,j,m} = fitgmdist(All_Scores,p.k_vec(i),...
                    'CovarianceType',p.Sigma{j},...
                    'SharedCovariance',p.SharedCovariance{m},...
                    'RegularizationValue',p.RegularizationValue,...
                    'Replicates',p.Replicate_Num,...
                    'Options',options);
                
                aic(i,j,m) = gm{i,j,m}.AIC;
                bic(i,j,m) = gm{i,j,m}.BIC;
                
                converged(i,j,m) = gm{i,j,m}.Converged;
                
                disp([m j i]);
                
            end
        end
    end
    
else % p.Parpool == 2 % parallel processing - off
    
    for m = 1:p.nSC
        for j = 1:p.nSigma
            for i = 1:p.nK
                
                rng(1); % for reproducability (always get the same random numbers)
                
                gm{i,j,m} = fitgmdist(All_Scores,p.k_vec(i),...
                    'CovarianceType',p.Sigma{j},...
                    'SharedCovariance',p.SharedCovariance{m},...
                    'RegularizationValue',p.RegularizationValue,...
                    'Replicates',p.Replicate_Num,...
                    'Options',options);
                
                aic(i,j,m) = gm{i,j,m}.AIC;
                bic(i,j,m) = gm{i,j,m}.BIC;
                
                converged(i,j,m) = gm{i,j,m}.Converged;
                
                disp([m j i]);
                
            end
        end
    end
    
end

allConverge = (sum(converged(:)) == p.nK*p.nSigma*p.nSC); % check that all solutions converged

if allConverge==1
   
    disp('All GMMs successfully converged');
    
else
    
    disp('Not all GMMs successfully converged');
    
end

% Stop parallel pool
if p.Parpool == 1 % parallel processing - on
    delete(gcp('nocreate'));
end


%% Find the optimum clustering

aic_opt_value = min(aic,[],'all');
bic_opt_value = min(bic,[],'all');

gmOpt_BIC      = gm{bic==bic_opt_value};
clusterOpt_BIC = cluster(gmOpt_BIC,All_Scores);
k_Opt_BIC      = gmOpt_BIC.NumComponents;
gmOpt_AIC      = gm{aic==aic_opt_value};
clusterOpt_AIC = cluster(gmOpt_AIC,All_Scores);
k_Opt_AIC      = gmOpt_AIC.NumComponents;

if p.Info_crit == 1
    % Use the best clustering option as judged by the BIC
    gmOpt      = gmOpt_BIC;
    clusterOpt = clusterOpt_BIC;
    k_Opt      = k_Opt_BIC;
elseif p.Info_crit == 2
    % Use the best clustering option as judged by the AIC
    gmOpt      = gmOpt_AIC;
    clusterOpt = clusterOpt_AIC;
    k_Opt      = k_Opt_AIC;
else % p.Info_crit == 3
    % Choose own number of clusters
    k_Opt      = p.Num_Clus;
    gmOpt      = gm{k_Opt,p.Sig_opt,p.Sig_share_opt}; % PAR Mod 29,11,2021 (was gm{k_Opt))
    clusterOpt = cluster(gmOpt,All_Scores);
end


%% Plot AICs and BICs for each combination or cluster numbers and covariance options

% Reshape so that:
% rows    = ks and 
% columns = cov: diagonal-shared,full-shared,diagonal-unshared,full-unshared
aic_mat = reshape(aic,p.nK,p.nSigma*p.nSC); 
bic_mat = reshape(bic,p.nK,p.nSigma*p.nSC);

figure;
subplot(2,1,1);
plot(p.k_vec,aic_mat(:,1)); hold on;
plot(p.k_vec,aic_mat(:,2));
plot(p.k_vec,aic_mat(:,3));
plot(p.k_vec,aic_mat(:,4));
plot(k_Opt_AIC ,aic_opt_value,'ro'); % plots the optimum soln
title('AIC for various $k$ and $\Sigma$ choices','Interpreter','latex');
xlabel('number of clusters, $k$','Interpreter','Latex');
ylabel('AIC');
legend({'diagonal-shared','full-shared','diagonal-unshared',...
    'full-unshared'},'Location','northeast');
set(gca,'FontSize',12);
subplot(2,1,2);
plot(p.k_vec,bic_mat(:,1)); hold on;
plot(p.k_vec,bic_mat(:,2));
plot(p.k_vec,bic_mat(:,3));
plot(p.k_vec,bic_mat(:,4));
plot(k_Opt_BIC ,bic_opt_value,'ro'); % plots the optimum soln
title('BIC for various $k$ and $\Sigma$ choices','Interpreter','latex')
xlabel('number of clusters, $k$','Interpreter','Latex');
ylabel('BIC');
set(gca,'FontSize',12);
set(gcf,'color','w');


