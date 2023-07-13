function [hier_clust_vec,hier_clust_num_clust,cophenetic_dist] = Hier_Clust_fn_2(All_Scores,p)
% 06,05,2021 Mod: Replacing dom axis angle NaNs with -1.

% Modify dom axis angle where it exists for clustering
if p.Obj_clust_vec(7) == 1
    All_Scores(isnan(All_Scores(:,end)),end) = -1;
end


if sum(isnan(All_Scores),'all') == 0
    pdist_vec   = pdist(All_Scores,p.pdist_opt); % pdist(All_Scores,p.pdist_opt), pdist(All_Scores,'cityblock')
else % If there are NaN entries in AllScores due to missing data
    pdist_vec   = pdist(All_Scores,@naneucdist); % pdist(All_Scores,p.pdist_opt), pdist(All_Scores,'cityblock')
end
%Input data, specified as a numeric matrix of size m-by-n. Rows correspond 
%to individual observations, and columns correspond to individual variables.
% rows = cells, columns = variables
% (some) dist options:
% euclidean (default)
% mahalanobis
% cityblock
% minkowski (default exponent = 2 (same as euclidean when = 2, same as city block when = 1, same as Chebychev when = inf))
% chebychev
% correlation
% spearman

linkage_mat = linkage(pdist_vec,p.linkage_opt); % linkage(pdist_vec,'average')
% linkage options:
% single (default)
% average
% centroid (with euclidean dist only)
% complete
% median   (with euclidean dist only)
% ward     (with euclidean dist only)
% weighted

% Plot dendrogram
% figure;
% subplot(1,2,1);
% dendrogram(linkage_mat); % Plots no more than 30 leaf nodes unless specify otherwise.
% subplot(1,2,2);
% dendrogram(linkage_mat,0);

% Verify Dissimilarity
cophenetic_dist = cophenet(linkage_mat,pdist_vec);

% Verify Consistency
inconsistency = inconsistent(linkage_mat);
% a link whose height differs noticeably from the height of the links below it indicates that the objects
% joined at this level in the cluster tree are much farther apart from each other than their components 
% were when they were joined. This link is said to be inconsistent with the links below it.
% In cluster analysis, inconsistent links can indicate the border of a natural division in a data set
% Links that join distinct clusters have a high inconsistency coefficient; links that join indistinct 
% clusters have a low inconsistency coefficient.

%cutoff = 1.165; % (1.165 <-> 10 clusters)

%distance_cutoff = 3;
%distance_cutoff = median([linkage_mat(end-9,3) linkage_mat(end-8,3)]);
%distance_cutoff = median([linkage_mat(end-(p.Num_Clus-1),3) linkage_mat(end-(p.Num_Clus-2),3)]);

% Cluster data
if     p.Info_crit == 1 % inconsistency cutoff
    hier_clust_vec = cluster(linkage_mat,'cutoff',p.inconsistency_cutoff);
elseif p.Info_crit == 2 % distance cutoff
    hier_clust_vec = cluster(linkage_mat,'cutoff',p.distance_cutoff,'criterion','distance');
else % p.Info_crit == 3 % prescribe number of clusters
    hier_clust_vec = cluster(linkage_mat,'MaxClust',p.Num_Clus);
end
%hier_clust_vec = cluster(linkage_mat,'cutoff',cutoff);
%hier_clust_vec  = cluster(linkage_mat,'cutoff',distance_cutoff,'criterion','distance');
%hier_clust_vec = cluster(linkage_mat,'MaxClust',p.Num_Clus);
%
%If the criterion for defining clusters is 'distance', then cluster groups all leaves at or below a node
%into a cluster, provided that the height of the node is less than C.
%
%If the criterion for defining clusters is 'inconsistent', then the inconsistent values of a node and all
%its subnodes must be less than C for cluster to group them into a cluster. cluster begins from the root
% of the cluster tree Z and steps down through the tree until it encounters a node whose inconsistent value
%is less than the threshold C, and whose subnodes (or descendants) have inconsistent values less than C. 
%Then cluster groups all leaves at or below the node into a cluster (or a singleton if the node itself is 
%a leaf). cluster follows every branch in the tree until all leaf nodes are in clusters.

hier_clust_unique_vec = unique(hier_clust_vec);        % unique cluster labels
hier_clust_num_clust  = length(hier_clust_unique_vec); % number of clusters

% Plot consistency dist
if p.Info_crit == 1 % inconsistency cutoff
    figure;
    histogram(inconsistency(:,4),200); hold on;
    incon_v1 = vline(p.inconsistency_cutoff, 'r','inconsistency threshold');
    xlabel('inconsistency');
    ylabel('density');
    set(incon_v1,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
end

% Plot dendrogram
% figure;
% dendrogram(linkage_mat,'ColorThreshold',distance_cutoff); % linkage_mat,0,'ColorThreshold',distance_cutoff
% There is no good way to colour the dendrogram if I use cutoff.

figure;
if p.Info_crit == 2 % distance cutoff
    subplot(1,2,1);
    dendrogram(linkage_mat,'ColorThreshold',p.distance_cutoff); hold on;% Plots no more than 30 leaf nodes unless specify otherwise.
    dist_h1 = hline(p.distance_cutoff,'b--','distance cutoff');
    title('first 30 leaf nodes');
    xlabel('branch');
    ylabel('distance');
    set(dist_h1,'LineWidth',1.5);
    set(gca,'FontSize',12);
    subplot(1,2,2);
    dendrogram(linkage_mat,0,'ColorThreshold',p.distance_cutoff); hold on;
    dist_h2 = hline(p.distance_cutoff,'b--','distance cutoff');
    title('full dendrogram');
    xlabel('branch');
    ylabel('distance');
    set(dist_h2,'LineWidth',1.5);
    set(gca,'FontSize',12);
    set(gcf,'color','w');
else % inconsistency cutoff or prescribe number of clusters
    subplot(1,2,1);
    dendrogram(linkage_mat);% Plots no more than 30 leaf nodes unless specify otherwise.
    title('first 30 leaf nodes');
    xlabel('branch');
    ylabel('distance');
    set(gca,'FontSize',12);
    subplot(1,2,2);
    dendrogram(linkage_mat,0);
    title('full dendrogram');
    xlabel('branch');
    ylabel('distance');
    set(gca,'FontSize',12);
    set(gcf,'color','w');
end


