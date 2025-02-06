clear, clc, close all
%#ok<*SAGROW>
warning off
rng(1959)

% Path to the folder with the video data in MAT format
VideoPath = 'DataRGB\';

% Path to MATLAB tracks
TrackPaths{1} = 'TrReMATLAB\';

% Path to BASIC tracks
TrackPaths{2} = 'TrReBASIC\';

% Path to FCG tracks
TrackPaths{3} = 'TrFCG\';

addpath(genpath('FINCH')) % to add the FINCH clustering algorithm


% addpath(genpath(['C:\Lucy\Documents\RESEARCH\CODE\Clustering',...
%     '\ConstrainedOnlineClustering\CleanCode\']));


Videos = {'EP000002',...
    'EP000005',...
    'EP000009',...
    'EP000010',...
    'EP000016',...
    'EP000028',...
    'EP000033',...
    'EP000036',...
    'EP000060',...
    'EP000078',...
    'Koi_5652_952_540','Pigeons_4927_960_540_600f',...
    'Pigeons_8234_1280_720',...
    'Pigeons_29033_960_540_300f',...
    'Pigs_49651_960_540_500f'};

LinkageMethods = {'single', 'complete', 'average','weighted', ...
    'centroid', 'median', 'ward'};

m = numel(LinkageMethods);


% (A) Tracks on their own -------------------------------------------------

for i = 1:numel(Videos)

    video = Videos{i};
    % Data filename
    ofnRGB = [VideoPath,video,'_RGB.mat'];
    load(ofnRGB)

    for j = 1:3
        filename_tr_ML = [TrackPaths{j},video,'.csv'];
        M = readmatrix(filename_tr_ML);
        TrackLabels = M(:,2);

        A(i,j) = adjusted_rand_index(TrackLabels,Labels);
    end
end
fprintf('TRACKS ONLY -- done.\n\n')
%--------------------------------------------------------------------------

% (B) Clustering of raw data (unconstrained) ------------------------------

for i = 1:numel(Videos)

    video = Videos{i};
    % Data filename
    ofnRGB = [VideoPath,video,'_RGB.mat'];
    load(ofnRGB)

    DataRGB = zscore(DataRGB);
    c = max(Labels); % true number of clusters

    for ilink = 1:m
        AL = run_linkage(DataRGB,c,ilink);
        R(i,ilink) = adjusted_rand_index(AL,Labels);
        Clusters_Raw(i,ilink) = numel(unique(AL));
    end

    AL = run_kmeans(DataRGB,c);
    R(i,m+1) = adjusted_rand_index(AL,Labels);
    Clusters_Raw(i,m+1) = numel(unique(AL));

    AL = run_gmm(DataRGB,c);
    R(i,m+2) = adjusted_rand_index(AL,Labels);
    Clusters_Raw(i,m+2) = numel(unique(AL));

    AL = run_finch(DataRGB,c);
    R(i,m+3) = adjusted_rand_index(AL,Labels);
    Clusters_Raw(i,m+3) = numel(unique(AL));

    AL = run_spectral(DataRGB,c);
    R(i,m+4) = adjusted_rand_index(AL,Labels);
    Clusters_Raw(i,m+4) = numel(unique(AL));

    AL = run_dbscan(DataRGB,c);
    R(i,m+5) = adjusted_rand_index(AL,Labels);
    Clusters_Raw(i,m+5) = numel(unique(AL));


    fprintf('Raw data Video %i %s done.\n',i,video)

end

B = R;
save Res_NonConstrainedAB

%--------------------------------------------------------------------------

% (C) RAW DATA - Constrained ----------------------------------------------

R = [];

for i = 1:numel(Videos)

    video = Videos{i};
    % Data filename
    ofnRGB = [VideoPath,video,'_RGB.mat'];
    load(ofnRGB)
    DataRGB = zscore(DataRGB);

    for j = 1:3
        filename_tr_ML = [TrackPaths{j},video,'.csv'];
        M = readmatrix(filename_tr_ML);
        TrackLabels = M(:,2);
        Frames = M(:,1);


        c = max(Labels); % true number of clusters

        PARAM.L = c; % number of clusters
        PARAM.EnsembleSize = 5;
        PARAM.vis = 0;
        PARAM.BaseClusterer = "average";
        CC.MaxClusters = c;

        [CL,ML] = find_all_constraints(Frames,TrackLabels);
        AL = ccen(DataRGB, ML, CL, PARAM, CC);
        R(i,j) = adjusted_rand_index(AL,Labels);
        Clusters_Raw_Constrained(i,j) = numel(unique(AL));
    end

    fprintf('Raw-constrained -> Video %i %s done.\n',i,video)

end

C = R;

save Res_ConstrainedC
%--------------------------------------------------------------------------

% (D) TRACK CENTROIDS (unconstrained) -------------------------------------


R = [];

for i = 1:numel(Videos)

    video = Videos{i};
    % Data filename
    ofnRGB = [VideoPath,video,'_RGB.mat'];
    for j = 1:3
        load(ofnRGB)
        DataRGB = zscore(DataRGB);
        c = max(Labels); % true number of clusters


        filename_tr_ML = [TrackPaths{j},video,'.csv'];
        M = readmatrix(filename_tr_ML);
        TrackLabels = M(:,2);
        Frames = M(:,1);

        % find track centroids
        u = unique(TrackLabels);
        numtracks = numel(u); % number of tracks

        if numtracks > c % do clustering
            TrackCentroids = [];
            for jc = 1:numel(u)
                TrackCentroids = [TrackCentroids;...
                    mean(DataRGB(TrackLabels == u(jc),:),1)]; %#ok<*AGROW>
            end
            DataRGB = TrackCentroids;

            for ilink = 1:m
                AL = run_linkage(DataRGB,c,ilink);
                AL = reassemble_labels(AL,TrackLabels);
                R(j,i,ilink) = adjusted_rand_index(AL,Labels);
                Clusters_Tr(j,i,ilink) = numel(unique(AL));
            end

            AL = run_kmeans(DataRGB,c);
            AL = reassemble_labels(AL,TrackLabels);
            R(j,i,m+1) = adjusted_rand_index(AL,Labels);
            Clusters_Tr(j,i,m+1) = numel(unique(AL));

            AL = run_gmm(DataRGB,c);
            AL = reassemble_labels(AL,TrackLabels);
            R(j,i,m+2) = adjusted_rand_index(AL,Labels);
            Clusters_Tr(j,i,m+2) = numel(unique(AL));

            AL = run_finch(DataRGB,c);
            AL = reassemble_labels(AL,TrackLabels);
            R(j,i,m+3) = adjusted_rand_index(AL,Labels);
            Clusters_Tr(j,i,m+3) = numel(unique(AL));

            AL = run_spectral(DataRGB,c);
            AL = reassemble_labels(AL,TrackLabels);
            R(j,i,m+4) = adjusted_rand_index(AL,Labels);
            Clusters_Tr(j,i,m+4) = numel(unique(AL));

            AL = run_dbscan(DataRGB,c);
            AL = reassemble_labels(AL,TrackLabels);
            R(j,i,m+5) = adjusted_rand_index(AL,Labels);
            Clusters_Tr(j,i,m+5) = numel(unique(AL));
        else
            % no clustering required
            R(j,i,1:m+5) = adjusted_rand_index(TrackLabels,Labels);
            Clusters_Tr(j,i,1:m+5) = numtracks;
        end


    end

    fprintf('Track centroids data,  Video %i %s done.\n',i,video)

end
D = R;

save Res_NonConstrainedD

%--------------------------------------------------------------------------

% (E) TRACK CENTROIDS (constrained) -------------------------------------


R = [];

for i = 1:numel(Videos)

    video = Videos{i};
    % Data filename
    ofnRGB = [VideoPath,video,'_RGB.mat'];
    for j = 1:3
        load(ofnRGB)
        DataRGB = zscore(DataRGB);
        c = max(Labels); % true number of clusters


        filename_tr_ML = [TrackPaths{j},video,'.csv'];
        M = readmatrix(filename_tr_ML);
        TrackLabels = M(:,2);
        Frames = M(:,1);

        % find track centroids
        u = unique(TrackLabels);
        numtracks = numel(u); % number of tracks

        if numtracks > c % do clustering
            TrackCentroids = [];
            for jc = 1:numel(u)
                TrackCentroids = [TrackCentroids;...
                    mean(DataRGB(TrackLabels == u(jc),:),1)]; %#ok<*AGROW>
            end
            DataRGB = TrackCentroids;

            PARAM.L = c; % number of clusters
            PARAM.EnsembleSize = 5;
            PARAM.vis = 0;
            PARAM.BaseClusterer = "average";
            CC.MaxClusters = c;

            CL = find_cl_tracks(Frames,TrackLabels);
            AL = ccen(DataRGB, [], CL, PARAM, CC);
            AL = reassemble_labels(AL,TrackLabels);
            R(i,j) = adjusted_rand_index(AL,Labels);
            Clusters_Track_Constrained(i,j) = numel(unique(AL));

        else
            % no clustering required
            R(i,j) = adjusted_rand_index(TrackLabels,Labels);
            Clusters_Track_Constraint(i,j) = numtracks;
        end


    end

    fprintf('Constrained track centroids,  Video %i %s done.\n',i,video)

end
E = R;

save Res_ConstrainedE

%--------------------------------------------------------------------------



% (F) CBC -----------------------------------------------------------------

R = [];

for i = 11:numel(Videos)

    video = Videos{i};
    % Data filename
    ofnRGB = [VideoPath,video,'_RGB.mat'];
    load(ofnRGB)
    DataRGB = zscore(DataRGB);

    for j = 1:3
        filename_tr_ML = [TrackPaths{j},video,'.csv'];
        M = readmatrix(filename_tr_ML);
        TrackLabels = M(:,2);
        Frames = M(:,1);

        c = max(Labels); % true number of clusters


        % find track centroids
        u = unique(TrackLabels);
        numtracks = numel(u); % number of tracks

        if c < numtracks
            AL = classifier_clustering(DataRGB, Frames, ...
                c, TrackLabels, Labels,1);

            R(i,j) = adjusted_rand_index(AL,Labels);
            Clusters_CBC(i,j) = numel(unique(AL));
        else

            % no clustering required
            R(i,j) = adjusted_rand_index(TrackLabels,Labels);
            Clusters_Track_Constraint(i,j) = numtracks;
        end
    end
    fprintf('CBC,  Video %i %s done.\n',i,video)
end

F = R;

save Res_ConstrainedF

% FUNCTIONS ---------------------------------------------------------------
%**************************************************************************
%--------------------------------------------------------------------------

function assigned_labels = run_linkage(data,c,mode)

% mode = choice of linkage method from 1 to 7

LinkageMethods = {'single', 'complete', 'average','weighted', ...
    'centroid', 'median', 'ward'};

Z = linkage(data,LinkageMethods{mode});
assigned_labels = cluster(Z,maxclust=c);

end
%--------------------------------------------------------------------------

function assigned_labels = run_kmeans(data,c)
assigned_labels = kmeans(data,c);
end
%--------------------------------------------------------------------------

function assigned_labels = run_gmm(data,c)

try
    options = statset('MaxIter',1000);
    gmfit = fitgmdist(data,c,'CovarianceType','diagonal', ...
        'SharedCovariance',true,'Options',options); % Fitted GMM
    assigned_labels = cluster(gmfit,data); % Cluster index
catch
    assigned_labels = ones(size(data,1),1);
end

end
%--------------------------------------------------------------------------

function assigned_labels = run_dbscan(data,c)


E = 0.5:0.5:4; % epsilon
Mi = 1:8; % minpts

for i = 1:numel(E)
    for j = 1:numel(Mi)
        AL = dbscan(data,E(i),Mi(j));
        C(i,j) = numel(unique(AL));
    end
end

% Find the closest to the desired number of classes
[~,index_cl] = min(abs(C(:)-c));
[ibest,jbest] = ind2sub([numel(E),numel(Mi)],index_cl);
assigned_labels = dbscan(data,E(ibest),Mi(jbest));

end
%--------------------------------------------------------------------------

function assigned_labels = run_spectral(data,c)

% PATCH to eliminate spurious error
flag = 1;
while flag
    try
        flag = 0;
        assigned_labels = spectralcluster(data,c);
    catch
        flag = 1;
        c = c - 1;
        if c == 2
            flag = 0; % end loop
        end
    end
end
% END OF PATCH
end
%--------------------------------------------------------------------------


function assigned_labels = run_finch(data,c)

[AL, num_clust]= FINCH(data,[], 0);
% May return several clustering results with different number of clusters.
% We will use the number of clusters closer to c.

[~,chosen] = min(abs(num_clust-c));
assigned_labels = AL(:,chosen);

end
%--------------------------------------------------------------------------

function [CL,ML] = find_all_constraints(Frames,Labels)

uu = unique(Frames);
CL = [];

for i = 1:numel(uu)
    Framesi = find(Frames == uu(i)); % objects in ith frame
    for j1 = 1:numel(Framesi)-1
        for j2 = j1+1:numel(Framesi)
            CL = [CL;Framesi(j1),Framesi(j2)];
        end
    end
end

% --
u = unique(Labels);
ntr = numel(u);

ML = [];
for i = 1:ntr % for each track/label
    Tracki = find(Labels == u(i)); % objects for track i
    % We shall assume that all objects in the track are linked timewise
    % (only for our Must Link approach!)

    % Niote that it is enough to generate one link for each object. The
    % transitive closure will take care of the other links.

    for j = 1:numel(Tracki)-1
        ML = [ML;Tracki(j), Tracki(j+1)];
    end

end

end

%--------------------------------------------------------------------------

function CL = find_cl_tracks(Frames,Labels)

u = unique(Labels);
ntr = numel(u);
CL = [];
for i = 1:ntr-1
    Framesi = Frames(Labels == u(i));

    for j = i+1:ntr
        Framesj = Frames(Labels == u(j));
        if ~isempty(intersect(Framesi,Framesj))
            CL = [CL;i, j];
        end
    end

end

end

%--------------------------------------------------------------------------