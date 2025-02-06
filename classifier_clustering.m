function [AL,zz1,zz2,zz3,zz4] = ...
    classifier_clustering(data, Frames, c, TrackLabels, Labels, flag)

utr = unique(TrackLabels);
ntr = numel(utr); % number of tracks
% N = size(data,1); % number of objects

if nargin == 5
    flag = 0;
end

CL = find_cl_tracks(Frames,TrackLabels);
fprintf('Number of CL constraints: %i\n',size(CL,1))
fprintf('Data size: %i %i\n',size(data))
fprintf('Number of tracks: %i\n',ntr)
fprintf('Desired number of clusters: %i\n',c)


AL = TrackLabels;
zz1 = 1; zz2 = adjusted_rand_index(AL,Labels);
zz3 = 1;
zz4 = ntr;

thr_resub = 0.7;

if ntr < c
    error('Too few tracks')
end

% iter = 1;
% while iter <= 20  % Combine tracks until there are c clusters
%     iter = iter + 1;
while ntr > c % Combine tracks until there are c clusters
    utr = unique(TrackLabels);
    ntr = numel(utr); % number of tracks

    if flag
        fprintf('%i clusters, %.4f\n',ntr,zz2(end))
    end

    clf = fitcdiscr(data,TrackLabels,"DiscrimType","diaglinear");

    cm = confusionmat(TrackLabels,predict(clf,data)); % re-sub

    % Mask the CL constraints
    cm(sub2ind([ntr,ntr],CL(:,1),CL(:,2))) = 0;
    cm(sub2ind([ntr,ntr],CL(:,2),CL(:,1))) = 0;

    % Rescale the confusion matrix to reflect which clusters are most
    % assigned to another
    cm = cm./repmat(sum(cm,2),1,ntr);

    cm = cm.*(1-eye(ntr)); % mask the diagonal

    [~,indmax] = max(cm(:)); % find the most confused classes (presumably
    % the same class label identity)
    [c1,c2] = ind2sub([ntr,ntr],indmax); % clusters to merge

    % Check if merge is possible and if so, - do it
    if cm(c1,c2) > 0 % no CL constraint
        TrackLabels(TrackLabels == utr(c1)) = utr(c2);
        % CL = find_cl_tracks(Frames,TrackLabels);
        CL = update_cl(CL,c1,c2);
        zz1 = [zz1;cm(c1,c2)]; %#ok<*AGROW>
        zz2 = [zz2,adjusted_rand_index(TrackLabels,Labels)];
        zz3 = [zz3,mean(AL==TrackLabels)]; % resubstitution error
        % if zz3(end) < thr_resub
        %     ntr = c;
        % end
        zz4 = [zz4, ntr];
    else
        ntr = c;
    end
end

AL = TrackLabels;
fprintf('Obtained number of clusters: %i\n',numel(unique(AL)))

end

function CL = update_cl(CL,c1,c2)
% Checked with random numbers
CL(CL == c1) = c2;
e = unique(sort(CL,2),'rows');
e(e>c1) =  e(e>c1) - 1;
CL = e;

end