function  assembled_labels = reassemble_labels(AL_raw,TrackLabels)
% Distribute the labels of the super-clusters to the original objects

assembled_labels = zeros(size(TrackLabels));
u = unique(TrackLabels); % This should be the number of elements of AL_raw

uu = unique(AL_raw);
for j = 1:numel(uu)
    SuperJ = find(AL_raw == uu(j));

    for k = 1:numel(SuperJ)
        % all these must be relabelled as j
        % the assigned labels will be 1, 2, ...

        assembled_labels(TrackLabels == u(SuperJ(k))) = j;
    end

end
end