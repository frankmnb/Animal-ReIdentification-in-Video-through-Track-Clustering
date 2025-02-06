function plot_realtime(points, labels, c, k)

    % hold(ax,'on')
    % hold on

    u = numel(unique(labels));
    markers = '+ox^v><dshp*';
    sizemark = length(markers);
    rep = ceil(u/sizemark);
    newM = repmat(markers,1,rep); newM = newM(1:u);

    if nargin == 3
        k = 10;
    end

    gscatter(points(:,1),points(:,2),labels,c,newM,k)