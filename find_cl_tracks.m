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