clear, clc, close all
T = readmatrix('Data_concentric_circles_3_100.csv');
d = T(:,[1,2]);
l = T(:,3);

figure
scatter(d(:,1),d(:,2),20,l,'filled')
axis equal

[c,n] = FINCH(d,[],1);

figure
scatter(d(:,1),d(:,2),20,c(:,end),'filled')
axis equal
colormap jet