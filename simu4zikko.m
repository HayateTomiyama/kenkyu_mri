%4分子分離シミュレーションの時はこれ


[chem]=simu4(32,5,4,3,2,0,20.2,3.8,12.7,7,1.9,0,0.01,0.01);
A=reshape(chem(1,:,:),32,32);
B=reshape(chem(2,:,:),32,32);
C=reshape(chem(3,:,:),32,32);
D=reshape(chem(4,:,:),32,32);

t = tiledlayout(2,2,'TileSpacing','Compact');

% Tile 1
nexttile
imagesc(abs(A))
axis image
colorbar
caxis([0 50000])
title('A')
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

% Tile 2
nexttile
imagesc(abs(B))
axis image
colorbar
caxis([0 50000])
title('B')
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';


% Tile 3
nexttile
imagesc(abs(C))
axis image
colorbar
caxis([0 50000])
title('C')
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';


% Tile 4
nexttile
imagesc(abs(D))
axis image
colorbar
caxis([0 50000])
title('D')
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
