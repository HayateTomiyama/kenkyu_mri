%2分子シミュレーションの実行はこれ
%function [chem,A,B,C]=simu2(MATRIX,P,Q,ppm_diff_P,ppm_diff_Q,ta,delta_t,initial_phase,delta_B_rate,abs_noise_max)
%P,Qは強度　〇 * 10^4 の丸の部分を入力
%ppm_diff_〇は   分子〇の中心周波数からのずれ
%ta 1個目のTE
%delta_B_rate 磁場不均一のずれの割合　0.01とか
close all;

[chem]=simu2(32,5,3,0,12,5.1,2.55,0,0.01,10);
A=reshape(chem(1,:,:),32,32);
B=reshape(chem(2,:,:),32,32);


t = tiledlayout(1,2,'TileSpacing','Compact');

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

