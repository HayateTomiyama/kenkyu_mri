%2���q�V�~�����[�V�����̎��s�͂���
%function [chem,A,B,C]=simu2(MATRIX,P,Q,ppm_diff_P,ppm_diff_Q,ta,delta_t,initial_phase,delta_B_rate,abs_noise_max)
%P,Q�͋��x�@�Z * 10^4 �̊ۂ̕��������
%ppm_diff_�Z��   ���q�Z�̒��S���g������̂���
%ta 1�ڂ�TE
%delta_B_rate ����s�ψ�̂���̊����@0.01�Ƃ�
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

