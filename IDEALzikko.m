
close all;
A=reshape(ans(1,:,:),32,32);
B=reshape(ans(2,:,:),32,32);
C=reshape(ans(3,:,:),32,32);
[chem]=IDEAL(5.1, 2.55, 0, 2, 12, A, B, C);
% varargin{1} = ta
% varargin{2} = delta_TE
% varargin{3} = ppm
% varargin{4} = chem_num
% varargin{5} 〜 varargin{4+ppm_num} = ppm_diff
% varargin{5+ppm_num} 〜 varargin{nargin} = image ( > chem_num)
% 3画像2分子分離の入力例 [chem]=IDEAL(5.0 (最初のTE), 2.5 (TEの間隔), 0 (中心周波数が合ってたら0), 2 (分子数), 12 (ppmの差), A (Teが最小のの画像), B (TEが2番目の画像), C(TEが最大の画像))

m1=reshape(chem(1,:,:),32,32);
m2=reshape(chem(2,:,:),32,32);

figure
I=imrotate(abs(m1),90);
imagesc(abs(m1))
figure
J=imrotate(abs(m2),90);
imagesc(abs(m2))

figure
imagesc(abs(A))
%figure
%imagesc(angle(A))
%figure
%imagesc(angle(B))
%figure
%imagesc(angle(C))