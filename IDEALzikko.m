
close all;
A=reshape(ans(1,:,:),32,32);
B=reshape(ans(2,:,:),32,32);
C=reshape(ans(3,:,:),32,32);
[chem]=IDEAL(5.1, 2.55, 0, 2, 12, A, B, C);
% varargin{1} = ta
% varargin{2} = delta_TE
% varargin{3} = ppm
% varargin{4} = chem_num
% varargin{5} �` varargin{4+ppm_num} = ppm_diff
% varargin{5+ppm_num} �` varargin{nargin} = image ( > chem_num)
% 3�摜2���q�����̓��͗� [chem]=IDEAL(5.0 (�ŏ���TE), 2.5 (TE�̊Ԋu), 0 (���S���g���������Ă���0), 2 (���q��), 12 (ppm�̍�), A (Te���ŏ��̂̉摜), B (TE��2�Ԗڂ̉摜), C(TE���ő�̉摜))

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