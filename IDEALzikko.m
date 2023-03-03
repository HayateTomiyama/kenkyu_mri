
close all;
A=reshape(ans(1,:,:),32,32);
B=reshape(ans(2,:,:),32,32);
C=reshape(ans(3,:,:),32,32);
[chem]=IDEAL(5.1, 2.55, 0, 2, 12, A, B, C);
% varargin{1} = ta
% varargin{2} = delta_TE
% varargin{3} = ppm
% varargin{4} = chem_num
% varargin{5} ` varargin{4+ppm_num} = ppm_diff
% varargin{5+ppm_num} ` varargin{nargin} = image ( > chem_num)
% 3‰æ‘œ2•ªq•ª—£‚Ì“ü—Í—á [chem]=IDEAL(5.0 (Å‰‚ÌTE), 2.5 (TE‚ÌŠÔŠu), 0 (’†Sü”g”‚ª‡‚Á‚Ä‚½‚ç0), 2 (•ªq”), 12 (ppm‚Ì·), A (Te‚ªÅ¬‚Ì‚Ì‰æ‘œ), B (TE‚ª2”Ô–Ú‚Ì‰æ‘œ), C(TE‚ªÅ‘å‚Ì‰æ‘œ))

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