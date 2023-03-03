function [chem]=IDEAL(varargin)
% varargin{1} = ta, varargin{2} = delta_TE, varargin{3} = ppm
% varargin{4} = chem_num
% varargin{5} ` varargin{4+ppm_num} = ppm_diff
% varargin{5+ppm_num} ` varargin{nargin} = image ( > chem_num)
% 3‰æ‘œ2•ªq•ª—£‚Ì“ü—Í—á [chem]=IDEAL(5.0 (Å‰‚ÌTE), 2.5 (TE‚ÌŠÔŠu), 0 (’†Sü”g”‚ª‡‚Á‚Ä‚½‚ç0), 2 (•ªq”), 12 (ppm‚Ì·), A (Te‚ªÅ¬‚Ì‚Ì‰æ‘œ), B (TE‚ª2”Ô–Ú‚Ì‰æ‘œ), C(TE‚ªÅ‘å‚Ì‰æ‘œ))

sz = size(varargin{nargin});
MATRIX_x = sz(1,2);
MATRIX_y = sz(1,1);
ta = varargin{1};
delta_TE = varargin{2};
ppm = varargin{3};
chem_num = varargin{4};
ppm_num = chem_num-1;
image_num = nargin-ppm_num-4;

if (image_num <= chem_num)
    disp('error, ‰æ‘œ”‚ª‘«‚è‚Ä‚¢‚È‚¢“™‚Ì‰Â”\«‚ ‚è');
    return
end

image = zeros(image_num,MATRIX_y,MATRIX_x);
for i = 1:image_num
    image(i,:,:) = varargin{i+ppm_num+4};
end

[image,mask] = noiseeliminator(MATRIX_x,MATRIX_y,image_num,image);

rot = zeros(1,ppm_num);

for k = 1:ppm_num
    rot(1,k) = 2*pi*16.342*varargin{4+k};
end

TE = zeros(1,image_num);
initial = zeros(1,image_num);

for k = 1:image_num
    TE(1,k) = 0.001*(ta+(k-1)*delta_TE);
    initial(1,k) = exp(-2i*pi*16.342*ppm*TE(1,k));
    image(k,:,:) = image(k,:,:).*initial(1,k);
end

P = eye(image_num);
phi = zeros(MATRIX_y,MATRIX_x);
Q = zeros(image_num,chem_num);

for k = 1:image_num
    for l = 1:chem_num
        if (l == 1)
            Q(k,l) = 1;
        else
            Q(k,l) = exp(1i*rot(1,l-1)*TE(1,k));
        end
    end
end

R=Q.';
S=zeros(image_num,1);
T = [ones(image_num,1) Q];

for k = 1:MATRIX_y
    for l = 1:MATRIX_x
        if (mask(k,l) == 1)
            RO_dash = 2;
            w=0;
            for m = 1:image_num
                S(m,1) = image(m,k,l);
            end
            while (abs(real(RO_dash(1,1))) > 1 || 2*pi*abs(imag(RO_dash(1,1))) > 1) && (w < 2000)
                if w == 0
                    klphi = 0;
                end
                for m = 1:image_num
                    P(m,m) = exp(-2i*pi*klphi*TE(1,m));
                end
                RO = ((R*Q)^(-1))*R*P*S;
                comp = RO(1,1);
                for m = 1:ppm_num
                    comp = comp + RO(m+1,1)*exp(1i*pi*rot(1,m));
                end
                for m = 1:image_num
                    T(m,1) = 2i*pi*TE(1,m)*comp;
                end
                U = T.';
                RO_dash = ((U*T)^(-1))*U*(P*S-Q*RO);
                klphi = klphi + RO_dash(1,1);
                klphi = ((real(klphi)*delta_TE*0.001-round(real(klphi)*delta_TE*0.001))*1000/delta_TE)+1i*imag(klphi);
                w = w + 1;
            end
            phi(k,l) = klphi;
        end
    end
end

smoothed_phi = boxcar_average(phi,MATRIX_y,MATRIX_x);

chem = zeros(chem_num,MATRIX_y,MATRIX_x);

W = zeros(image_num,image_num);

for k = 1:MATRIX_y
    for l = 1:MATRIX_x
            for m = 1:image_num
                S(m,1) = image(m,k,l)*exp(-2i*pi*real(smoothed_phi(k,l))*TE(1,m));
                W(m,m) = exp(2*imag(smoothed_phi(k,l))*TE(1,m));
            end
            RO = ((R*W*Q)^(-1))*R*W*S;
            for m = 1:chem_num
                chem(m,k,l) = RO(m,1);
            end
    end
end

chem = noiseeliminator_z(chem,chem_num,MATRIX_x,MATRIX_y);

end

function [image,mask]=noiseeliminator(MATRIX_x,MATRIX_y,image_num,image)

mask = ones(MATRIX_y,MATRIX_x);

A=reshape(image(1,:,:),MATRIX_y,MATRIX_x);

for k = 1:MATRIX_y
    for l = 1:MATRIX_x
        if abs(A(k,l)) < 0.08*max(abs(A),[],'all')
            mask(k,l) = 0;
            for m = 1:image_num
                image(m,k,l) = 0;
            end
        end
    end
end
end

function chem = noiseeliminator_z(chem,chem_num,MATRIX_x,MATRIX_y)

for m =1:chem_num
    A = reshape(chem(m,:,:),MATRIX_y,MATRIX_x);
    for k = 1:MATRIX_y
        for l = 1:MATRIX_x
            if abs(A(k,l)) < 0.2*max(abs(A),[],'all')
                chem(m,k,l) = 0;
            end     
        end
    end
end    

end

function [F]=boxcar_average(A,MATRIX_x,MATRIX_y)

F=zeros(MATRIX_y,MATRIX_x);

for i=3:MATRIX_y-2
    for j=3:MATRIX_x-2
        Q=0;
        for p=i-2:i+2
            for q=j-2:j+2
                Q=Q+A(p,q);
            end
        end
        F(i,j)=Q/25;
    end
end
end