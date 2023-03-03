function [image,mask]=noiseeliminator(MATRIX_x,MATRIX_y,image_num,image)

mask = ones(MATRIX_y,MATRIX_x);

A=reshape(image(1,:,:),MATRIX_y,MATRIX_x);

for k = 1:MATRIX_y
    for l = 1:MATRIX_x
        if abs(A(k,l)) < 0.2*max(abs(A),[],'all')
            mask(k,l) = 0;
            for m = 1:image_num
                image(m,k,l) = 0;
            end
        end
    end
end

end