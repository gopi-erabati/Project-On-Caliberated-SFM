function img = createImage(points)
% a function to create image from 2d points

img = zeros(512,512);

for i= 1:size(points,2)
    img(ceil(points(2,i)), ceil(points(1,i))) = 255;
end

end



