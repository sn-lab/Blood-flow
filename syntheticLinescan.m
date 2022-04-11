% TODO: make (-/+) 75 and 45 degree data
%%
I = zeros(3072, 3072, 'uint16');
%%
stripSize = 4;
for i = stripSize:stripSize*2:3072
    I(i:i+stripSize, :) = 32767;
end

%%
angle = 15;
I = imrotate(I, angle, 'bilinear', 'crop');
I = I(1025:2048, 1025:2048);

%%
figure; imshow(I, [])