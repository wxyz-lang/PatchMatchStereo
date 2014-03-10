img = imread('1.png');
img = medfilt2(img, [3 3]);
imwrite(img, '1_medfilt.png');

img = imread('2.png');
img = medfilt2(img, [3 3]);
imwrite(img, '2_medfilt.png');

img = imread('3.png');
img = medfilt2(img, [3 3]);
imwrite(img, '3_medfilt.png');