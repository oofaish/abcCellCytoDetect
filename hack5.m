close all;
I = target.c;
%I = rgb2gray( RGB );

% Extract edges.
BW = edge(I,'canny');
[H,T,R] = hough(BW,'RhoResolution',0.5,'Theta',-90:0.5:89.5);

% Display the original image.
subplot(2,1,1);
imshow(BW);
title('Gantrycrane Image');

% Display the Hough matrix.
subplot(2,1,2);
imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
      'InitialMagnification','fit');
title('Hough Transform of Gantrycrane Image');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot);