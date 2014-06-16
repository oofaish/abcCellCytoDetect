%close all;
params = abcParams( );
ts = [];
ls = [];

first = struct();
final = struct();

x = [ 0, 0, pi /3, 1];

gBlur = fspecial('gaussian',[10 10],0.5);  
cellOnlyImageCopy = imfilter(cellOnlyImage,gBlur,'same');


[ final.canvas, final.m, final.n ] = abcTransformCell( cellOnlyImageCopy, cellOnlyMask, nucleusOnlyMask, x, params );

first =final;

for theta = 0:0.1:pi * 2
    x = [ 0, 0, theta, 1];
    [ final.canvas, final.m, final.n ] = abcTransformCell( cellOnlyImageCopy, cellOnlyMask, nucleusOnlyMask, x, params );
    
    ts( end + 1 ) = theta;
    ls( end + 1 ) = abcLossFnSimple( final, first, op );                       
    disp( theta );
    %imshow( final.c );
    plot( ts, ls, 'k' );
    pause( 0.01 );
end

plot( ts, ls, 'k' );