%close all;
params = abcParams( );
original = abcTransformCell( cellOnlyImage, cellOnlyMask, nucleusOnlyMask, [ 0, 0, 0, 1 ], params, [] );
stat = regionprops( double( original.n ), 'Centroid');
%nCentroid = [ floor( stat(1,1).Centroid(1,2) ), floor( stat(1,1).Centroid(1,1) ) ] + [ 150, 400 ] - [ 256 256 ];
nCentroid = [ round( stat(1,1).Centroid(1,2) ), round( stat(1,1).Centroid(1,1) ) ];


for theta = 0:0.1:pi
    x = [ 0, 0, theta, 1];
    final = abcTransformCell( cellOnlyImage, cellOnlyMask, nucleusOnlyMask, x, params, nCentroid );

    final.c( nCentroid( 1 ) - 5:1:nCentroid( 1 ) + 5, nCentroid( 2 ) - 5:1:nCentroid( 2 ) + 5 ) = 1;
    imshow( final.c );
    
    pause( 0.01 );
end