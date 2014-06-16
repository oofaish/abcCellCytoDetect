function [ tx, ty ] = abcTranslationForCellFromNucleusPosition( nCentroid, cellId, theta, delta, allCellInfos, params );
%need a function, when I give it:

%- the position of the nucleus (from Carlos's code
%- cell number
%- cell rotation and scale

%it gives me an appropriate x and y
    x = [ 0, 0, theta, delta ];
    
    final = abcTransformCell( allCellInfos.c{ cellId }, allCellInfos.m{ cellId }, allCellInfos.n{ cellId }, x, params );
    
    stat = regionprops(final.n, 'Centroid');

    tx = nCentroid( 1 ) - stat(1,1).Centroid(1,2);
    ty = nCentroid( 2 ) - stat(1,1).Centroid(1,1);
end