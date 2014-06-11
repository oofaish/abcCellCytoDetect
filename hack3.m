if exist( 'homoSpaceTransform', 'file' ) ~= 2
    addpath( '/Users/oofaish/Projects/Matlab/cellGenerator/Funs' );
end
if exist( 'Free_lyingCellMaskSet', 'var' ) ~= 1
    load('/Users/oofaish/Projects/Matlab/cellGenerator/N3_R_040_050/Variables/EDF_AllFreelyingCellPool.mat');
    disp( 'Loading variables...' );
end

if exist( 'abcGenerateImage', 'file' ) ~= 2
    addpath( '../abcCellGenerator/' );
end


params = abcParams( );
%FIXME - doing a smaller canvas for now.
params.canvasSize = [ 250, 250 ];
i = 1;

%these are the smaller images, so not on canvas
cellOnlyImage   = Free_lyingCellSet{i};
cellOnlyMask    = Free_lyingCellMaskSet{i};
nucleusOnlyMask = Free_lyingCellNucleiMaskSet{i};

%newCanvas  = abcEmptyCanvas( params.canvasSize,  true );

x = [ 30, -50, pi / 4, 1.1 ];

original = struct();
new = struct();


if( 0 )
    [ original.c, original.m, original.n ] = abcTransformCell( cellOnlyImage, cellOnlyMask, nucleusOnlyMask, [ 0, 0, 0, 1 ], params );
    [ new.c, new.m, new.n ] = abcTransformCell( cellOnlyImage, cellOnlyMask, nucleusOnlyMask, x, params );

    close all;
    figure
    subplot( 121 );
    imshow( original.c );
    subplot( 122 );
    imshow( new.c );
end



