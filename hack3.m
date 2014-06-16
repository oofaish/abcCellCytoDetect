close all;

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

saveMovie = false;

global writerObj;
if saveMovie
    writerObj = VideoWriter( strcat( movieName, '.avi' ) );
else
    writerObj = -1;
end

if( writerObj ~= -1 )
    open(writerObj);
end



params = abcParams( );
%FIXME - doing a smaller canvas for now.
%params.canvasSize = [ 250, 250 ];
i = 7;

%these are the smaller images, so not on canvas
allCellInfos   = struct();
allCellInfos.c = Free_lyingCellSet;
allCellInfos.m = Free_lyingCellMaskSet;
allCellInfos.n = Free_lyingCellNucleiMaskSet;

cellOnlyImage   = allCellInfos.c{ i };
cellOnlyMask    = allCellInfos.m{ i };
nucleusOnlyMask = allCellInfos.n{ i };

if( 0 )
    gBlur = fspecial('gaussian',[3 3],2);    
    cellOnlyImage = imfilter(cellOnlyImage,gBlur,'same');
end


%newCanvas  = abcEmptyCanvas( params.canvasSize,  true );

%FIXME - what is the origin of rotation?
xTarget = [ 0, 0,  pi * 9 / 10, 1.5 ];

op = struct();
%op.x0 = [ 0, 0, 0, 1 ]; 
op.x0 = [ 0, 1 ];
op.xMin = [ 0, 0.8 ];
op.xMax = [ 2 * pi, 1.2 ];
op.xStep = [ 0.2, 0.2 ];
op.c = cellOnlyImage;
op.m = cellOnlyMask;
op.n = nucleusOnlyMask;
op.FinDiffRelStep = [ 0.01, 0.1 ]; %gradstep
op.diffMaxChange = 1;
op.diffMinChange = 0.4;
op.method = 'gridSearch';
op.lossMethod = 'double';
%op.method = 'lm';

hCounter = 1:params.canvasSize( 1 );
wCounter =  1:params.canvasSize( 2 );
    
%maxAllowedDist = 17;
maxAllowedDist = 170;
coverageMask   = bsxfun( @plus, ( hCounter - 256 ).^2, ( wCounter' - 256 ).^2 ) < maxAllowedDist^2;

original = abcTransformCell( cellOnlyImage, cellOnlyMask, nucleusOnlyMask, [ 0, 0, 0, 1 ], params, [] );

%get the target and the center of the nucleus (which will actually come from Carlos
target   = abcTransformCell( cellOnlyImage, cellOnlyMask, nucleusOnlyMask, xTarget, params, [] );
stat = regionprops( target.n, 'Centroid');
target.nCentroid = [ round( stat(1,1).Centroid(1,2) ), round( stat(1,1).Centroid(1,1) ) ];

target.coverageMask = coverageMask;

if( 1 )
    disp( 'Gridsearch First' );
    [ x2, minLoss ] = abcMinLossHack( target, params, op );
end
op.x0 = [ x2( 3 ) x2( 4 ) ];
disp( 'Now LM' );
op.method = 'lm';
op.lossMethod = 'vector';

[ x3, minLoss ] = abcMinLossHack( target, params, op );

disp( xTarget );
disp( x3 );


if( 1 )
    
    final = abcTransformCell( cellOnlyImage, cellOnlyMask, nucleusOnlyMask, x3, params, target.nCentroid );
    close all;
    scrsz = get(0,'ScreenSize');
    h1 = figure('Position',[0 scrsz( 4 ) / 1.5 scrsz(3) scrsz(4) * 3 / 5]);

    subplot( 131 );
    imshow( target.c );
    subplot( 132 );
    imshow( final.c );
    
    canvas       = double( final.c ) .* target.coverageMask;
    targetCanvas = double( target.c ) .* target.coverageMask;
    tmp0 = ( targetCanvas - canvas ) / 255.0;
    
    subplot( 133 );
    imshow( ( tmp0 + 1 ) ./ 2 );
    
end



