rng('shuffle');

close all;

cellToTest = -1;%if this is set a specific cell, we dont create a new image
                 %just work on the last image.
testingOneCell = ( cellToTest > 0 );

if exist( 'abcGenerateImage', 'file' ) ~= 2
    addpath( '../abcCellGenerator/' );
    addpath( '../fminsearchbnd/' );
end

regenerateImage = false;

%generate a random cells image
if( regenerateImage || ~exist( 'canvas', 'var' ) )
    params = struct();
    params.totalNumberOfCells = 5;
    params.randomiseAlpha = false;
    params = abcParams( params );
    [ canvas, nucleusInfos, cellInfos ] = abcGenerateImage( false, false, '', params );
    
    if( 1 )
        gBlur = fspecial('gaussian',[3 3],2);    
        canvas = imfilter(canvas,gBlur,'same');
    end

    
end

optimisationParams = abcOptimisationParams();

[ splits, masks ] = abcSplitImage( nucleusInfos, optimisationParams, params.canvasSize );

assert( size( splits, 2 ) == max( max( splits ) ), 'Something went wrong, I missed some cells in the split' );

pSpace = struct();
pSpace.radius = 10;
pSpace.majorVsMinor = 1;
pSpace.majorVsMinorAngle =  0;
 
cellParams = abcCellParams();
%FIXME - actually implement this. Maybe you should use the nucleus mask to
%zero out those parts of the edge (but then what if you zero out real edges
%by accident?
cellParams.drawNucleus = false;
 
optimisationParams = abcOptimisationParams();

%imshow( canvas );
scrsz = get(0,'ScreenSize');
h1 = figure('Position',[0 scrsz( 4 ) / 1.5 scrsz(3) scrsz(4) * 3 / 5]);
if( 1 || testingOneCell )
    subplot( 131 );
    imshow( canvas );
    pause( 0.01 );
end

%h2 = subplot( 122 );
%hold on;
newMask = canvas * 0;
newCanvas  = abcEmptyCanvas( params.canvasSize,  true );

for i = 1:size( splits, 1 )
    
    if testingOneCell && i ~= cellToTest
        continue
    end

    
    cells = splits(i, :);
    cells = cells( cells ~= 0 );
    
    %if testingOneCell && i ~= cellToTest
    %    continue
    %end

    [ x, minLoss, thisCanvas, mask ] = abcMinLoss( cells, cellParams, params,  canvas, nucleusInfos, cellInfos, optimisationParams, masks{ i } );    
    %[ pSpace, score, cellMask ] = abcMinLossForCell( pSpace, cellParams, canvas, nucleusInfos{i}, cellInfos{ i }, optimisationParams );
    
     newCanvas = abcAddNextCellToCanvas( newCanvas, thisCanvas );
     numberOfCells = numel( cells );
     actualX = x * 0;
     if( 1 )
         for i = 0:(numberOfCells-1)
            actualX( i * 3 + 1 ) = cellInfos{ cells( i + 1 ) }.pSpace.radius;
            actualX( i * 3 + 2 ) = cellInfos{ cells( i + 1 ) }.pSpace.majorVsMinor;
            actualX( i * 3 + 3 ) = cellInfos{ cells( i + 1 ) }.pSpace.majorVsMinorAngle;
         end
     else
         for i = 1:numberOfCells
            actualX( i ) = cellInfos{ cells( i ) }.pSpace.radius;
         end
         
     end

    abcFancyPrintCellResults( i, minLoss, x, actualX );
    
    newMask( mask == 1 ) = 1;
    subplot( 132 );
    imshow( newCanvas );
    pause( 0.01 );
end
  
 
 hold off;