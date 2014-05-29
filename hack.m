close all

cellToTest = -1;%if this is set a specific cell, we dont create a new image
                 %just work on the last image.
testingOneCell = ( cellToTest > 0 );

if exist( 'abcGenerateImage', 'file' ) ~= 2
    addpath( '../abcCellGenerator/' );
    addpath( '../fminsearchbnd/' );
end

%generate a random cells image
if( 0 && ~testingOneCell )
    params = struct();
    params.totalNumberOfCells = 10;
    params.randomiseAlpha = false;
    params = abcParams( params );
    [ canvas, nucleusInfos, cellInfos ] = abcGenerateImage( false, false, '', params );
end

%go through each detected nucleus, and search for ideal parameters
 
%Set the parameters you are changing (currently only radius)
pSpace = struct();
pSpace.radius = 10;
pSpace.majorVsMinor = 1;
pSpace.majorVsMinorAngle =  0;
 
cellParams = abcCellParams();
%FIXME - actually implement this. Maybe you should use the nucleus mask to
%zero out those parts of the edge (but then what if you zero out real edges
%by accident?
cellParams.drawNucleus = false;
 
optimisationParams = struct();
optimisationParams.minRadius = cellParams.nucleusRadius * 2;
optimisationParams.maxRadius = 100;
%optimisationParams.maxAngle = 0;
%optimisationParams.minAngle = pi;

optimisationParams.checkUniformIntensity = false;%makes no difference

optimisationParams.x0         = [ optimisationParams.minRadius * 2, 1.3, pi/4 ];
%optimisationParams.x0         = [ 60, 1.34, 2.0364 / 2 ];
optimisationParams.lowerBound = [ optimisationParams.minRadius, 1, 0 ];
optimisationParams.upperBound = [ optimisationParams.maxRadius, 3, pi ];

optimisationParams.visuallyVerbose = false && testingOneCell;
%optimisationParams.algorithm = 'gridsearch';
optimisationParams.algorithm = 'simplex';
%optimisationParams.lossFunction = 'exactNonZeros';
optimisationParams.lossFunction = 'truncatedSquareLoss';
optimisationParams.truncation = optimisationParams.maxRadius;
optimisationParams.epsilon = 1e-6;
    

h1 = figure;
if( true || ~testingOneCell )
    subplot( 121 );
    imshow( canvas );
    pause( 0.01 );
end

%h2 = subplot( 122 );
%hold on;
newMask = canvas * 0;
for i = 1:numel( nucleusInfos )
    if testingOneCell && i ~= cellToTest
        continue
    end
    
    optimisationParams.actualRadius = cellInfos{ i }.pSpace.radius;
    [ pSpace, score, cellMask ] = abcMinLossForCell( pSpace, cellParams, canvas, nucleusInfos{i}, cellInfos{ i }, optimisationParams );
    
    abcFancyPrintCellResults( i, score, pSpace, cellInfos{ i }.pSpace );
    
    newMask( cellMask == 1 ) = 1;
    subplot( 122 );
    imshow( newMask );
    pause( 0.01 );
end
 
 subplot( 122 );
 imshow( newMask );
 
 
 hold off;