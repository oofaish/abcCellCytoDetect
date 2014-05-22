close all

cellToTest = -1;%if this is set a specific cell, we dont create a new image
                 %just work on the last image.
testingOneCell = ( cellToTest > 0 );

if exist( 'abcGenerateImage', 'file' ) ~= 2
    addpath( '../abcCellGenerator/' );
end

%generate a random cells image
if( ~testingOneCell )
    params = struct();
    params.totalNumberOfCells = 20;
    params.randomiseAlpha = true;
    params = abcParams( params );
    [ canvas, nucleusInfos, cellInfos ] = abcGenerateImage( false, false, '', params );
end

%go through each detected nucleus, and search for ideal parameters
 
%create a structure with parameters that are staying fixed
%Set the parameters you are changing (currently only radius)
 
pSpace = struct();
pSpace.radius = 10;
 
cellParams = abcCellParams();
%FIXME - actually implement this. Maybe you should use the nucleus mask to
%zero out those parts of the edge (but then what if you zero out real edges
%by accident?
cellParams.drawNucleus = false;
 
optimisationParams = struct();
optimisationParams.minRadius = cellParams.nucleusRadius * 2;
optimisationParams.x0        = optimisationParams.minRadius * 3;
optimisationParams.maxRadius = 100;
optimisationParams.visuallyVerbose = false || testingOCell;
%optimisationParams.algorithm = 'gridsearch';
optimisationParams.algorithm = 'simplex';
%optimisationParams.lossFunction = 'exactNonZeros';
optimisationParams.lossFunction = 'truncatedSquareLoss';
optimisationParams.truncation = optimisationParams.maxRadius;
    

h1 = figure;
if( ~testingOneCell )
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
    
    disp( [ 'Cell: ', num2str( i ) ] );
    disp( [ 'Cell: ', num2str( i ), ', Actual Radius= ', num2str( cellInfos{ i }.pSpace.radius ), ', Found Radius= ', num2str( pSpace.radius ) ] );
    newMask( cellMask == 1 ) = 1;
    subplot( 122 );
    imshow( newMask );
    pause( 0.01 );
end
 
 subplot( 122 );
 imshow( newMask );
 
 
 hold off;