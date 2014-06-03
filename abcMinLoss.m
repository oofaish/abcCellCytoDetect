function [ x, minLoss, canvas, mask ] = abcMinLoss( cells, cellParams, allParams,  canvas, nucleusInfos, cellInfos, optimisationParams, coverageMask )
%ABCMINLOSS change parameters in pSpace to minimise the loss
%function

    numberOfCells = numel( cells );
    
    targetEdges = edge( canvas );
    yStruct = struct();
    yStruct.targetEdges = targetEdges;
    yStruct.targetMask  = canvas ~= 1;
    yStruct.target      = canvas;
    yStruct.coverageMask = coverageMask;

    
    cellParamsArray = cell( numberOfCells, 1 );
    
    for j = 1:numberOfCells
        i = cells( j );
        cellParams.nucleusX = nucleusInfos{i}.centroid.x;
        cellParams.nucleusY = nucleusInfos{i}.centroid.y;

        %FIXME - for now we are assuming the cell center is the same as the nucleus
        cellParams.x = nucleusInfos{i}.centroid.x;
        cellParams.y = nucleusInfos{i}.centroid.y;
        
        cellParamsArray{ j } = cellParams;
    end
    
    if strcmp( optimisationParams.algorithm, 'simplex' )
        [ x, minLoss, canvas, mask ] = simplex( yStruct, cellParamsArray, allParams, optimisationParams );
    else
        error( [ 'Unknown method ', optimisationParams.algorithm ] );
    end
end

function [ x, minLoss, canvas, mask ] = simplex( yStruct, cellParamsArray, allParams, optimisationParams )
    f = @( x )( simplexLossFn( x, yStruct, cellParamsArray, allParams, optimisationParams ) );
    
    if optimisationParams.visuallyVerbose %'PlotFcns', @optimplotfval
        options = optimset( 'Display', 'iter', 'TolFun', optimisationParams.epsilon, 'TolX', optimisationParams.epsilon );
    else
        %options = optimset( 'DiffMinChange', 0.0, 'TolFun', optimisationParams.epsilon, 'TolX', optimisationParams.epsilon );
    end
    
    %options = optimoptions( 'fmincon', 'Display', 'Iter', 'DiffMinChange', 0.1 );
    options = optimoptions( 'fmincon', 'Diagnostics', 'on', 'DiffMinChange', 0.1 );
    numberOfCells = numel( cellParamsArray );
    x0 = repmat( optimisationParams.x0, 1, numberOfCells );
    lb = repmat( optimisationParams.lowerBound, 1, numberOfCells );
    ub = repmat( optimisationParams.upperBound, 1, numberOfCells );
    
    
    if( 0 )
        [ x, minLoss, exitFlag ] = fminsearch( f, x0, options );
        %[ x, minLoss, exitFlag ] = fminsearchbnd( f, x0, lb, ub, options );
        disp( [ 'Exit Flag Was ', num2str( exitFlag ) ] );
    else
        A = eye( numel( lb ), numel( lb ) );
        %A( 1:3, 1:3 ) = A( 1:3, 1:3 ) * -1;
        [ x, minLoss, exitFlag ] = fmincon( f, x0, A, ub );
    end
    pSpace        = struct();
    
    for i = 0:( numberOfCells - 1 )
        pSpace.radius                   = x( i * 3 + 1 );
        pSpace.majorVsMinor             = x( i * 3 + 2 );
        pSpace.majorVsMinorAngle        = x( i * 3 + 3 );
        cellParams                      = abcStructureUnion( pSpace, cellParamsArray{ i + 1 } );
        cellParamsArray{ i + 1 }        = cellParams;
    end

    [ canvas, mask ] = abcDrawCells( cellParamsArray, allParams );
end

function loss = simplexLossFn( x, yStruct, cellParamsArray, allParams, optimisationParams )
    if( 1 || optimisationParams.visuallyVerbose )
        disp( x );
    end
    
    numberOfCells = length( x ) / 3;
    %cellParamsArray = cell( numberOfCells, 1 );
    pSpace                          = struct();
    for i = 0:( numberOfCells - 1 )
        pSpace.radius                   = x( i * 3 + 1 );
        pSpace.majorVsMinor             = x( i * 3 + 2 );
        pSpace.majorVsMinorAngle        = x( i * 3 + 3 );
        cellParams                      = abcStructureUnion( pSpace, cellParamsArray{ i + 1 } );
        cellParamsArray{ i + 1 }        = cellParams;
    end
    
    [ canvas, mask ]  = abcDrawCells( cellParamsArray, allParams );
    xStruct.canvas    = canvas;
    xStruct.mask      = mask;
    loss              = lossFn2( xStruct, yStruct, optimisationParams );                       
    
end


function l = lossFn2( xStruct, yStruct, optimisationParams )
    canvas  = xStruct.canvas .* yStruct.coverageMask;
    targetCanvas = yStruct.target .* yStruct.coverageMask;
    
    tmp0 = ( targetCanvas - canvas );
    
    if( optimisationParams.visuallyVerbose )
        imshow( abs( tmp0 ) );
        pause
        %pause( 0.1 );
    end
    
    tmp1 = tmp0( yStruct.coverageMask ); 
    l    = sum( sum( tmp1 .* tmp1 ) );
    %disp( l );
    %pause;
end

function l = lossFn( xStruct, yStruct, optimisationParams )
%create a distance matrix for the original canvas
%then for the points on the edge of our new specific cell, sum up the total
%distance, making sure you dont pick ones that are too far away
%('truncated')

    isTruncatedSquareLoss = strcmp( optimisationParams.lossFunction, 'truncatedSquareLoss' );
    targetEdges = yStruct.targetEdges;
    cellCanvas  = xStruct.canvas;
            
    if optimisationParams.checkUniformIntensity
        canvas      = yStruct.target;
        intensities = canvas( xStruct.cellMask );
    
        %avgInt = mean( intensities );
        modInt = mode( intensities );
        Q90    = quantile( intensities, 0.95 );
        intensityScore =( Q90 - modInt );
    else
        intensityScore = 0;
    end
    
    cellEdge = edge( cellCanvas );

    distanceToEdge = bwdist( cellEdge, 'euclidean' );

    if( 0 && optimisationParams.visuallyVerbose )
        tmp1 = distanceToEdge .* targetEdges;
        imagesc( tmp1 );
        pause( 0.1 );
    end

    edgeToEdge = distanceToEdge( targetEdges );
    edgeToEdgeSorted = edgeToEdge(:);
    edgeToEdgeSorted = sort( edgeToEdgeSorted );

    nonZeros = sum( sum( cellEdge ) );
    nonZeros = min( nonZeros, length( edgeToEdgeSorted ) );    

    tmp0 = edgeToEdgeSorted( 1:nonZeros );    
    if( isTruncatedSquareLoss )
        tmp1 = tmp0( tmp0 < optimisationParams.truncation );
    else
        tmp1 = tmp0;
    end
    
    l = ( tmp1' * tmp1 ) / nonZeros;
    
    %disp( [ 'IntensityScore=', num2str( Q90 ), '---',  num2str( modInt ), ', DistanceScore = ', num2str( l ) ] );
    %disp( [ 'IntensityScore=', num2str( intensityScore ), ', DistanceScore = ', num2str( l ) ] );
    
    l = l + intensityScore;
end 