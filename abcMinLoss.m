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
    
    %if strcmp( optimisationParams.algorithm, 'simplex' )
    [ x, minLoss, canvas, mask ] = simplex( yStruct, cellParamsArray, allParams, optimisationParams );
    %else
    %    error( [ 'Unknown method ', optimisationParams.algorithm ] );
    %end
end

function [ x, minLoss, canvas, mask ] = simplex( yStruct, cellParamsArray, allParams, optimisationParams )
    f = @( x )( simplexLossFn( x, yStruct, cellParamsArray, allParams, optimisationParams ) );
    
    if optimisationParams.visuallyVerbose %'PlotFcns', @optimplotfval
        options = optimset( 'Display', 'iter', 'TolFun', optimisationParams.epsilon, 'TolX', optimisationParams.epsilon );
    else
        %options = optimset( 'DiffMinChange', 0.0, 'TolFun', optimisationParams.epsilon, 'TolX', optimisationParams.epsilon );
    end
    
    numberOfCells = numel( cellParamsArray );
    x0 = repmat( optimisationParams.x0, 1, numberOfCells );
    lb = repmat( optimisationParams.lowerBound, 1, numberOfCells );
    ub = repmat( optimisationParams.upperBound, 1, numberOfCells );
    gradStep = repmat( [ 1, 0.05, 0.05 ], 1, numberOfCells );
    %FIXME - remember you are blurring.
    %DiffMaxChange = 0.5;%'Algorithm', 'sqp'    
    DiffMaxChange = 1;%'Algorithm', 'sqp'    
    if( 0 )
        [ x, minLoss, exitFlag ] = fminsearch( f, x0, options );
        %[ x, minLoss, exitFlag ] = fminsearchbnd( f, x0, lb, ub, options );
        disp( [ 'Exit Flag Was ', num2str( exitFlag ) ] );
    elseif( strcmp( optimisationParams.lossMethod, 'vector' )  )
        options = optimoptions( 'lsqnonlin', ...
            'Algorithm', 'levenberg-marquardt', ...
            'TolX', 1e-5, ...
            'TolFun', 1e-4,...
            'Diagnostics', 'on',...
            'Display', 'Iter', ...
            'FinDiffRelStep', gradStep,...
            'DiffMaxChange', DiffMaxChange, ...
            'FinDiffType', 'central' );%, 'PlotFcns', @optimplotfval );'FinDiffType', 'central'
            lb = [];
            ub = [];
        [x,minLoss,residualVector,exitFlag ] = lsqnonlin( f, x0, lb, ub, options );
    elseif( strcmp( optimisationParams.lossMethod, 'double' ) )
        options = optimoptions( 'fmincon','Algorithm', 'sqp', 'TolX', 1e-4, 'TolFun', 1e-4, 'Diagnostics', 'on', 'Display', 'Iter', 'FinDiffRelStep', gradStep, 'DiffMaxChange', DiffMaxChange,'FinDiffType', 'central' );%, 'PlotFcns', @optimplotfval );'FinDiffType', 'central'
        A = eye( numel( lb ), numel( lb ) );
        %A( 1:3, 1:3 ) = A( 1:3, 1:3 ) * -1;
        [ x, minLoss, exitFlag, output, lambda, grad ] = fmincon( f, x0, [], [], [], [], lb, ub, [], options );
    else
        error( 'Unkown lossmethod' );
    end
    pSpace        = struct();
    
    if( 1 )
        for i = 0:( numberOfCells - 1 )
            pSpace.radius                   = x( i * 3 + 1 );
            pSpace.majorVsMinor             = x( i * 3 + 2 );
            pSpace.majorVsMinorAngle        = x( i * 3 + 3 );
            cellParams                      = abcStructureUnion( pSpace, cellParamsArray{ i + 1 } );
            cellParamsArray{ i + 1 }        = cellParams;
        end
    else %FIXME - 1 parameter hack
        for i = 1:numberOfCells
            pSpace.radius                   = x( i);
            cellParams                      = abcStructureUnion( pSpace, cellParamsArray{ i } );
            cellParamsArray{ i }            = cellParams;
        end        
    end

    [ canvas, mask ] = abcDrawCells( cellParamsArray, allParams );
end

function loss = simplexLossFn( x, yStruct, cellParamsArray, allParams, optimisationParams )    
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
    
    if( 1 || optimisationParams.visuallyVerbose )
        subplot( 133 );
        %colormap( 'default' );
        %axis image;
        imshow( ( tmp0 + 1 ) ./ 2 );
        
        pause( 0.001 );
        global writerObj;
        if( writerObj ~= -1 )
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        

    end
    
    tmp1 = tmp0( yStruct.coverageMask ); 
    
    if strcmp( optimisationParams.lossMethod, 'vector' ) 
        l = tmp1;
    else
        l    = sum( sum( tmp1 .* tmp1 ) );
    end
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
    
    l = double( sum( l + intensityScore ) );
end 