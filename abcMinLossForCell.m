function [ pSpace, minLoss, cellMask ] = abcMinLossForCell( pSpace, cellParams, canvas, nucleusInfo, cellInfo, optimisationParams )
%ABCMINLOSSFORCELL change parameters in pSpace to minimise the loss
%function


    cellParams.nucleusX = nucleusInfo.centroid.x;
    cellParams.nucleusY = nucleusInfo.centroid.y;

    %FIXME - for now we are assuming the cell center is the same as the nucleus
    cellParams.x = nucleusInfo.centroid.x;
    cellParams.y = nucleusInfo.centroid.y;

    targetEdges = edge( canvas );
    yStruct = struct();
    yStruct.targetEdges = targetEdges;
    
    %this is a bit naughty - I am mixing cellParams and allParams
    %will probably blow in my face later.
    allParams = abcParams( cellParams );
    %ceParams = abcStructureUnion( cellParams, abcParams() );
    if strcmp( optimisationParams.algorithm, 'gridsearch' )
        [ pSpace, minLoss, cellMask ] = gridSearch( pSpace, yStruct, allParams, optimisationParams );
    elseif strcmp( optimisationParams.algorithm, 'simplex' )
        [ pSpace, minLoss, cellMask ] = simplex( pSpace, yStruct, allParams, optimisationParams );
    end
end

function [ pSpace, minLoss, cellMask ] = gridSearch( pSpace, yStruct, allParams, optimisationParams )
    rs = optimisationParams.minRadius:optimisationParams.maxRadius;%ALI - just for funs
    losses = zeros( size( rs ) );
    cellMasks = cell( size( losses ) );
    xStruct = struct();
    i = 0;
    for radius = rs
        i = i + 1;
        %disp( radius );
        pSpace.radius = radius;
        cellParams = abcStructureUnion( pSpace, allParams );
        [ cellCanvas, cellMask, ~, ~ ] = abcDrawCell( cellParams, allParams.canvasSize );
        xStruct.cellCanvas = cellCanvas;
        losses( i ) = lossFn( xStruct, yStruct, optimisationParams );                       
        cellMasks{ i } = cellMask;
    end
    
    [ minLoss, minIndex ] = min( losses );
    
    if( optimisationParams.visuallyVerbose )
        h = figure;
        plot( rs, losses );
        ms2 = linspace( minLoss * 0.9, max( losses ) * 1.1, 10 );
        hold on;
        if( isfield( optimisationParams, 'actualRadius' ) )
            hold on;
            plot( zeros( 10, 1) + optimisationParams.actualRadius, ms2,'r' );
            hold off
        end
        pause;
        try
            close( h );
        catch
            %Dont bother doing anything - this might happen is window is
            %closed.
        end
        
    end
    
    minRadius = rs( minIndex );
    cellMask = cellMasks{ minIndex };
    pSpace.radius = minRadius;

end

function [ pSpace, minLoss, cellMask ] = simplex(  pSpace, yStruct, allParams, optimisationParams )
    
    x = [ pSpace.radius ];
    f = @( x )( simplexLossFn( x, yStruct, allParams, optimisationParams ) );
    
    if optimisationParams.visuallyVerbose
        options = optimset( 'Display', 'iter', 'PlotFcns', @optimplotfval );
    else
        options = optimset( );
    end

    [ x, minLoss ] = fminsearch( f, optimisationParams.x0, options );
    
    pSpace.radius = x( 1 );
    
    %get the mask
    cellParams  = abcStructureUnion( pSpace, allParams );
    [ ~, cellMask, ~, ~ ] = abcDrawCell( cellParams, allParams.canvasSize );
end

function loss = simplexLossFn( x, yStruct, allParams, optimisationParams )
    if x(1) > optimisationParams.minRadius
        pSpace                          = struct();
        pSpace.radius                   = x( 1 );
        cellParams                      = abcStructureUnion( pSpace, allParams );
        [ cellCanvas, ~, ~, ~ ]         = abcDrawCell( cellParams, allParams.canvasSize );
        xStruct.cellCanvas              = cellCanvas;
        loss                           = lossFn( xStruct, yStruct, optimisationParams );                       
    else
        loss = 10000;
    end
   
end

function l = lossFn( xStruct, yStruct, optimisationParams )
%create a distance matrix for the original canvas
%then for the points on the edge of our new specific cell, sum up the total
%distance, making sure you dont pick ones that are too far away
%('truncated')

    isTruncatedSquareLoss = strcmp( optimisationParams.lossFunction, 'truncatedSquareLoss' );
    targetEdges = yStruct.targetEdges;
    cellCanvas  = xStruct.cellCanvas;

    cellEdge = edge( cellCanvas );

    distanceToEdge = bwdist( cellEdge, 'euclidean' );

    if( optimisationParams.visuallyVerbos )
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
    
    l = tmp1' * tmp1;
end 