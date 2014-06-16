function [ x, minLoss ] = abcMinLossHack( yStruct, allParams, optimisationParams )
%ABCMINLOSSHACK change parameters in pSpace to minimise the loss
%function    
    if strcmp( optimisationParams.method, 'gridSearch' )
        [ x, minLoss] = gridSearch( yStruct, allParams, optimisationParams );
    elseif strcmp( optimisationParams.method, 'lm' )
        [ x, minLoss] = lm( yStruct, allParams, optimisationParams );
    else
        error( 'unknown method' );
    end
end

function [ x, minLoss ] = gridSearch( yStruct, allParams, op )

    f = @( x )( topLossFn( x, yStruct, allParams, op ) );

    xSpace1 =  op.xMin( 1 ):op.xStep( 1 ):op.xMax( 1 );
    xSpace2 =  op.xMin( 2 ):op.xStep( 2 ):op.xMax( 2 );
    
    xGrid = allcomb( xSpace1, xSpace2 );
    minLoss = f( xGrid(1,:) );
    x = xGrid( 1,:);
    for row=2:size( xGrid, 1 )
        newLoss = f( xGrid( row, :) );
        if newLoss < minLoss
            x = xGrid( row, : );
            minLoss = newLoss;
            disp( 'New One' );
        end
        disp( xGrid( row, : ) );
        disp( num2str( newLoss ) );
    end
    
    x = [ 0 0 x( 1 ) x( 2 ) ];
end

function [ x, minLoss ] = lm( yStruct, allParams, optimisationParams )
    f = @( x )( topLossFn( x, yStruct, allParams, optimisationParams ) );
    
    %x0 = optimisationParams.x0;
    %gradStep = repmat( [ 2, 2, 0.1, 0.1 ], 1, 1 );
    gradStep = optimisationParams.FinDiffRelStep;
    x0 = optimisationParams.x0;
    %FIXME - remember you are blurring.
    %DiffMaxChange = 0.5;%'Algorithm', 'sqp'    
    %DiffMaxChange = 10;%'Algorithm', 'sqp'    
        %'FinDiffRelStep', gradStep,...
        %'DiffMaxChange', DiffMaxChange, ...
    
    options = optimoptions( 'lsqnonlin', ...
        'Algorithm', 'levenberg-marquardt', ...
        'FinDiffRelStep', gradStep,...
        'TolX', 1e-5, ...
        'TolFun', 1e-4,...
        'Diagnostics', 'on',...
        'Display', 'Iter', ...
        'FinDiffType', 'central' );%, 'PlotFcns', @optimplotfval );'FinDiffType', 'central'
        lb = [];
        ub = [];
        
    if( isfield( optimisationParams, 'DiffMinChange' ) )
        options = optimoptions( options, 'DiffMinChange', optimisationParams.DiffMinChange );
    end
    
    if( isfield( optimisationParams, 'DiffMaxChange' ) )
        options = optimoptions( options, 'DiffMaxChange', optimisationParams.DiffMaxChange );
    end

    
    [x,minLoss,residualVector,exitFlag ] = lsqnonlin( f, x0, lb, ub, options );
    x = [ 0 0 x( 1 ) x( 2 ) ];
end

function loss = topLossFn( x, yStruct, allParams, optimisationParams ) 
    %x( 1 ) = round( x( 1 ) );
    %x( 2 ) = round( x( 2 ) );
    x = [ 0 0 x( 1 ) x( 2 ) ];
    result = abcTransformCell( optimisationParams.c, optimisationParams.m, optimisationParams.n, x, allParams, yStruct.nCentroid );
    xStruct.canvas    = result.c;
    xStruct.mask      = result.m;
    loss              = LossFn( xStruct, yStruct, optimisationParams );                       
end

function l = LossFn( xStruct, yStruct, optimisationParams )
    %canvas  = double( edge( xStruct.canvas, 'canny' ) );% .* yStruct.coverageMask;
    %targetCanvas = double( edge( yStruct.c, 'canny' ) );% .* yStruct.coverageMask;
    %tmp0 = targetCanvas - canvas;
    
    
    canvas       = double( xStruct.canvas ) .* yStruct.coverageMask;
    targetCanvas = double( yStruct.c ) .* yStruct.coverageMask;
    tmp0 = ( targetCanvas - canvas ) / 255.0;
    
    if( 1 || optimisationParams.visuallyVerbose )
        %subplot( 133 );
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
    
    %tmp1 = tmp0(:); 
    tmp1 = tmp0( yStruct.coverageMask );
    
    if strcmp( optimisationParams.lossMethod, 'vector' ) 
        l = tmp1;
    else
        l    = sum( sum( tmp1 .* tmp1 ) );
    end
    %pause;
end