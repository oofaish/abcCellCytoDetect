function l = LossFnSimple( xStruct, yStruct, optimisationParams )
    canvas       = double( xStruct.canvas );% .* yStruct.coverageMask;
    targetCanvas = double( yStruct.canvas );% .* yStruct.coverageMask;
    tmp0 = ( targetCanvas - canvas ) / 255.0;
    
    tmp1 = tmp0(:); 
    %tmp1 = tmp0( yStruct.coverageMask );
    
    %if 1 || strcmp( optimisationParams.lossMethod, 'vector' ) 
    %    l = tmp1;
    %else
    l    = sum( sum( tmp1 .* tmp1 ) );
    %end
    %pause;
end