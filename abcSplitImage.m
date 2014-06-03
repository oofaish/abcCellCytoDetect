function [ splits, masks ] = abcSplitImage( nucleusInfos, optimisationParams, canvasSize )

    %now start from the centre of each cell and go as far as maxRadius.
    %count how many cells are in that area?

    splits = zeros( numel( nucleusInfos ), numel( nucleusInfos ) );

    hCounter = 1:canvasSize( 1 );
    wCounter =  1:canvasSize( 2 );
    
    maxAllowedDist = optimisationParams.maxRadius;

    cellMasks = cell( numel( nucleusInfos ), 1 );
    
    for i = 1:numel( nucleusInfos )
        nucleusInfo = nucleusInfos{i};
        x1 = nucleusInfo.centroid.x;
        y1 = nucleusInfo.centroid.y;
        mask   = bsxfun( @plus, ( hCounter - x1 ).^2, ( wCounter' - y1 ).^2 ) < maxAllowedDist^2;
        cellMasks{ i } = mask;
        
        [ row, column ] = find( splits == i );
        if( isempty( row ) )
            iToUse = i;
            totalInThisSplit = 1;
            splits( i, totalInThisSplit ) = i;
        else
            iToUse = row;
            totalInThisSplit = sum( splits( row, :) ~= 0 );
        end

        for j = i + 1:numel( nucleusInfos )
            nucleusInfo = nucleusInfos{j};
            x2 = nucleusInfo.centroid.x;
            y2 = nucleusInfo.centroid.y;

            dist = ( ( x1 - x2 ) ^ 2 + ( y1  - y2 ) ^ 2 ) ^ 0.5;
            
            if dist < maxAllowedDist && isempty( find( splits( iToUse, : ) == j, 1 ) )
                totalInThisSplit = totalInThisSplit + 1;
                splits( iToUse, totalInThisSplit ) = j;
            end
        end
    end

    splits = splits(max(splits')>0,:);
    masks = cell( size( splits, 1 ), 1 );
    for i = 1:size( splits, 1 )
        for j = 1:size( splits, 2 )
            if splits( i, j ) > 0
                mask = cellMasks{ splits( i, j ) };
                mask( masks{ i } == 1 ) = 1;
                masks{ i } = mask;
            end
        end
    end
    
end



