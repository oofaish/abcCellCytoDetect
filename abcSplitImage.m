function [ splits, masks ] = abcSplitImage( nucleusInfos, optimisationParams, canvasSize )

    %now start from the centre of each cell and go as far as maxRadius.
    %count how many cells are in that area?

    splits = zeros( numel( nucleusInfos ), numel( nucleusInfos ) );

    hCounter = 1:canvasSize( 1 );
    wCounter =  1:canvasSize( 2 );
    
    maxAllowedDist = optimisationParams.maxRadius * 1;

    cellMasks = cell( numel( nucleusInfos ), 1 );
    
    for i = 1:numel( nucleusInfos )
        nucleusInfo = nucleusInfos{i};
        x1 = nucleusInfo.centroid.x;
        y1 = nucleusInfo.centroid.y;
        mask   = bsxfun( @plus, ( hCounter - x1 ).^2, ( wCounter' - y1 ).^2 ) < maxAllowedDist^2;
        cellMasks{ i } = mask;
        
        %[ row1, column1 ] = find( splits == i );

        %if( isempty( row ) )
        %    iToUse = i;
        %    totalInThisSplit = 1;
        %    splits( i, totalInThisSplit ) = i;
        %else
        %    iToUse = row;
        %    totalInThisSplit = sum( splits( row, :) ~= 0 );
        %end

        for j = 1:numel( nucleusInfos )
            if j == i 
                if i == 1
                    splits( 1, 1 ) = 1;
                end
                continue
            end
            
            nucleusInfo = nucleusInfos{j};
            x2 = nucleusInfo.centroid.x;
            y2 = nucleusInfo.centroid.y;
            dist = ( ( x1 - x2 ) ^ 2 + ( y1  - y2 ) ^ 2 ) ^ 0.5;
            
            if( dist < maxAllowedDist )
                [ row1, column1 ] = find( splits == i );
                [ row2, column2 ] = find( splits == j );
            
                if ~isempty( row2 ) 
                    if isempty( find( splits( row2, : ) == i, 1 ) )
                        totalInThisSplit = sum( splits( row2, :) ~= 0 );
                        splits( row2, totalInThisSplit + 1 ) = i;
                    end
                elseif ~isempty( row1 )
                    if isempty( find( splits( row1, : ) == j, 1 ) )
                        totalInThisSplit = sum( splits( row1, :) ~= 0 );
                        splits( row1, totalInThisSplit + 1 ) = j;
                    end
                else
                    row = sum( splits( :, 1) ~= 0 );%find first empty row
                    splits( row + 1, 1 ) = i;
                    splits( row + 1, 2 ) = j;
                end
            end
                
        end
        
        %if at the end of all of this, i is not anywhere, stick it
        %somewhere
        if isempty( find( splits == i ) )
            row = sum( splits( :, 1) ~= 0 );%find first empty row
            splits( row + 1, 1 ) = i;
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



