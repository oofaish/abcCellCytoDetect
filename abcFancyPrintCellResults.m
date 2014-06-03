function  abcFancyPrintCellResults( i, score, x, actualX  )

    fieldNames = { 'radius', 'axis', 'angle' };
    
    printLine = [ 'Cell: ', num2str( i ), ', Score=', num2str( score ) ];
    
    
    %for i = 1:numel( fieldNames )
    %    fieldName = fieldNames( i );
    %    fieldName = fieldName{1};
    %    printLine = [ printLine, [ ', ', fieldName{i}, ' diff=', num2str( abs( pSpace.(fieldName ) - actualPSpace.(fieldName ) ) ) ] ];
    %end
    
    disp( [ num2str( i ), '-', num2str( score ) ] );
    disp( x - actualX );
    
    %disp( printLine );

%disp( [ 'Cell: ', num2str( i ), 'score= ', num2str( score ), ', Actual= ' ] );
%    disp( cellInfos{ i }.pSpace );
%    disp( 'Found=' );
%    disp( pSpace );
