x5 = 1:1e7;
x6 = 1:1e4;
counter2 = arrayfun( @(x) ( builtin('_ismemberoneoutput',x, x5 ) ), x6 );
counter3 = binarySearch( x5, x6,[], 0 );
