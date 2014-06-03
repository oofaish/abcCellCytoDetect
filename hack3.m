f = @( x ) (  ( x(1) - 3 ) .^ 2  + ( x(2) - 5 ) .^ 2 );

A = eye( 2, 2);
b = [ 1000, 1000 ];
x0 = [ 15, 7 ];
%options = optimoptions( 'fmincon', 'Display', 'Iter', 'DiffMinChange', 0.1 );
[ x, fval, flag ] = fmincon( f, x0, A, b );

x