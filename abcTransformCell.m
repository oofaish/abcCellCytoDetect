function [ cellImage2, cellMask2, nucleusMask2 ] = abcTransformCell( cellImage, cellMask, nucleusMask, x, params )

canvasSize = params.canvasSize( 1 );
assert( params.canvasSize( 1 ) == params.canvasSize( 1 ), 'I am assuming height and width are the same' );

dx = x( 1 );
dy = x( 2 );

theta = x( 3 );
delta = x( 4 );

tx = canvasSize / 2 + dx;
ty = canvasSize / 2 + dy;
%theta = pSpace.theta;
%delta = pSpace.delta;

[Y_corners,M_homoTrans] = homoSpaceTransform(tx, ty, theta, delta, cellImage );

%+----------------------------------------+
%| Draw the cell mask on 512 x 512 canvas |
%+----------------------------------------+
Canvas = zeros( canvasSize, canvasSize );
Canvas = roipoly(Canvas, Y_corners(:,1), Y_corners(:,2));
[y_canvas,x_canvas]=find(Canvas==1);

Y_tilde = round(inv(M_homoTrans)*[x_canvas';y_canvas';ones(1,size(x_canvas,1))]);     
idxCellMask = find(cellMask == 1);
idxYTilde = (sub2ind(size(cellMask), round(Y_tilde(1,:)), round(Y_tilde(2,:))))';

%+-----------------------------------+
%| Draw the cell on 512 x 512 canvas |
%+-----------------------------------+
tmp_Img = uint8(zeros(canvasSize,canvasSize));
%tmp_Img2 = uint8(zeros(canvasSize,canvasSize));
tmp_pickup_Cell_Mask = false(canvasSize,canvasSize);
tmp_pickup_Cell_NucleiMask = false(canvasSize,canvasSize);

%I could not get this mother fucker to work
%counter = (1:size( y_canvas, 1 ) );

%counter2 = arrayfun( @(x) ( ismember( x, idxCellMask ) ), idxYTilde );
counter2 = arrayfun( @(x) ( builtin('_ismemberoneoutput',x,idxCellMask ) ), idxYTilde );

%counter3 = counter2;
%counter4 = counter( counter2 ~= 0 );
%counter51 = Y_tilde(1, counter4 );
%counter52 = Y_tilde(2, counter4 );
%counter6 = x_canvas( counter4, 1 );
%counter7 = y_canvas( counter4, 1 );
%tmp_Img2( counter6, counter7 ) = cellImage( counter51, counter52 );

for k = 1:size(y_canvas,1)
    idx = idxYTilde(k,1);
    isCellPixel = counter2( k );
    if isCellPixel == 1
        tmp_Img(x_canvas(k,1),y_canvas(k,1)) = cellImage(Y_tilde(1,k), Y_tilde(2,k));
        tmp_pickup_Cell_Mask(x_canvas(k,1),y_canvas(k,1)) = logical(cellMask(Y_tilde(1,k), Y_tilde(2,k)));
        tmp_pickup_Cell_NucleiMask(x_canvas(k,1),y_canvas(k,1)) = logical(nucleusMask(Y_tilde(1,k), Y_tilde(2,k)));
    end
end

%+----------------------------------+
%| Translate tmp_Img* to the center |
%|           if j == 1              |
%+----------------------------------+
img_1stCell_stat = regionprops(tmp_pickup_Cell_Mask, 'Centroid');
dx = tx - img_1stCell_stat(1,1).Centroid(1,2);
dy = ty - img_1stCell_stat(1,1).Centroid(1,1);

[xTmpImg, yTmpImg] = find(tmp_pickup_Cell_Mask == 1);
[xTmpNucleiMask, yTmpNucleiMask] = find(tmp_pickup_Cell_NucleiMask == 1);

cellImage2   = uint8(zeros(canvasSize,canvasSize));
cellMask2    = false(canvasSize,canvasSize);
nucleusMask2 = false(canvasSize,canvasSize);

%+--------------------------+
%| Move image and cyto mask |
%+--------------------------+
for k = 1:size(xTmpImg,1)
    cellImage2(round(xTmpImg(k,1) + dx), round(yTmpImg(k,1) + dy)) = ...
        tmp_Img(round(xTmpImg(k,1)), round(yTmpImg(k,1)));
    cellMask2(round(xTmpImg(k,1) + dx), round(yTmpImg(k,1) + dy)) = 1;
end

%+------------------+
%| Move Nuclei Mask |
%+------------------+
for k = 1:size(xTmpNucleiMask,1)
    nucleusMask2(round(xTmpNucleiMask(k,1) + dx), round(yTmpNucleiMask(k,1) + dy)) = 1;
end



end

