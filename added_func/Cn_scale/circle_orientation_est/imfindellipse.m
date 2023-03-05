function [ellipsestack,centers,propsstack] = imfindellipse( projections,thr,minWidth,showfig )
if ~exist('showfig','var')
    showfig=0;
end
if ~exist('minWidth','var')
    minWidth=2;
end
[rows, ~, nImages] = size(projections);
ellipsestack = zeros(size(projections));
centers = zeros(nImages,2);
for i=1:nImages
    Image = projections(:,:,i);
% Binarize the image
% Get the mask where the region is solid.
thr_i = thr * max(Image(:));
binaryImage = Image > thr_i;
% Fill it and take the largest blob:
binaryImage = imfill(binaryImage, 'holes');
binaryImage = bwmorph(binaryImage, 'bridge');
binaryImage = bwareafilt(binaryImage, 1);

% Get a new binary image where it only includes "wide" rows.
regionWidth = zeros(rows, 1);
for row = 1 : rows
  thisWidth = find(binaryImage(row, :), 1, 'last') - find(binaryImage(row, :), 1, 'first');
  if ~isempty(thisWidth)
    regionWidth(row) = thisWidth;
  end
end
% Define what "too narrow" is
binaryImage(regionWidth < minWidth, :) = false; % Erase too narrow rows.

% Measure elliptical properties
props = regionprops(binaryImage,'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Centroid', 'Perimeter','Area');
[~,m] = max([props.Area]);
props = props(m);
%------------------------------------------------------------------
% Create an ellipse with specified
% semi-major and semi-minor axes, center, and image size.
% From the FAQ: https://matlab.wikia.com/wiki/FAQ#How_do_I_create_an_ellipse.3F
xCenter = props.Centroid(1);
yCenter = props.Centroid(2);
xRadius = props.MinorAxisLength / 2;
yRadius = props.MajorAxisLength / 2;
% Make an angle array of about the same number of angles as there are pixels in the perimeter.
theta = linspace(0, 2*pi, ceil(props.Perimeter));
x = xRadius * cos(theta);
y = yRadius * sin(theta);
% Now we might need to rotate the coordinates slightly.  Make a rotation matrix
% Reference: https://en.wikipedia.org/wiki/Rotation_matrix
angleInDegrees = props.Orientation - 90;
rotationMatrix = [cosd(angleInDegrees), -sind(angleInDegrees); sind(angleInDegrees), cosd(angleInDegrees)];
% Now do the rotation
xy = [x', y'];
xyRotated = xy * rotationMatrix;
x = round(xyRotated(:, 1)+xCenter);
y = round(xyRotated(:, 2)+yCenter);
ellipse = zeros(size(binaryImage));
linearidx = sub2ind(size(ellipse),y,x);
ellipse(linearidx) = 1; 
ellipsestack(:,:,i) = ellipse;
centers(i,:) = [xCenter,yCenter];
propsstack(i) = props;
%plot(xCenter, yCenter, 'r+', 'LineWidth', 2, 'MarkerSize', 13); % Put a cross at teh ellipse center.
end

if showfig
figure;viewellipse(projections,1,4,0,centers,ellipsestack);
fprintf('X Center = column #%.1f\n', xCenter);
fprintf('Y Center = row #%.1f\n', yCenter);
fprintf('Minor Axis Length = %.1f\n', xRadius);
fprintf('Major Axis Length = %.1f\n', yRadius);
fprintf('Orientation = %.1f degrees\n', angleInDegrees);
end
end

