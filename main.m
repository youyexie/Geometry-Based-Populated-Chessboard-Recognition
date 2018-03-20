%% Load the board
clear variables
close all
clc

% The location of the chessboard image
dirNameBoard = strcat(pwd,'\Boards\');
Board = '1.jpg';

% Show the intermediate recognition process
Plotimg = true; 

%% Implement the chessboard recognition
I = imread(sprintf('%s/%s', dirNameBoard, Board));
Itemp = I;
if size(I,3) > 1
    I = rgb2gray(I);
end

% Find the populated chessboard
[corners, nMatches, avgErr] = findBoard(I,Plotimg);

%% Show the recognition result
% Get the pieces intersections
nx=9; % Intersection points in x direction
ny=9; % Intersection points in y direction
Intersections = getIntersection(corners,nx,ny,0);
Squarecenters = findSquarecenters(Intersections);

% Draw the boundary of the chessboard
figure(8),imshow(Itemp),title('Board detection result')
line([corners(1,1) corners(2,1)],[corners(1,2) corners(2,2)],  'Color', 'g','LineWidth',2);
line([corners(2,1) corners(3,1)],[corners(2,2) corners(3,2)],  'Color', 'g','LineWidth',2);
line([corners(3,1) corners(4,1)],[corners(3,2) corners(4,2)],  'Color', 'g','LineWidth',2);
line([corners(1,1) corners(4,1)],[corners(1,2) corners(4,2)],  'Color', 'g','LineWidth',2);

% Draw the intersections
minD = min(size(I));
DMIN = minD/50;
DMIN=DMIN;
figure(8)
for i=1:size(Intersections(1,:)',1)
    rectangle('Position', [Intersections(1,i)-DMIN/2 Intersections(2,i)-DMIN/2 DMIN DMIN], 'EdgeColor', 'g','Curvature',[1 1],'LineWidth',2);
end