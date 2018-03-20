function [corners, nMatches, avgErr] = findBoard(I,Plotimg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  I - The chessboard image
%         Plotimg - Whether to show the intermediate steps
% Output: corners - The locations of the four outer corners as a 4x2 array, in
%                  the form [ [x1,y1]; [x2,y2]; ... ].
%         nMatches - Number of matching points found (ideally is 81)
%         avgErr - The average reprojection error of the matching points
%
% Find a 8x8 checkerboard in the image I.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the output
corners = [];
nMatches = [];
avgErr = [];

% Convert RGB image to grayscale
if size(I,3)>1
    I = rgb2gray(I);
end

% Record the processing time
tic

% Get the edge image
[~,thresh] = edge(I, 'canny');      % First get the automatic threshold
E = edge(I, 'canny',5*thresh);     % Raise the threshold

if Plotimg
    figure(1), imshow(E), title('Edges');
end

%% Implement the Hough transform to find lines.
[H,thetaValues,rhoValues] = hough(E);

% Extract peaks from the Hough array H.  Parameters for this:
%   houghThresh:  Minimum value to be considered a peak.  Default
%                 is 0.5*max(H(:))
%   NHoodSize:  Size of suppression neighborhood.  Default is
%               [size(H,1)/50, size(H,2)/50].  Must be odd numbers.
myThresh = ceil(0.05*max(H(:)));
NHoodSize = ceil([size(H,1)/70, size(H,2)/70]);
NumLines = 35;

% Force odd size
if mod(NHoodSize(1),2)==0  NHoodSize(1) = NHoodSize(1)+1;  end
if mod(NHoodSize(2),2)==0  NHoodSize(2) = NHoodSize(2)+1;  end
peaks = houghpeaks(H, ...
    NumLines, ...             % Maximum number of peaks to find
    'Threshold', myThresh, ...      % Threshold for peaks
    'NHoodSize', NHoodSize);    % Default = floor(size(H)/50);

% Display Hough array and draw peaks on Hough array.
if Plotimg
    figure(2), imshow(H, []), title('Hough Transform'), impixelinfo;
    for i=1:size(peaks,1)
        rectangle('Position', ...
            [peaks(i,2)-NHoodSize(2)/2, peaks(i,1)-NHoodSize(1)/2, ...
            NHoodSize(2), NHoodSize(1)], 'EdgeColor', 'r');
    end
end

%% Lines clustering in the scaled H-space using k-means

% Scale the H-space to weight theta and rho equally
minwidth = min(length(thetaValues),length(rhoValues));
Tkm = fitgeotrans([1 1; length(thetaValues) 1;length(thetaValues) length(rhoValues); 1 length(rhoValues)], [1 1; minwidth 1; minwidth minwidth; 1 minwidth], 'projective');

shrinkpeaks = transformPointsForward(Tkm, ...
    [peaks(:,2) peaks(:,1)]);

lines = [thetaValues(peaks(:,2));rhoValues(peaks(:,1))]';

% Shift the H-space if one set of lines is close to 90 degree
shift = 0;
if min(shrinkpeaks(:,1))<10 && sum(shrinkpeaks(:,1)<10)>=3
    
    % Sort the lines based on its angles
    [~,d2] = sort(shrinkpeaks(:,1),1);
    sortpeaks=shrinkpeaks(d2,:);
    
    % Find the angle needs to be shifted
    for i = 1: round(length(sortpeaks)/2)
        if (sortpeaks(i+1,2) - sortpeaks(i,2)) > 20
            shift = sortpeaks(i,1) + 1;
        else
            continue
        end
    end
end

Shiftshrinkpeaks = shrinkpeaks;
Shiftshrinkpeaks(:,1) = Shiftshrinkpeaks(:,1) + shift;
for k=1:length(Shiftshrinkpeaks)
    if Shiftshrinkpeaks(k,1)>180
        Shiftshrinkpeaks(k,1)=Shiftshrinkpeaks(k,1)-180;
    end
end
index = kmeans(Shiftshrinkpeaks,2);

% If bad initialization happens, do it again
while (sum(index==1)< min(round(NumLines/3),10) ) | (sum(index==2)< max(round(NumLines/3),10) )
    fprintf('Bad k-means initialization, we will do it again.\n');
    index = kmeans(Shiftshrinkpeaks,2);
end

lines1 = [];
lines2 = [];

% Store the lines in two matrices based on the clustering result
for k=1:size(index,1)
    if index(k)==1
        lines1 = [lines1 lines(k,:)'];
    else
        lines2 = [lines2 lines(k,:)'];
    end
end

lines1srhinkindex=[];
lines2srhinkindex=[];
for k=1:size(index,1)
    if index(k)==1
        lines1srhinkindex = [lines1srhinkindex ;shrinkpeaks(k,:)];
    else
        lines2srhinkindex = [lines2srhinkindex ;shrinkpeaks(k,:)];
    end
end

Hsrhink = imwarp(H, Tkm);

% Draw the k-means clustering result in the scaled H-space
if Plotimg
    figure(3), imshow(Hsrhink, []), title('K-means result in the scaled H-space'), impixelinfo;
    for i=1:size(lines1srhinkindex,1)
        rectangle('Position', ...
            [lines1srhinkindex(i,1)-8/2, lines1srhinkindex(i,2)-8/2, ...
            8, 8], 'EdgeColor', 'r');
    end
    for i=1:size(lines2srhinkindex,1)
        rectangle('Position', ...
            [lines2srhinkindex(i,1)-8/2, lines2srhinkindex(i,2)-8/2, ...
            8, 8], 'EdgeColor', 'g');
    end
end

%% Outliers elimination
% For each point in the H-space, we calculate the distance between itself
% and its nearest neighbour and assign to that point. If a point's assigned
% distance is way too large compared with other points in the same group,
% it will be marked as an outlier and eliminated.

distanceVec1 = zeros(1,size(lines1,2));

for k=1:size(lines1srhinkindex,1)
    nearpoints=knnsearch(lines1srhinkindex, lines1srhinkindex(k,:), 'k', 2);
    distanceVec1(k)=norm( lines1srhinkindex(nearpoints(1),:) - lines1srhinkindex(nearpoints(2),:) );
end
avgDis1 = mean(distanceVec1);
delete1 = find( distanceVec1 >= (avgDis1*3) );
lines1srhinkindex(delete1,:)=[];

distanceVec2 = zeros(1,size(lines2,2));
for k=1:size(lines2srhinkindex,1)
    nearpoints=knnsearch(lines2srhinkindex, lines2srhinkindex(k,:), 'k', 2);
    distanceVec2(k)=norm( lines2srhinkindex(nearpoints(1),:) - lines2srhinkindex(nearpoints(2),:) );
end
avgDis2 = mean(distanceVec2);
delete2 = find( distanceVec2 >= (avgDis2*3) );
lines2srhinkindex(delete2,:)=[];

% Plot in the scaled H-space with outlier elimination
if Plotimg
    figure(4), imshow(Hsrhink, []), title('K-means and outlier elimination result'), impixelinfo;
    for i=1:size(lines1srhinkindex,1)
        rectangle('Position', ...
            [lines1srhinkindex(i,1)-8/2, lines1srhinkindex(i,2)-8/2, ...
            8, 8], 'EdgeColor', 'r');
    end
    for i=1:size(lines2srhinkindex,1)
        rectangle('Position', ...
            [lines2srhinkindex(i,1)-8/2, lines2srhinkindex(i,2)-8/2, ...
            8, 8], 'EdgeColor', 'g');
    end
end

%% Sort the lines and find their corresponding intersections
lines1(:,delete1) = [];
lines2(:,delete2) = [];

lines1 = sortLines(lines1, size(E));
lines2 = sortLines(lines2, size(E));

% Sort the lines, from top to bottom (for horizontal lines) and left to
% right (for vertical lines).

fprintf('Most common angles range: [%.1f %.1f] and [%.1f %.1f] in degree\n', min(lines1(1,:)), max(lines1(1,:)),min(lines2(1,:)), max(lines2(1,:)));

% Intersect every pair of lines, one from set 1 and one from set 2.
% Output is the x,y coordinates of the intersections:
%   xIntersections(i1,i2): x coord of intersection of i1 and i2
%   yIntersections(i1,i2): y coord of intersection of i1 and i2
[xIntersections, yIntersections] = findIntersections(lines1, lines2);

% Plot all measured intersection points.
if Plotimg
    
    figure(5), imshow(E), title('Sorted lines and their intersections');
    drawLines( ...
        lines1(2,:), ...    % rhos for the lines
        lines1(1,:), ...    % thetas for the lines
        size(E), ...        % size of image being displayed
        'r');               % color of line to display
    drawLines( ...
        lines2(2,:), ...    % rhos for the lines
        lines2(1,:), ...    % thetas for the lines
        size(E), ...        % size of image being displayed
        'g');               % color of line to display
    
    hold on
    plot(xIntersections(:),yIntersections(:),'yd');
    hold off
end

%% Define a "reference" image.

% Reference image is (IMG_SIZE_REF+1) x (IMG_SIZE_REF+1)
IMG_SIZE_REF = 319;

% Get predicted intersections of lines in the reference image.
[xIntersectionsRef, yIntersectionsRef] = createReference(IMG_SIZE_REF,Plotimg);

%% Geometric projection and reference matching
% Find the best correspondence between the points in the input image and
% the points in the reference image.  If found, the output is the four
% outer corner points from the image, represented as a a 4x2 array, in the
% form [ [x1,y1]; [x2,y2]; ... ].

[corners, nMatches, avgErr] = findCorrespondence( ...
    xIntersections, yIntersections, ...         % Input image points
    xIntersectionsRef, yIntersectionsRef, ...   % Reference image points
    I,Plotimg);

% Print the processing time
fprintf('Finding the board takes %f seconds\n',toc)
end
