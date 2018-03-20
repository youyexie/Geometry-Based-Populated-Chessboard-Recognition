function [corners, nMatchesBest, avgErrBest] = findCorrespondence( ...
    xIntersections, yIntersections, ...         % Input image points
    xIntersectionsRef, yIntersectionsRef, ...   % Reference image points
    E,Plotimg)                                     % Edge image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  xIntersections,yIntersections - The coordinates of intersections
%         xIntersectionsRef, yIntersectionsRef - The reference coordinates
%                                                of intersections
%         E - The edge image
%         Plotimg - Whether to show the intermediate steps
% Output: corners - The locations of the four outer corners as a 4x2 array, in
%                  the form [ [x1,y1]; [x2,y2]; ... ].
%         nMatchesBest - Number of matching points found (ideally is 81)
%         avgErrBest - The average reprojection error of the matching points
%
% Find the best correspondence between the points in the input image and
% the points in the reference image.  If found, the output is the four
% outer corner points from the image, represented as a a 4x2 array, in the
% form [ [x1,y1]; [x2,y2], ... ].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the coordinates of the four outer corners of the reference image,
% in clockwise order starting from the top left.
pCornersRef = [ ...
    xIntersectionsRef(1,1), yIntersectionsRef(1,1);
    xIntersectionsRef(1,end), yIntersectionsRef(1,end);
    xIntersectionsRef(end,end), yIntersectionsRef(end,end);
    xIntersectionsRef(end,1), yIntersectionsRef(end,1) ];

M = 6;      % Number of lines to search in each direction
DMIN = 10;   % radius for the intersection detection

nMatchesBest = 0;   % Number of matches of best candidate found so far
avgErrBest = 1e10;   % The average error of the best candidate

N1 = size(xIntersections,1);
N2 = size(xIntersections,2);

M1 = round(N1*0.6);
M2 = round(N2*0.6);

for i1a=1:min(M1,N1)
    for i1b=N1:-1:max(N1-M1,i1a+8)
        for i2a=1:min(M2,N2)
            for i2b=N2:-1:max(N2-M2,i2a+8)
                
                if (i1a == i1b) || (i2a == i2b)
                    continue
                end
                
                % Get the four corners corresponding to the intersections
                % of lines (1a,2a), (1a,2b), (1b,2b, and (1b,2a).
                pCornersImg = zeros(4,2);
                pCornersImg(1,:) = [xIntersections(i1a,i2a) yIntersections(i1a,i2a)];
                pCornersImg(2,:) = [xIntersections(i1a,i2b) yIntersections(i1a,i2b)];
                pCornersImg(3,:) = [xIntersections(i1b,i2b) yIntersections(i1b,i2b)];
                pCornersImg(4,:) = [xIntersections(i1b,i2a) yIntersections(i1b,i2a)];
                
                % Make sure that points are in clockwise order.
                % If not, exchange points 2 and 4.
                v12 = pCornersImg(2,:) - pCornersImg(1,:);
                v13 = pCornersImg(3,:) - pCornersImg(1,:);
                if v12(1)*v13(2) - v12(2)*v13(1) < 0
                    temp = pCornersImg(2,:);
                    pCornersImg(2,:) = pCornersImg(4,:);
                    pCornersImg(4,:) = temp;
                end
                
                % Let the points start from left top
                sortcorners = pCornersImg;
                sortcorners=sortcorners(:,1).^2+sortcorners(:,2).^2;
                index = find(sortcorners==min(sortcorners));
                if length(index)>1
                    index = index(1);
                end
                pCornersImg = circshift(pCornersImg, [-(index-1), 0]);
                
                % Fit a homography using those four points.
                T = fitgeotrans(pCornersImg,pCornersRef, 'projective');
                
                % Transform all intersections points to the reference image.
                pIntersections = transformPointsForward(T, ...
                    [xIntersections(:) yIntersections(:)]);
                                
                pIntersectionsRef = [xIntersectionsRef(:) yIntersectionsRef(:)];
                
                % For each predicted reference point, find the closest
                % detected image point.
                dPts = 1e6 * ones(size(pIntersectionsRef,1),1);
                for i=1:size(pIntersectionsRef,1)
                    x = pIntersectionsRef(i,1);
                    y = pIntersectionsRef(i,2);
                    d = ((x-pIntersections(:,1)).^2 + (y-pIntersections(:,2)).^2).^0.5;
                    dmin = min(d);
                    dPts(i) = dmin;
                end
                
                % If the distance is less than DMIN, count it as a match.
                nMatches = sum(dPts < DMIN);
                
                % Calculate the avg error of the matched points.
                avgErr = mean(dPts(dPts < DMIN));
                
                % Keep the best combination found so far, in terms of
                % the number of matches and the minimum error.
                if nMatches < nMatchesBest
                    continue;
                end
                
                if (nMatches == nMatchesBest) && (avgErr > avgErrBest)
                    continue;
                end
                
                % Got a better combination; save it.
                avgErrBest = avgErr;
                nMatchesBest = nMatches;
                corners = pCornersImg;
                
                Ewarp = imwarp(E, T, 'OutputView', ...
                    imref2d(size(zeros(max(xIntersectionsRef(:)),max(yIntersectionsRef(:)))),...
                    [1 max(xIntersectionsRef(:))], [1 max(yIntersectionsRef(:))]));
                
                if Plotimg
                    % Display the predicted and measured points.
                    figure(7),imshow(Ewarp,[]);
                    title('Predicted and measured points');
                    hold on
                    plot(pIntersections(:,1), pIntersections(:,2), 'yd');
                    hold off
                    
                    for k = 1:size(pIntersectionsRef,1)
                        rectangle('Position', [pIntersectionsRef(k,1)-DMIN pIntersectionsRef(k,2)-DMIN 2*DMIN 2*DMIN], ...
                            'Curvature', [1 1], 'EdgeColor', 'w', 'LineWidth', 2);
                    end
                    
                    newpCornersImg = transformPointsForward(T, ...
                        [pCornersImg(:,1) pCornersImg(:,2)]);
                    
                    rectangle('Position', [newpCornersImg(1,1)-10 newpCornersImg(1,2)-10 20 20], ...
                        'Curvature', [1 1], 'EdgeColor', 'r', 'LineWidth', 2);
                    rectangle('Position', [newpCornersImg(2,1)-10 newpCornersImg(2,2)-10 20 20], ...
                        'Curvature', [1 1], 'EdgeColor', 'g', 'LineWidth', 2);
                    rectangle('Position', [newpCornersImg(3,1)-10 newpCornersImg(3,2)-10 20 20], ...
                        'Curvature', [1 1], 'EdgeColor', 'b', 'LineWidth', 2);
                    rectangle('Position', [newpCornersImg(4,1)-10 newpCornersImg(4,2)-10 20 20], ...
                        'Curvature', [1 1], 'EdgeColor', 'y', 'LineWidth', 2);
                    
                end
                fprintf(' Found %d matches, average error = %f\n', ...
                    nMatchesBest, avgErrBest);
            end
            
        end
        
    end
    
end

end