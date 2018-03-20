function [xIntersections, yIntersections] = findIntersections(lines1, lines2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  lines1, lines2 - Two groups of lines
% Output: xIntersections, yIntersections - Coordinates of intersections
%
% Intersect every pair of lines, one from set 1 and one from set 2.
% Output arrays contain the x,y coordinates of the intersections of lines.
%   xIntersections(i1,i2): x coord of intersection of i1 and i2
%   yIntersections(i1,i2): y coord of intersection of i1 and i2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1 = size(lines1,2);
N2 = size(lines2,2);

xIntersections = zeros(N1,N2);
yIntersections = zeros(N1,N2);
for i1=1:N1
    % Extract rho, theta for this line
    r1 = lines1(2,i1);
    t1 = lines1(1,i1);
    
    % A line is represented by (a,b,c), where ax+by+c=0.
    % We have r = x cos(t) + y sin(t), or  x cos(t) + y sin(t) - r = 0.
    l1 = [cosd(t1); sind(t1); -r1];
    
    for i2=1:N2
        % Extract rho, theta for this line
        r2 = lines2(2,i2);
        t2 = lines2(1,i2);
        
        l2 = [cosd(t2); sind(t2); -r2];
        
        % Two lines l1 and l2 intersect at a point p where p = l1 cross l2
        p = cross(l1,l2);
        p = p/p(3);
        
        xIntersections(i1,i2) = p(1);
        yIntersections(i1,i2) = p(2);
    end
end

end
