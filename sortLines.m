function lines = sortLines(lines, sizeImg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  lines - The unsorted input lines
%         sizeImg - The size of the image
% Output: lines - The sorted lines
%
% If the lines are mostly horizontal, sort on vertical distance from yc.
% If the lines are mostly vertical, sort on horizontal distance from xc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xc = sizeImg(2)/2;  % Center of image
yc = sizeImg(1)/2;

t = lines(1,:);     % Get all thetas
r = lines(2,:);     % Get all rhos

% If most angles are between -45 .. +45 degrees, lines are mostly
% vertical.
nLines = size(lines,2);
nVertical = sum(abs(t)<45);
if nVertical/nLines > 0.5
    % Mostly vertical lines.
    dist = (-sind(t)*yc + r)./cosd(t) - xc;  % horizontal distance from center
else
    % Mostly horizontal lines.
    dist = (-cosd(t)*xc + r)./sind(t) - yc;  % vertical distance from center
end

[~,indices] = sort(dist, 'ascend');
lines = lines(:,indices);

end
