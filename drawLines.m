function drawLines(rhos, thetas, imageSize, color)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  rhos, thetas -  Hough transform parameters of the lines
%         imageSize - Image size
%         color - color of the line
% Output: None
%
% This function draws lines on whatever image is being displayed.
% Input parameters:
%   rhos,thetas: representation of the line (theta in degrees)
%   imageSize:  [height,width] of image being displayed
%   color:  color of line to draw
%
% Equation of the line is rho = x cos(theta) + y sin(theta), or
%    y = (rho - x*cos(theta))/sin(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:length(thetas)
    if abs(thetas(i)) > 45
        % Line is mostly horizontal.  Pick two values of x,
        % and solve for y = (-ax-c)/b
        x0 = 1;
        y0 = (-cosd(thetas(i))*x0+rhos(i))/sind(thetas(i));
        x1 = imageSize(2);
        y1 = (-cosd(thetas(i))*x1+rhos(i))/sind(thetas(i));
    else
        % Line is mostly vertical. Pick two values of y,
        % and solve for x = (-by-c)/a
        y0 = 1;
        x0 = (-sind(thetas(i))*y0+rhos(i))/cosd(thetas(i));
        y1 = imageSize(1);
        x1 = (-sind(thetas(i))*y1+rhos(i))/cosd(thetas(i));
    end
    
    line([x0 x1], [y0 y1], 'Color', color);
    text(x0,y0,sprintf('%d', i), 'Color', color);
end

end