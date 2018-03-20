function [xIntersectionsRef, yIntersectionsRef] = createReference(sizeRef,Plotimg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  sizeRef - The size of the reference (sizeRef+1) x (sizeRef+1)
%         Plotimg - Whether to show the intermediate steps
% Output: xIntersectionsRef - The x-coordinate of the reference
%         yIntersectionsRef - The y-coordinate of the reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizeSquare = sizeRef/8;     % size of one square

% Predict all line intersections.
[xIntersectionsRef, yIntersectionsRef] = meshgrid(1:9, 1:9);
xIntersectionsRef = (xIntersectionsRef-1)*sizeSquare + 1;
yIntersectionsRef = (yIntersectionsRef-1)*sizeSquare + 1;

if Plotimg
    % Draw reference image.
    Iref = zeros(sizeRef+1, sizeRef+1);
    figure(6), imshow(Iref), title('Reference image');
    
    % Show all reference image intersections.
    hold on
    plot(xIntersectionsRef, yIntersectionsRef, 'y+');
    hold off
end
end