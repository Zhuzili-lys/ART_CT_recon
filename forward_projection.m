%Yusheng Liu Implementation of ray-driven forward projection
%Assumptions: 
%1. Linear Interpolation
%2. Average distribution of picture in each pixel
%3. The rotation of beams is concentric with the region to be reconstructed in the
%   beginning
%4. The ray at the center always passes the rotation center of the
%   projections
%5. The distance between pixel centers is 1
%6. At least 2 sampling beams

%Initialization
% IMAGE_PATH = "dataPath";
% IMAGE_NAME= {"Image001.tif"};
% PROJECTION_PATH = "projectionPath";
% PROJECTION_NAME = {"Projection.tif"};
% 
% SINO_PATH = "sinogramPath";
% RESULT_PATH = "result";

function A = forward_projection(IMAGE_NAME, PIXEL_DISTANCE, PROJECTION_SIZE, PROJECTION_ANGLE_SIZE, PROJECTION_RANGE, RAY_DISTANCE)

    %Parameter of vehicles
    arguments
        IMAGE_NAME = "window1.png";
        PIXEL_DISTANCE = 0.6;
        PROJECTION_SIZE = 40; %number of rays
        PROJECTION_ANGLE_SIZE = 200;
        PROJECTION_RANGE = [-pi/4, 3*pi/4];
        RAY_DISTANCE = 1;
    end
    
    originalImg = imread(IMAGE_NAME);
    grayImg = rgb2gray(originalImg);
    IMAGE_SIZE = size(grayImg);
    angleRangeLength = abs(PROJECTION_RANGE(2)-PROJECTION_RANGE(1));
    projectionAngleList = linspace(PROJECTION_RANGE(1)+angleRangeLength/PROJECTION_ANGLE_SIZE,PROJECTION_RANGE(2),PROJECTION_ANGLE_SIZE); %angle range -pi/4 to 3*pi/4
    imageCenter = IMAGE_SIZE/2;
    %Ray-driven forward projection with linear interpolation
    %iterate through all the projections
    pixelSize = IMAGE_SIZE(1)*IMAGE_SIZE(2);
    projectionSize = PROJECTION_ANGLE_SIZE*PROJECTION_SIZE;
    A = zeros(projectionSize, pixelSize); %pre-allocate

    for angleIndex = 1:PROJECTION_ANGLE_SIZE
        theta = projectionAngleList(angleIndex);
        %iterating through the rays of the same angle
        for rayIndex = 1:PROJECTION_SIZE
            distance = (rayIndex-PROJECTION_SIZE/2)*RAY_DISTANCE;
            %when the angle is no larger than pi/4, track first zero of x
            if(theta>-pi/4 && theta <= pi/4)
                for fixedVar = 1:IMAGE_SIZE(1)
                    [ptrZero1,ptrZero2] = searchZero((1:IMAGE_SIZE(2)), 1, IMAGE_SIZE(2), theta, distance, imageCenter, fixedVar, true); 
                    %if we have found interception point
                    if(ptrZero1~=-1)
                        lineValue1 = abs(evaluate_line(theta, distance, imageCenter, ptrZero1, fixedVar));
                        lineValue2 = abs(evaluate_line(theta, distance, imageCenter, ptrZero2, fixedVar));
                        projectionPosition = (angleIndex-1)*PROJECTION_SIZE+rayIndex;
                        pixelPosition1 = ptrZero1+(fixedVar-1)*IMAGE_SIZE(2);
                        pixelPosition2 = ptrZero2+(fixedVar-1)*IMAGE_SIZE(2);
                        if(lineValue1+lineValue2~=0)
                            A(projectionPosition, pixelPosition1) =  lineValue2/(lineValue1+lineValue2);
                            A(projectionPosition, pixelPosition2) =  lineValue1/(lineValue1+lineValue2);  
                        else
                            A(projectionPosition, pixelPosition1) =  1;
                            A(projectionPosition, pixelPosition2) =  1;
                        end
                    end
                end
            %when angle between pi/4 and pi*3/4, track y
            elseif(theta>pi/4 && theta<pi*3/4)
                for fixedVar = 1:IMAGE_SIZE(2)
                    [ptrZero1,ptrZero2] = searchZero((1:IMAGE_SIZE(1)),1, IMAGE_SIZE(1), theta,distance,imageCenter,fixedVar,false);
                    %if interception exists
                    if(ptrZero1~=-1)
                        lineValue1 = abs(evaluate_line(theta, distance, imageCenter, fixedVar, ptrZero1));
                        lineValue2 = abs(evaluate_line(theta, distance, imageCenter, fixedVar, ptrZero2));
                        projectionPosition = (angleIndex-1)*PROJECTION_SIZE+rayIndex;
                        pixelPosition1 = (ptrZero1-1)*IMAGE_SIZE(2)+fixedVar;
                        pixelPosition2 = (ptrZero2-1)*IMAGE_SIZE(2)+fixedVar;
                        if(lineValue1+lineValue2~=0)
                            A(projectionPosition, pixelPosition1) = lineValue2/(lineValue1+lineValue2);
                            A(projectionPosition, pixelPosition2) =  lineValue1/(lineValue1+lineValue2);
                        else
                            A(projectionPosition, pixelPosition1) = 1;
                            A(projectionPosition, pixelPosition2) = 1;
                        end
                    end
                end
            end
        end
    end
end

function [ptrZero1, ptrZero2]=searchZero(array, ptr1, ptr2, theta, d, center, fixedVar, isHorizontal)
    middle = floor((ptr1+ptr2)/2);
    if(isHorizontal)
        lineValue1 = evaluate_line(theta,d, center, ptr1, fixedVar);
        lineValue2 = evaluate_line(theta,d, center, ptr2, fixedVar);
        lineValueMid = evaluate_line(theta, d, center, middle, fixedVar);
    elseif(~isHorizontal)
        lineValue1 = evaluate_line(theta,d, center, fixedVar,ptr1);
        lineValue2 = evaluate_line(theta,d, center, fixedVar,ptr2);
        lineValueMid = evaluate_line(theta, d, center, fixedVar, middle);
    end
    if(lineValue1*lineValue2>0)
        ptrZero1 = -1;
        ptrZero2 = -1;
    elseif(lineValue1 == 0)
        ptrZero1 = ptr1;
        ptrZero2 = ptr1;
    elseif(lineValue2 == 0)
        ptrZero1 = ptr2;
        ptrZero2 = ptr2;
    elseif(ptr2 == ptr1+1 || ptr2==ptr1)
        ptrZero1 = ptr1;
        ptrZero2 = ptr2;
    elseif(lineValue1*lineValueMid<=0)
        [ptrZero1, ptrZero2] = searchZero(array, ptr1, middle, theta, d, center, fixedVar,isHorizontal);
    elseif(lineValue2*lineValueMid<=0)
        [ptrZero1, ptrZero2] = searchZero(array, middle, ptr2, theta, d, center, fixedVar,isHorizontal);
    end
end

%function giveing the signed distance to the line specified by theta and d
%from the center
function lineValue = evaluate_line(theta, d, center, x, y)
    lineValue = sin(theta)*(y-center(2))-cos(theta)*(x-center(1))-d;
end