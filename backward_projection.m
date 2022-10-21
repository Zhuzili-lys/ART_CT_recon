%Yusheng Liu Implementation of ray-driven backward projection
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
function reconImg = backward_projection(A,imgName, PROJECTION_SIZE,PROJECTION_ANGLE_SIZE)
arguments
    A = [];
    imgName = "window1.png"
    PROJECTION_SIZE = 40;
    PROJECTION_ANGLE_SIZE = 200;
end
    LEARNING_RATE = 1.05;
    THRESHOLD = 100;

    %regarding the gtound truth of the projection
    originalImg = imread(imgName);
    grayImg = rgb2gray(originalImg);
    imageSize = size(grayImg);
    DISTRIBUTION_TRUTH= double(reshape(grayImg,[imageSize(1)*imageSize(2), 1]));
    projectionSize = PROJECTION_SIZE*PROJECTION_ANGLE_SIZE;
    PROJECTION_TRUTH = A*DISTRIBUTION_TRUTH;
    distribution = ones(imageSize(1)*imageSize(2),1);
    
    %calculate the projection data
    %compare the projection with the projection ground truth data
    errorSum = THRESHOLD+1;
    iter = 0;
    errorList = [];

    while(errorSum>THRESHOLD)
        iter = iter + 1;
        errorSum = 0;
        projection = A*distribution;
        for projectionIndex = 1:PROJECTION_SIZE*PROJECTION_ANGLE_SIZE
            %iterate through the projection to give the error
            errorNew = (projection(projectionIndex)-PROJECTION_TRUTH(projectionIndex));
            errorNew = errorNew * errorNew;
            errorSum = errorSum + errorNew;    
        end
        errorSum = sqrt(errorSum);
        %subtract each row to get a_i
        for ii = 1:projectionSize
%             if(ii == 9951)
%                 ii;
%             end
            aRow = A(ii,:);
            aNorm = aRow*aRow';
            %pixel-wise correction of the distribution
            if(aNorm~=0)
                distribution = distribution-LEARNING_RATE*(aRow*distribution-PROJECTION_TRUTH(ii))*aRow'/aNorm;
            end
        end
        errorList = [errorList,errorSum];
        plot(1:iter,errorList);
%         pause(0.01);
        hold on;
        ITER = num2str(iter);
        ERROR = num2str(errorSum);
        sprintf(ITER+","+ERROR)
    end
    reconImg = reshape(distribution,imageSize);
    reconImg = uint8(reconImg);
    montage({grayImg,reconImg});
end

%make correction to the distribution


