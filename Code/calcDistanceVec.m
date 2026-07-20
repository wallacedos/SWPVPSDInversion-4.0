function distanceVec = calcDistanceVec(staCoor,mId)
%   Summary of this function goes here.
%   distanceVec = calcDistanceVec(staCoor,mId)
%   Detailed explanation goes here.
%   The function is for calculating the distance vector from the
%   coordinates of the seismic stations.
%
%
%   IN      
%       staCoor: describes the coordinates of all stations (m), its size 
%                is N*axisNum, axisNum is the number of the axis of the 
%                N is the number of the seismic stations.
%           mId: If it is "FirstRef", only the first station is refered;
%                "Pair": the distance between any two-staion pairs
%              
%  OUT   
%   distanceVec: the distance vector, m.
%
%  Author(s): Yan Yingwei
%  Copyright: 2023-2030 
%  Version: 1.0  Date: 9/8/2023
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

if nargin==1
    mId = "FirstRef";
end

[N,axisNum] = size(staCoor);

if strcmp(mId,"FirstRef")
    distanceVec = zeros(1,N-1);
    for j=2:N
        disSquare = 0.0;
        for i=1:axisNum
            dis = staCoor(j,i)-staCoor(1,i);
            disSquare = disSquare+dis*dis;
        end
        distanceVec(j-1) = sqrt(disSquare);
    end
elseif strcmp(mId,"Pair")
    num = round(((N*(N-1))/2+0.1));
    distanceVec = zeros(1,num);
    count = 0;
    for k=1:N-1
        for j=k+1:N
            disSquare = 0;
            for i=1:axisNum
                dis = staCoor(j,i)-staCoor(k,i);
                disSquare = disSquare+dis*dis;
            end
            count = count+1;
            distanceVec(count) = sqrt(disSquare);
        end
    end
end
end

