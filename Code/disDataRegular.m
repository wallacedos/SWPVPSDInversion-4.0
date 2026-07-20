function  pv = disDataRegular(pv,maxModeSN)
%   Summary of this function goes here.
%   pv = disDataRegular(pv,maxModeSN)
%   Detailed explanation goes here.
%   The function is for regularizing the dispersion data by the maximum
%   mode-order.
%
%   IN
%             pv: the phase velocity values matrix, each column represents 
%                 a mode,.
%      maxModeSN: the maximum mode-order.
%
%
%
%  OUT
%             pv: the output phase velocity values matrix, it is
%                 regularized.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2024-2030
%  Revision:  1.0  Date: 3/1/2024
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).


[M1,N1] = size(pv);
pvPredicted = zeros(M1,maxModeSN);
if (N1<maxModeSN)
    pvPredicted(1:M1,1:N1) = pv;
else
    pvPredicted(1:M1,1:maxModeSN) = pv(:,1:maxModeSN);
end
pv = pvPredicted;
end

