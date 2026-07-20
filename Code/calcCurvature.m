function [curv] = calcCurvature(x,y)
%   Summary of this function goes here.
%   [curv] = calcCurvature(x,y)
%   Detailed explanation goes here.
%   The function is for calculating the curvature.
%
%  In
%      x: vector of x value;
%      y: vector of y value;
%
%  Out
%
%    curv: curvature vector
%
%  Author(s): Yan Yingwei
%  Copyright: 2022-2030 
%  Revision:  1.0  Date: 11/28/2022
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

l = length(x);
curv = zeros(1,l-2);
for i=1:l-2
    a = x(i:i+2);
    b = y(i:i+2);
    [kappa,~] = PJcurvature(a,b);
    curv(i) = kappa;
end
end

