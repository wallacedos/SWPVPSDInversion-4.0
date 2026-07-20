function [vsSta, depthSta] = calcStairsData(vs,h, maxDepth)
%   Summary of this function goes here.
%   [vsSta, depthSta] = calcStairsData(vs,h, maxDepth)
%   Detailed explanation goes here.
%   The function is for calculating the stairs data of the input S-wave
%   velocity model.
%
%   IN
%             vs: the S-wave velocity of layered model, it is a row vector.
%              h: the layer thickness of layered model, it is a row vector.
%       maxDepth: the maximum depth of the output stairs data of the output
%                 S-wave velocity.
%
%
%
%  OUT
%          vsSta: the output stairs data of the S-wave velocity.
%       depthSta: the output stairs data of the depth.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2024-2030
%  Revision:  1.0  Date: 2/26/2024
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

if nargin==2
   maxDepth = 1.5*sum(h);
end

n = length(vs);
vs = [vs(1) vs];
depth = zeros(1,n+1);
for i=2:n
    depth(i) = sum(h(1:i-1));
end
depth(end) = maxDepth;
[vsSta,depthSta] = stairs(vs,depth);
end

