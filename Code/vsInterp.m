function [vsi,depthi] = vsInterp(vs,depth,dh,maxdepth)
%   Summary of this function goes here.
%   [vsi,depthi] = vsInterp(vs,depth,dh,maxdepth)
%   Detailed explanation goes here.
%   The function is for getting the 1-D S-wave velocity after
%   interpolation.
%
%   IN      
%        vs: the S-wave velocity vector, a row vector.
%     depth: the depth vector of the input S-wave velocity, the first
%            element is 0, m.
%        dh: the interval (m) of the interpolation.
%  maxdepth: the max depth of the interpolation
%                
%
%  OUT   
%       vsi: the S-wave velocity vector after the interpolation.
%    depthi: the depth vecotor of the output S-wave velocity.
% 
% 
%  Author(s): Yan Yingwei
%  Copyright: 2021-2030 
%  Revision:  1.0  Date: 11/22/2021
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

depthq = 0:dh:depth(end);
vsq = interp1(depth,vs,depthq,"linear");
depthi = 0:dh:maxdepth;
vsi = zeros(length(depthi),1);
if depthq(end)>maxdepth
    l = length(depthi);
    vsi = vsq(1:l);
else
    l = length(depthq);
    vsi(1:l) = vsq;
    vsi(l+1:end)=vsq(end);
end
end

