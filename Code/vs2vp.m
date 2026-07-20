function vp = vs2vp(vs,vpdvs)
%   Summary of this function goes here.
%   vp = vs2vp(vs,vpdvs)
%   Detailed explanation goes here.
%   The function is for getting the P-wave velocity according the ratio
%   value vector of VP/VS (input parameter vpdvs)
%
%   IN      
%        vs: the S-wave velocity vector, a row vector.
%     vpdvs: the ratio vector between the P-wave velocity and S-wave
%            velocity,it's a row vector.
%                
%
%  OUT   
%        vp: P-wave velocity vector
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

vp = vpdvs.*vs;
end

