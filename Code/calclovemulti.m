function pv = calclovemulti(f,VS,H,den,n_mode)
%   Summary of this function goes here.
%   pv = calclovemulti(f,VS,H,den,n_mode)
%   Detailed explanation goes here.
%   The function is for calculating multi-mode of Love wave.
%
%   IN    
%         f: row vector of frequency£¬ Hz.
%        VS: row vector of shear wave velocity, m/s.
%         H: row vector of thickness of layer, m/s.
%       den: row vector of density, m/s.
%    n_mode: the mode number needs to compute.
%
%
%  OUT   
%        pv: the phase velocities Love waves, m/s.
%
%  Example:
%  f=5:100;VS=[200 400];H=[10];den=[2000 2000];
%  pv=calclovemulti(f,VS,H,den);
%
% History:
%          1) Use parallel computing to solve for the phase velocitie of 
%          Love waves. 4/9/2025
%
%  Author(s): Yan Yingwei
%  Copyright: 2017-2030
%  Revision: 2.0  Date: 9/4/2025
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

if nargin<=4
    n_mode = 1;
end

if(nargin<=3)
    [~,c]=size(VS);
    den=2000*ones(1,c);
end

cmax=max(VS);
cmin=min(VS);
dc1=(cmax-cmin)/400;
cc=cmin:dc1:cmax;

[~,nf]=size(f);
m = [VS,H,den];

pv = zeros(nf,n_mode);
parfor i=1:nf
    ccc = LoveSearchRoot(f(i),m,cc,n_mode);
    pv(i,:) = ccc;
end
if n_mode~=1
    [~,N] = size(pv);
    efficientModeNum = 1;
    for k=N:-1:1
        if sum(pv(:,k),"all")~=0
            efficientModeNum = k;
            break;
        end
    end
    pv = pv(:,1:efficientModeNum);
end
end

