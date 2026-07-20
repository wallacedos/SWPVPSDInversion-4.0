function pv = calcmulti(f,VS,H,VP,den,n_mode)
%   Summary of this function goes here.
%   pv = calcmulti(f,VS,H,VP,den,n_mode)
%   Detailed explanation goes here.
%   The function is for calculating multi-mode of Rayleigh wave.
%
%   IN    f: row vector of frequency.
%        VS: row vector of shear wave velocity.
%         H: row vector of thickness of layer.
%        VP: row vector of primary wave velocity.
%       den: row vector of density.
%    n_mode: the mode number need to compute.
%
%  OUT   pv: the phase velocity of multi-mode Rayleigh wave.
%
%  Example:
%  f=5:100;VS=[200 400];VP=[400 800];H=[10];den=[2000 2000];
%  pv=calcmulti(f,VS,H,VP,den);
%
%  Author(s): Yan Yingwei
%  Copyright: 2017-2030
%  Revision: 1.0  Date: 2/27/2017
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

if nargin<=5
    n_mode = 1;
end
if(nargin<=4)
    [~,c]=size(VS);
    den=2*ones(1,c);
end
if(nargin<=3)
    VP=2*VS;
end

cmax=max(VS);
cmin=0.88*min(VS);
dc1=(cmax-cmin)/400;
cc=cmin:dc1:cmax;

[~,nf]=size(f);
m = [VS,H,VP,den];

pv = zeros(nf,n_mode);
parfor i=1:nf
    ccc = RayleighSearchRoot(f(i),m,cc,n_mode);
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
