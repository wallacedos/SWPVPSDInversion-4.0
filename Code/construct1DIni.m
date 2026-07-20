function mIni = construct1DIni(f,pvObs,depth,alp,mId)
%   Summary of this function goes here.
%   mIni = construct1DIni(f,pvObs,depth,alp, mId)
%   Detailed explanation goes here.
%   The function is for generating the initial S-wave velocity model from
%   the observed Rayleigh-wave data. The initial S-wave velocity model can
%   be utilized for the inversion of surface-wave dispersion curves.
%
%   IN
%              f: the frequency vector of the observed Rayleigh-wave
%                 dispersion curve, it is a row vector.
%          pvObs: the obseved phase velocities of the fundamental-mode
%                 Rayleigh-wave, it is a vector.
%          depth: the depth vector of the initial S-wave velocity, it is a 
%                 row vector.
%            alp: the empirical parameter of the penetrating depth of
%                 Rayleigh-wave, default as 0.5 times of wavelength.
%            mId: the selection for describing the construction method of
%                 initial model, 0 means the common construction method
%                 based on the wavelength, 1 means the construction method
%                 based on the travel-time.
%
%
%  OUT
%            mIni: the initial S-wave velocity model, such as
%                  [S-wave velocity, layer thickness], it is a row vector.
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

if nargin==4
    mId = 0;
end
if nargin==3
    alp = 0.5;
    mId = 0;
end

lambdaFundamental = pvObs(:)./(f(:));
[Y,I] = sort(lambdaFundamental);
x = alp*Y;
lDepth = length(depth);
mIni = zeros(1,lDepth);
l_x = length(x);

if mId==0
    mIniSeed = pvObs(I)./0.88;
end
if mId==1
    if l_x~=1
        mIniSeed = zeros(l_x,1);
        for i=1:l_x-1
            t2 = 1/f(i+1); t1 = 1/f(i);
            v2 = pvObs(i+1); v1 = pvObs(i);
            VS_i  = abs(((t2*v2^4-t1*v1^4)/(t2-t1)))^(1/4);
            mIniSeed(i) = VS_i;
        end
        mIniSeed(l_x) = mIniSeed(l_x-1);
        mIniSeed = mIniSeed(I);
    end
end

if (l_x)~=1
    loc1 = find(depth>x(1));
    if ~isempty(loc1)
        loc1 = loc1(1);
        loc2 = find(depth<x(end));
        loc2 = loc2(end);
        xq = depth(loc1:loc2);
        xq = xq(:);
        mIni(loc1:loc2) = interp1(x,mIniSeed,xq,"PCHIP");
        mIni(1:loc1-1) = mIni(loc1);
        mIni(loc2+1:end) = mIni(loc2);
    else
        mIni = mIniSeed(1)*ones(1,lDepth);
    end
else
    mIni = pvObs(1)/0.88*ones(1,lDepth);
end
end

