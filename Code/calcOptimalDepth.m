function optimalDepth = calcOptimalDepth(ObsInfo,ModelInfo,maxSearchDepth,misfitType,waveType,num)
%   Summary of this function goes here.
%   optimalDepth = calcOptimalDepth(ObsInfo,ModelInfo,maxSearchDepth,misfitType,waveType,num)
%   Detailed explanation goes here.
%   The function is for calculating the optimal depth for the initial
%   S-wave velocity model.A simple linear search method is adopted.
%   
%
%   IN      
%        ObsInfo: the struct about the observed information, it consists of
%                 following fields, 'f' (frequency axis vector), 'pv'
%                 (observed phase velocity values matrix, each column represents
%                 a Rayleigh wave mode, 0-elements are the unobserved data point)
%                 'Wd' (it shares the same shape as the observed data pv),
%                 'modeSN' (it is a row vector, it represents the Rayleigh 
%                 wave mode numerical order).
%      ModelInfo: the struct about the model information during the inversion,
%                 it consists of following fields,'Bound' (searching bound
%                 for each paramter, it's n*2 matrix, n is the number of
%                 reconstructed parameter), 'Ini' (initial model of the inverion)
%                 'den' (density of the inversion, it remains a constant),
%                 'vpdvs' (the ratio of P-wave and S-wave velocity for each
%                 layer, it's also a row vector.),'vpdvs' only for rayleigh
%                 wave.
% maxSearchDepth: the maximum search depth for the initial model.
%     misfitType: the misfit function type string, "L1" and "L2" can be selected.
%       waveType: the surface wave type string,'rayleigh' and 'love' is
%                 selected.
%            num: the number of the trial depth point.
%    
%                
%  OUT   
% optimalDepth: the optimal detecting depth.
% 
%
%  Author(s): Yan Yingwei
%  Copyright: 2021-2030 
%  Revision:  1.0  Date: 11/23/2021
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

if nargin==5
   num = 50;
elseif nargin==4
   num = 50;
   waveType = "rayleigh";
elseif nargin==3
   num = 50;
   waveType = "rayleigh";
   misfitType = "L2";
end

initialModel = ModelInfo.Ini;
den = ModelInfo.den;
numLayer = floor((length(initialModel)+0.1+1)/2);

f = ObsInfo.f;
pvObserved = ObsInfo.pv;
modeSN = ObsInfo.modeSN;

[M,N] = size(pvObserved);
maxModeSN = max(modeSN);

WdState = isfield(ObsInfo,"Wd");
if WdState
    Wd = ObsInfo.Wd;
else
    Wd = ones(M,N);
end

VS = initialModel(1:numLayer);
trialDepth = linspace(0.01*maxSearchDepth,maxSearchDepth,num);
fValue = zeros(num,1);
pvPredicted = zeros(M,N);

if strcmp(waveType,"rayleigh")
    vpdvsState = isfield(ModelInfo,"vpdvs");
    if vpdvsState
        vpdvs = ModelInfo.vpdvs;
    else
        vpdvs = 2.45*ones(1,numLayer);
    end
    for i=1:num
        pvPredicted(:,:) = 0;
        H = trialDepth(i)/(numLayer-1)*ones(1,numLayer-1);
        try
            pvPredictedCal = calcmulti(f,VS,H, vpdvs.*VS,den);
            [M1,N1] = size(pvPredictedCal);
            
            if (N1<maxModeSN)
                pvPredicted(1:M1,1:N1) = pvPredictedCal;
            else
                pvPredicted(1:M1,1:maxModeSN) = pvPredictedCal(:,1:maxModeSN);
            end
            fValue(i) = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,misfitType);
        catch
            fValue(i) = Inf;
        end
    end
    [~,I] = min(fValue);
    optimalDepth = trialDepth(I);
elseif strcmp(waveType,"love")
    for i=1:num
        pvPredicted(:,:) = 0;
        H = trialDepth(i)/(numLayer-1)*ones(1,numLayer-1);
        try
            pvPredictedCal = calclovemulti(f,VS,H,den);
            [M1,N1] = size(pvPredictedCal);
            if (N1<maxModeSN)
                pvPredicted(1:M1,1:N1) = pvPredictedCal;
            else
                pvPredicted(1:M1,1:maxModeSN) = pvPredictedCal(:,1:maxModeSN);
            end
            fValue(i) = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,misfitType);
        catch
            fValue(i) = Inf;
        end
    end
    [~,I] = min(fValue);
    optimalDepth = trialDepth(I);
end
end

