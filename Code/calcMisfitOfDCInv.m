function fValue = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,epsilon,m,mReference,misfitType)
%   Summary of this function goes here.
%   fValue = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,epsilon,m,misfitType)
%   Detailed explanation goes here.
%   The function is for calculating the misfit function value.
%
%
%   IN
%     pvObserved: the obeserved phase velocity matrix, each column
%                 represents the surface wave phase velocity of specific
%                 mode,the observed phase velocity values of all mode share
%                 the same frequency axis, the 0-elements in the matrix are
%                 the unobserved data point.
%    pvPredicted: it is the predicted phase velocity value calculated from
%                 the trial model.
%         modeSN: the numerical order vector of all surface wave mode, it's
%                 a row vector, the order counts from 1, "1" represents the
%                 fundamental mode.
%             Wd: the weighted matrix for the observed phase velocity
%                 values, it has the same shape with pvObserved.
%        epsilon: the hyper-parameter for model punishment.
%              m: the model parameter of inversion,it consists of [S-wave
%                 velocity layer-thickness].
%      mReference: reference model, full 0 vector is most smooth
%      constraint,otherwise is the minimal L2 norm constraint between the
%      model and reference model.
%     misfitType: the misfit function type string, "L1" and "L2" can be selected.
%
%
%  OUT
%         fValue: the misfit function value.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2021-2030
%  Revision:  2.0  Date: 11/28/2022
%
%       History
%          1) Adding the most smooth model difference £¨the difference of
%          reference model and the model in iteration) punishment to the
%          calculation of misfit function value,
%          11/28/2022.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

if nargin==4
    misfitType = "L2";
elseif nargin==3
    misfitType = "L2";
    [M,N] = size(pvObserved);
    Wd = ones(M,N);
end

fValue = 0;
weightSum = 0;
[~,N] = size(pvObserved);
if strcmp(misfitType,"L2")
    for j=1:N
        [I,~] = find(pvObserved(:,j)~=0);
        weightSum = weightSum+sum(Wd(I,j).*Wd(I,j));
        fValue = fValue+sum((pvPredicted(I,modeSN(j))-pvObserved(I,j)).*(pvPredicted(I,modeSN(j))-pvObserved(I,j)).*Wd(I,j).*Wd(I,j));
    end
    fValue = sqrt(fValue/weightSum);
    PTerm = calcModelPenaltyTerm(epsilon,m,mReference);
    fValue = fValue+PTerm;
elseif strcmp(misfitType,"L1")
    for j=1:N
        [I,~] = find(pvObserved(:,j)~=0);
        weightSum = weightSum+sum(Wd(I,j));
        fValue = fValue+sum(abs(pvPredicted(I,modeSN(j))-pvObserved(I,j)).*Wd(I,j));
        fValue = fValue/weightSum;
    end
    PTerm = calcModelPenaltyTerm(epsilon,m,mReference);
    fValue = fValue+PTerm;
end
end

