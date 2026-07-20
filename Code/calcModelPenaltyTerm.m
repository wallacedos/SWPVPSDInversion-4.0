function PTerm = calcModelPenaltyTerm(epsilon,m,mReference)
% Summary of this function goes here.
% PTerm = calcModelPenaltyTerm(epsilon,m)
% Detailed explanation goes here.
% The function is for calculating the model penalty term of
% the inversion method of preconditioned steepest-descent (PSD)
% algorithm.
%
%   IN
%        epsilon: predefined hyper-parameter.
%              m: the model parameter of inversion,it consists of [S-wave
%                 velocity layer-thickness].
%      mReference: reference model, full 0 vector is most smooth
%      constraint,otherwise is the minimal L2 norm constraint between the
%      model and reference model.
%
%  OUT
%          PTerm:  the value of penalty term.
%
%  Author(s): Yan Yingwei
%  Copyright: 2022-2030
%  Revision:  1.0  Date: 11/28/2022
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

ll = length(m);
n = round((ll+1)/2+0.1);

s1 = sum(mReference,"all");

eeps = 1e-6;
if s1<eeps
    
    if n<3
        PTerm = 0.0;
        return;
    end
    
    vsDifference = m(1:n);
    hDifference = m(n+1:end);
    
    % gradient of S-wave velocity model
    gradVSDifference = zeros(1,n-1);
    for i=1:n-1
        gradVSDifference(i) = (vsDifference(i+1)-vsDifference(i))/hDifference(i);
    end
    
    % second derivatives of S-wave velocity model
    secondDiffVS = zeros(1,n-2);
    for i=1:n-2
        secondDiffVS(i) = (gradVSDifference(i+1)-gradVSDifference(i))/(0.5*hDifference(i)+0.5*hDifference(i+1));
    end
    
    PTerm = epsilon*norm(secondDiffVS')/sqrt(n-2);
else
    vs = m(1:n);
    h = m(n+1:end);
    meanVs = mean(vs);
    meanH = mean(h);
    r1 = meanVs/meanH;
    m = [vs r1*h];
    
    vsR = mReference(1:n);
    hR = mReference(n+1:end);
    meanVsR = mean(vsR);
    meanH = mean(hR);
    r2 = meanVsR/meanH;
    mReference = [vsR r2*hR];
    
    deltaM = m-mReference;
    PTerm = 0.5*epsilon*norm(deltaM(:))/(2*n-1);
end
end

