function [g,Hess,HessDiagElement] = calcLoveDCInvGradHess(ObsInfo,ModelInfo,m0,fValue0,epsilon, mReference,misfitType, stateId)
%   Summary of this function goes here.
%   [g,Hess,HessDiagElement] = calcLoveDCInvGradHess(ObsInfo,ModelInfo,m0,fValue0,epsilon, mReference,misfitType, stateId)
%   Detailed explanation goes here.
%   The function is for calculating the gradient and Hessian matrix of
%   the misfit function of love wave dispersion inversion.
%   For common computing mode:
%   The simple difference is used to calculate the gradient and Hessian
%   matrix of the misfit function. The Hessian matrix is assumed to be a
%   diagonally-dominant matrix.
%   For fast computing mode:
%   The gradient and Hessian matrix of the misfit function is calculated by
%   the implicit function theorem.
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
%                 'den' (density of the inversion, it remains a constant).
%             m0: the current model point, it's a row vector.
%        fValue0: the misfit function value of the current model point m0.
%        epsilon: the hyper-parameter for model punishment.
%      mReference: reference model, full 0 vector is most smooth
%                  constraint,otherwise is the minimal L2 norm constraint 
%                  between the model and reference model.
%     misfitType: the misfit function type string, "L1" and "L2" can be selected.
%        stateId: the string for describing the computing mode, "Fast"
%                 means the fast computing mode, "Common" means the common
%                 computing mode.
%
%
%  OUT
%              g: the gradient of the misfit function at current model point m0.
%           Hess: the Hessian matrix of the misfit function at current model point
%                 m0,it is considered as a diagonally-dominant matrix.
% HessDiagElement:the diagonal elements of the Hessian matrix.
%
%
% History:
%          1£© Adding the fast computing method, the computational
%          performance of the inversion program has been significantly
%          improved, this enhancement is achieved through the application
%          of the implicit function theorem, 8/30/2025
%
%  Author(s): Yan Yingwei
%  Copyright: 2021-2030
%  Revision:  2.0  Date: 8/30/2025
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

modeSN = ObsInfo.modeSN;
f = ObsInfo.f;
pvObserved = ObsInfo.pv;
Wd = ObsInfo.Wd;
den = ModelInfo.den;

maxModeSN = max(modeSN);
n = length(m0);
g = zeros(n,1);
Hess = zeros(n,n);
numLayer = floor((n+1+0.2)/2);
delta = 0.005;
HessDiag = zeros(n,1);
if strcmp(stateId,"Common")
    parfor i=1:n
        mTemp = m0;
        mTemp(i) = m0(i)+delta*m0(i);
        VSTemp = mTemp(1:numLayer);
        HTemp = mTemp(numLayer+1:end);
        pvPredicted = calclovemulti(f,VSTemp,HTemp,den,maxModeSN);
        pvPredicted = disDataRegular(pvPredicted ,maxModeSN);
        fValueR = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,epsilon,mTemp,mReference,misfitType);
        
        mTemp = m0;
        mTemp(i) = m0(i)-delta*m0(i);
        VSTemp = mTemp(1:numLayer);
        HTemp = mTemp(numLayer+1:end);
        pvPredicted = calclovemulti(f,VSTemp,HTemp,den,maxModeSN);
        pvPredicted = disDataRegular(pvPredicted,maxModeSN);
        fValueL = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,epsilon,mTemp,mReference,misfitType);
        g(i) = (fValueR-fValueL)/(2*delta*m0(i));
        HessDiag(i) = (fValueR+fValueL-2*fValue0)/(delta*m0(i)*delta*m0(i));
    end
    HessDiagElement = HessDiag;
end

if strcmp(stateId,"Fast")  && strcmp(misfitType,"L2")
    modeNum = length(modeSN);
    IndCell = ObsInfo.IndCell;
    pvPred0 = ObsInfo.pvPred0;
    WdSqaure_reci = ObsInfo.WdSquare_reci;
    [JMatrix,JMatrix_trans] = calcLoveJacobiMatrix(ObsInfo,ModelInfo,m0);
    fValueModel0 = calcModelPenaltyTerm(epsilon,m0,mReference);
    fValueData0 = fValue0-fValueModel0;
    
    J_transWdDeltaPV = zeros(n,1);
    
    for k=1:modeNum
        J_transWdDeltaPV = J_transWdDeltaPV+JMatrix_trans(:,IndCell{k},k)*(Wd(IndCell{k},k).*(pvObserved(IndCell{k},k)-pvPred0(IndCell{k},modeSN(k))));
    end
    g = -WdSqaure_reci/(fValueData0)*J_transWdDeltaPV;
    
    J_transWdJ = zeros(n,n);
    for k=1:modeNum
        WdMatrix = diag(Wd(IndCell{k},k));
        J_transWdJ = J_transWdJ+JMatrix_trans(:,IndCell{k},k)*WdMatrix*JMatrix(IndCell{k},:,k);
    end
    
    J_transWdJ_diag = diag(J_transWdJ);
    
    HessDiagElement = WdSqaure_reci/(fValueData0*fValueData0)*g.* J_transWdDeltaPV+WdSqaure_reci/(fValueData0)*J_transWdJ_diag;
    
    g_m = zeros(n,1);
    H_m = zeros(n,1);
    for i=1:n
        m = m0;
        m(i) = m(i)+delta*m(i);
        fValueModel0_right = calcModelPenaltyTerm(epsilon,m,mReference);
        m = m0;
        m(i) = m(i)-delta*m(i);
        fValueModel0_left = calcModelPenaltyTerm(epsilon,m,mReference);
        g_m(i) = (fValueModel0_right-fValueModel0_left)/(2*delta*m0(i));
        H_m(i) = (fValueModel0_right+fValueModel0_left-2*fValueModel0)/(delta*m0(i)*delta*m0(i));
    end
    g = g+g_m;
    HessDiagElement = HessDiagElement+ H_m;
end
end

