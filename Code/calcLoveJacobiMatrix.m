function [JMatrix,JMatrix_trans] = calcLoveJacobiMatrix(ObsInfo,ModelInfo,m0)
%   Summary of this function goes here.
%   [JMatrix,JMatrix_trans] = calcLoveJacobiMatrix(ObsInfo,ModelInfo,m0)
%   Detailed explanation goes here.
%   The function is for calculating the Jacobi matrix of the misfit
%   function of Love wave dispersion inversion.
%
%
%   IN
%        ObsInfo: the struct about the observed information, it consists of
%                 following fields, 'f' (frequency axis vector), 'pv'
%                 (observed phase velocity values matrix, each column represents
%                 a Love wave mode, 0-elements are the unobserved data point)
%                 'Wd' (it shares the same shape as the observed data pv),
%                 'modeSN' (it is a row vector, it represents the Love
%                 wave mode numerical order).
%      ModelInfo: the struct about the model information during the inversion,
%                 it consists of following fields,'Bound' (searching bound
%                 for each paramter, it's n*2 matrix, n is the number of
%                 reconstructed parameter), 'Ini' (initial model of the inverion)
%                 'den' (density of the inversion, it remains a constant).
%             m0: the current model point, it's a row vector.
%
%
%  OUT
%         JMatrix: the Jacobi matrix, its size is M*paraNum*modeNum, M is
%                  the number of frequency points, paraNum is the number of
%                  inversion parameters, modeNum is the number of rayleigh
%                  wave dispersion mode.
%   JMatrix_trans: the transposition of Jacobi matrix, its size is 
%                  paraNum*M*modeNum.
%
%
%  Author(s): Yan Yingwei
%  Copyright: 2025-2030
%  Revision:  1.0  Date: 8/30/2025
%
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

modeSN = ObsInfo.modeSN;
f = ObsInfo.f;
pvPred0 = ObsInfo.pvPred0;
DEN0 = ModelInfo.den;

paraNum = length(m0);

modeNum = length(modeSN);
M = length(f);
JMatrix = zeros(M,paraNum,modeNum);
JMatrix_trans = zeros(paraNum,M,modeNum);

JMat = zeros(M*paraNum,modeNum);

IndCell = ObsInfo.IndCell;
paraIndSequence = 1:paraNum;
m = [m0,DEN0];
for k=1:modeNum
    Ind_k = IndCell{k};
    lInd_k = length(Ind_k);
    
    iStart_k = Ind_k(1);
    elementNum = lInd_k*paraNum;
    freqSequence = zeros(1,elementNum);
    pvPredSequence = zeros(1,elementNum);
    DValue0Sequence = zeros(1,elementNum);
    jSequence = zeros(1,elementNum);
    
    for i=Ind_k
        freq = f(i);
        pvPred = pvPred0(i,modeSN(k));
        DValue0 = fastlovecalc(pvPred,freq,m);
        freqSequence((i-iStart_k)*paraNum+1:(i-iStart_k+1)*paraNum) = freq;
        pvPredSequence((i-iStart_k)*paraNum+1:(i-iStart_k+1)*paraNum) = pvPred;
        DValue0Sequence((i-iStart_k)*paraNum+1:(i-iStart_k+1)*paraNum) = DValue0;
        jSequence((i-iStart_k)*paraNum+1:(i-iStart_k+1)*paraNum) = paraIndSequence;
    end
    
    parfor iEle=1:elementNum
        freq = freqSequence(iEle);
        pvPred = pvPredSequence(iEle);
        DValue0 = DValue0Sequence(iEle);
        j = jSequence(iEle);
        JMat(iEle,k)= calcLoveJacobiMatrix_Ele(freq,pvPred,DValue0,m0,DEN0,j);
    end
end
for k=1:modeNum
    Ind_k = IndCell{k};
    iStart_k = Ind_k(1);
    for i=Ind_k
        JMatrix(i,:,k) = JMat(1+(i-iStart_k)*paraNum:(i-iStart_k+1)*paraNum,k);
        JMatrix_trans(:,i,k) = JMatrix(i,:,k);
    end
end
end

