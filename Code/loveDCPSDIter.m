function [mInv,fValueSequence,fValueSequenceData,fValueSequenceModel,mAccepted] = loveDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo)
%   Summary of this function goes here.
%  [mInv,fValueSequence,fValueSequenceData,fValueSequenceModel,mAccepted] = loveDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo)
%   Detailed explanation goes here.
%   The function is for inverting the fundamental/multimodal Love-wave
%   dispersion curves by preconditioned steepest-descent method.
%
%   IN
%        ObsInfo: the struct about the observed information, it consists of
%                 following fields, 'f' (frequency axis row vector), 'pv'
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
%                 'mReference' is the predifined reference model, default as 0.
%                  full 0 vector is most smooth
%                 constraint,otherwise is the minimal L2 norm constraint
%                 between the model and reference model.
%InvParameterInfo: it's a struct about the inversion parameter setting, it
%                consists of following fields, 'theta' (balance operator vector
%                about the Hessian operator, it is a row vector), 'alpha',
%                (it is trial step-length vector for inversion, it reflects
%                the percent perturbation about the current model, step-length),
%                'maxIter' (the maximum iterations of inversion),'misfitType'
%                (the choice of type for misfit function, 'L1' and 'L2' are
%                the alternative choice),'epsilon' is a hyper-parameter for model
%                punishment, the most smooth constraint of the difference
%                of the reference model and the model in iteration, "stateId"
%                represents the computing mode settings, there are two choices,
%                "Fast" means the fast computing mode, "Common" means the 
%                common computing mode.
%
%
%  OUT
%                mInv: the inverted result after iterations.
%      fValueSequence: the misfit function values of the iterations.
%  fValueSequenceData: the misfit function values of data fitting term.
% fValueSequenceModel: the misfit function values of the model punishment
%                      term.
%           mAccepted: the accepted models of each iteration.
%
%
% EXAMPLE
%
%  f = 5:100;
%  m = [150 450 600 4 8];den = [2000 2000 2000];pv = calclovemulti(f,m(1:3),m(4:5),den);
%  pvObserved = pv(:,1);
%  Wd = ones(length(f),1);
%
%  ObsInfo.f = f;
%  ObsInfo.pv = pvObserved;
%  ObsInfo.Wd = Wd;
%  ObsInfo.modeSN = [1];
%  ModelInfo.Bound = [50 500;50 1000;50 1000;0.5 30;0.5 30];
%  ModelInfo.Ini = [400 400 410 6 6];
%  ModelInfo.den = den;
%  InvParameterInfo.theta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1.1 10];
%  InvParameterInfo.alpha = [0.01 0.02 0.04 0.08 0.16 0.32];
%  InvParameterInfo.maxIter = 30;
%  InvParameterInfo.misfitType = 'L2';
%  InvParameterInfo.stateId = "Fast";
%  ModelInfo.Reference = zeros(1,length(m));
%  InvParameterInfo.epsilon = 0.0; % 0.0 means only fitting the observed data
%  global n_mode;n_mode = max(ObsInfo.modeSN);
%  [mInv,fValueSequence,fValueSequenceData,fValueSequenceModel,mAccepted] = loveDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo);
%
%
%  References:
%  Yan, Y., Wang, Z., Li, J., Huai, N., Liang, Y., Song, S., Zhang, J.,
%  & Zhang, L. 2020. Elastic SH- and Love-wave Full-Waveform Inversion for
%  shallow shear wave velocity with a preconditioned technique, Journal of
%  Applied Geophysics, 173, 103947,
%  https://doi.org/10.1016/j.jappgeo.2020.103947.
%
%  Yan, Y., Chen, X., Huai, N., Guan, J.2022.Modern inversion workflow of
%  the multimodal surface wave dispersion curves: Staging strategy and Pattern
%  search with embedded Kuhn-Munkres algorithm, Geophysical Journal
%  International,231(01), 47-71,
%  https://doi.org/10.1093/gji/ggac178.
%
%
%  You can find the principles of our method in any book on optimization
%  methods.
%
%  History:
%          1) Adopting the percent step-length for updating the model
%          parameters at 2/10/2022.
%          2) Adding the most smooth model difference £¨the difference of
%          reference model and the model in iteration) punishment,
%          11/28/2022.
%          3£© Adding the fast computing method, the computational 
%          performance of the inversion program has been significantly 
%          improved, this enhancement is achieved through the application 
%          of the implicit function theorem, 8/30/2025
%
%  Author(s): Yan Yingwei
%  Copyright: 2021-2030
%  Revision:  4.0  Date: 8/30/2025
%
%  College of Geoexploration Science and Technology, Jilin University.
%
%  Department of Earth and Space Sciences, Southern University of Science
%  and Technology (SUSTech).

%% input parameter verification, when the input parameter is 2
if nargin==2
    InvParameterInfo.theta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1.1 10];
    InvParameterInfo.alpha = [0.01 0.02 0.04 0.08 0.16 0.32];
    InvParameterInfo.maxIter = 50;
    InvParameterInfo.misfitType = "L2";
    
    % 0.0 means only fitting the observed data
    InvParameterInfo.epsilon = 0.0;
    
    % "Fast" means the fast computing mode, "common" means the common
    % computing mode
    InvParameterInfo.stateId = "Fast";
end
%% read the input parameter
freq = ObsInfo.f;         % frequency axis vector (Hz), a row vector
pvObserved = ObsInfo.pv;  % the observed phase velocity values matrix
[M,N] = size(pvObserved); % the shape of variable 'pvObserved'

modeSNState = isfield(ObsInfo,"modeSN");
if modeSNState
    modeSN= ObsInfo.modeSN; % the numerical order of Rayleigh wave mode
else
    modeSN = 1:N;
    ObsInfo.modeSN = modeSN;
end

WdInputState = isfield(ObsInfo,"Wd");
if WdInputState
    Wd = ObsInfo.Wd;   % the weighted matrix of observed data,size is [M,N]
else
    Wd = ones(M,N);
    ObsInfo.Wd = Wd;
end

maxModeSN = max(modeSN);

ModelBoundConstraint = ModelInfo.Bound; % searching bound of parameters
mCurrIni = ModelInfo.Ini; % intial model, as [S-wave velocity vector, thickness vector]

denInputState = isfield(ModelInfo,"den");
if denInputState
    den = ModelInfo.den; % constant density vector during the inversion
else
    nn = floor((length(mCurrIni)+1+0.1)/2);
    den = 2000*ones(1,nn);
    ModelInfo.den = den;
end

mReferenceState = isfield(ModelInfo,"mReference");

if mReferenceState
    mReference = ModelInfo.mReference;
else
    nn = length(mCurrIni);
    mReference = zeros(1,nn);
end


theta = InvParameterInfo.theta;   % balance factor of Hessian matrix
alpha = InvParameterInfo.alpha;  % trial S-wave velocity step-length
maxIter = InvParameterInfo.maxIter; % maximum iterations
misfitType = InvParameterInfo.misfitType; % misfit function type,'L1' or 'L2'

% equals to 0.0 means only fitting the observed data
epsilon = InvParameterInfo.epsilon;

stateId = InvParameterInfo.stateId;

IndCell = cell(N,1);
for j=1:N
    Ind = find(pvObserved(:,j)~=0);
    IndCell{j} = Ind';
end
ObsInfo.IndCell = IndCell;

WdSquare = zeros(N,1);
for j=1:N
    WdSquare(j) = sum(Wd(IndCell{j},j).*Wd(IndCell{j},j),"all");
end
ObsInfo.WdSquare_reci = 1./sum(WdSquare,"all");
%% inversion by iterations
lTheta = length(theta);
lAlpha = length(alpha);
n = length(mCurrIni);
mAccepted = zeros(maxIter,n);
mCurrSearch = zeros(lTheta*lAlpha,n);
fCurrSearchValue = zeros(lTheta*lAlpha,1);
numLayer = floor((n+1+0.2)/2);   % the number of layers

fValueSequence = zeros(maxIter+1,1);
fValueSequenceData = zeros(maxIter+1,1);
fValueSequenceModel = zeros(maxIter+1,1);
k = 0;
VSIni = mCurrIni(1:numLayer);
HIni = mCurrIni(numLayer+1:end);
pvPredictedCal = calclovemulti(freq,VSIni,HIni,den,maxModeSN);


pvPredicted = zeros(M,maxModeSN);
[M1,N1] = size(pvPredictedCal);

if (N1<maxModeSN)
    pvPredicted(1:M1,1:N1) = pvPredictedCal;
else
    pvPredicted(1:M1,1:maxModeSN) = pvPredictedCal(:,1:maxModeSN);
end
% mDifference = mCurrIni-mReference;
fCurrValue = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,epsilon, mCurrIni,mReference,misfitType);
fValueSequence(k+1) = fCurrValue;
fValueSequenceModel(k+1) = calcModelPenaltyTerm(epsilon,mCurrIni,mReference);
fValueSequenceData(k+1) = fValueSequence(k+1)-fValueSequenceModel(k+1);

fprintf("The %d iteration: current objective function value, data misfit, model punish term,are %f, %f, %f.\n",k,fCurrValue,fValueSequenceData(k+1),fValueSequenceModel(k+1));

balanceFactor = 0.0;
trialNum = lTheta*lAlpha;
pvPredictedCurr = pvPredicted;

pvPredictedMatrix = zeros(M,maxModeSN,trialNum);

thetaTrial = zeros(1,trialNum);
alphaTrial = zeros(1,trialNum);
for iTrial=1:trialNum
    iTheta = floor((iTrial-0.5)/lAlpha)+1;
    jAlpha = iTrial-(iTheta-1)*lAlpha;
    thetaTrial(iTrial) = theta(iTheta);
    alphaTrial(iTrial) = alpha(jAlpha);
end

for k=1:maxIter
    % calculate the gradient and the diagonal elements of Hessian matrix at
    % current model point
    ObsInfo.pvPred0 = pvPredictedCurr;
    pvPredictedMatrix(:,:,:) = 0;
    [gk,~,HessDiagElementK] = calcLoveDCInvGradHess(ObsInfo,ModelInfo,mCurrIni,fCurrValue,epsilon,mReference,misfitType,stateId);
    
    % the maximum absolute value of the diagonal elements of Hessian matrix
    maxAbsHessDiagElement = max(abs(HessDiagElementK));
    fCurrSearchValue(:,:) = 0;
    currMinVS = min(mCurrIni(1:numLayer));
    currMinH = min(mCurrIni(numLayer+1:end));
    
    % loop for trial models
    parfor iTrial=1:trialNum
        dk = -gk./(HessDiagElementK+thetaTrial(iTrial)*maxAbsHessDiagElement);
        dk(1:numLayer) = dk(1:numLayer)./(max(abs(dk(1:numLayer)))+balanceFactor);
        dk(numLayer+1:end) = dk(numLayer+1:end)./(max(abs(dk(numLayer+1:end)))+balanceFactor);
        dk(isnan(dk))=0;
        dk = dk';
        mTemp = mCurrIni;
        mTemp(1:numLayer) = mTemp(1:numLayer)+alphaTrial(iTrial)*currMinVS*dk(1:numLayer);
        mTemp(numLayer+1:end) = mTemp(numLayer+1:end)+alphaTrial(iTrial)*currMinH*dk(numLayer+1:end);
        mTemp = rebound(mTemp,ModelBoundConstraint);
        mCurrSearch(iTrial,:) = mTemp;
        
        VSTemp = mTemp(1:numLayer);HTemp = mTemp(numLayer+1:end);
        pvPredicted = calclovemulti(freq,VSTemp,HTemp,den,maxModeSN);
        pvPredicted = disDataRegular(pvPredicted,maxModeSN);
        pvPredictedMatrix(:,:,iTrial) = pvPredicted;
        fCurrSearchValue(iTrial) = calcMisfitOfDCInv(pvObserved,pvPredicted,modeSN,Wd,epsilon, mTemp,mReference,misfitType);
    end
    
    % choose the optimal model of the current iteration
    [~,I] = min(fCurrSearchValue);
    if fCurrSearchValue(I)<fCurrValue  % current iteration is successful
        mCurrIni = mCurrSearch(I,:);
        fCurrValue = fCurrSearchValue(I);
        fValueSequence(k+1) = fCurrValue;
        fValueSequence(k+1) = fCurrValue;
        fValueSequenceModel(k+1) = calcModelPenaltyTerm(epsilon,mCurrIni,mReference);
        fValueSequenceData(k+1) = fValueSequence(k+1)-fValueSequenceModel(k+1);
        mAccepted(k,:) = mCurrIni;
        pvPredictedCurr = pvPredictedMatrix(:,:,I);
    else  % current iteration has failed, the inversion is terminated, and the program returns.
        mInv = mCurrIni;
        fValueSequence = fValueSequence(1:k);
        fValueSequenceData = fValueSequenceData(1:k);
        fValueSequenceModel = fValueSequenceModel(1:k);
        mAccepted = mAccepted(1:k-1,:);
        return;
    end
    fprintf("The %d iteration: current objective function value, data misfit, model punish term,are %f, %f, %f.\n",k,fCurrValue,fValueSequenceData(k+1),fValueSequenceModel(k+1));
end
mInv = mCurrIni;
end

