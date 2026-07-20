%% 基于预处理最速下降法反演实测的瑞雷波基阶频散曲线
clear;               % 清空matlab工作空间

%% 加载测量的基阶面波相速度
currDir = cd; filePath = strcat(currDir,"\Data\GoafSWDC\");
fileName = strcat(filePath,"L3Merging(10).txt");
obsData = load(fileName); 
f = obsData(:,1); f = f';            % 输入的频率一定是一个行向量!!!
pvObserved = obsData(:,2);           % 基阶瑞雷波相速度，列向量

% 反演中的数据权重矩阵，默认所有数据权重相同，都为 1
Wd = ones(size(pvObserved)); 
modeSN = [1];                   % 阶次标识，1代表只反演基阶
%% 设置初始模型
numLayer = 20;                  % 地层划分为20层
den = 2000*ones(1,numLayer);    % 密度向量
vpdvs = 2.45*ones(1,numLayer);  % 反演过程中每层的纵横波速比
lambdaFundamentalMax = pvObserved(1)/f(1);  % 最大波长

alp = 0.63;

% 建立初始模型和设置模型参数的边界
mIni = [(pvObserved(1)+pvObserved(end))/0.88/2*ones(1,numLayer),0.5*lambdaFundamentalMax/(numLayer-1)*ones(1,numLayer-1)];

vsUBound = max(max(pvObserved))/0.88*2.5;
vsDBound = 50;
Bound = [vsDBound*ones(numLayer,1),vsUBound*ones(numLayer,1);0.005*lambdaFundamentalMax*ones(numLayer-1,1),alp*lambdaFundamentalMax*ones(numLayer-1,1)];

% 反演过程中的参考模型，无先验信息时，设为0，只对反演模型做最光滑约束
mReference = zeros(1,2*numLayer-1);
%% 设置反演参数
% 预处理最速下降法中的平衡因子，不要取1和0！！！
theta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1.1 10]; % 比较合理的设置
% 百分比扰动步长的候选参数，横波速度与厚度同时反演
alpha = [0.005 0.01 0.02 0.04 0.08];
maxIter = 100;       % 最大迭代次数为100
misfitType = "L2";   % 目标函数类型选为基于L2范数类型，也可以选择L1型
stateId = "Fast";    % 计算模式选为快速计算模式

% 迭代中的反演模型与参考模型最光滑约束的超参数epsilon，通过试触测试确定
% 实践上看，从1e-3或者1e-2开始递增比较适合。注意，当设置为0时，就是只拟合数据
epsilonSequence = [1e-3 1e-2 1e-1 1 10 100];
%% 设置用以反演的结构体
ObsInfo.f = f;
ObsInfo.pv = pvObserved;
ObsInfo.Wd = Wd;
ObsInfo.modeSN = modeSN;
ModelInfo.Bound = Bound;
ModelInfo.Ini = mIni;
ModelInfo.den = den;
ModelInfo.vpdvs = vpdvs;
ModelInfo.mReference = mReference;
InvParameterInfo.theta = theta;
InvParameterInfo.alpha = alpha;
InvParameterInfo.maxIter = maxIter;
InvParameterInfo.misfitType = misfitType;
InvParameterInfo.stateId = stateId;
%% 反演
lEpsilon = length(epsilonSequence);

% 不同的epsilon时的反演结果和总目标函数、数据目标函数、模型目标函数的序列
mInvSequence = zeros(lEpsilon,2*numLayer-1);
fValueSequence = zeros(maxIter+1,lEpsilon);
fValueDataSequence = zeros(maxIter+1,lEpsilon);
fValueModelSequence = zeros(maxIter+1,lEpsilon);

tic
for k=1:lEpsilon
    fprintf("k = % d  epsilon = %f \n", k, epsilonSequence(k));
    InvParameterInfo.epsilon = epsilonSequence(k);
    [mInv,fValue,fValueData,fValueModel,~] = rayleighDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo);
    mInvSequence(k,:) = mInv;
    ll = length(fValue);
    fValueSequence(1:ll,k) = fValue;
    fValueDataSequence(1:ll,k) = fValueData;
    fValueModelSequence(1:ll,k) = fValueModel;
end
t0 = toc;
t0_Each = t0/lEpsilon;
fprintf("Total inversion time for all trial tests is % f s.\n", t0);
fprintf("The inversion time for a trial is %f s.\n", t0_Each);
%% L曲线确定最后被选中的反演模型
dataMisfit = zeros(1,lEpsilon);
modelMisfit = zeros(1,lEpsilon);
for k=1:lEpsilon
    Ind = find(fValueDataSequence(:,k)~=0);
    dataMisfit(k) = fValueDataSequence(Ind(end),k);
    modelMisfit(k) = fValueModelSequence(Ind(end),k)/epsilonSequence(k);
end
% 计算对数域的曲率
curv = calcCurvature(log(dataMisfit), log(modelMisfit));

% 曲率绝对值最大的epsilon会被选中，对应的反演结果也会被选中
% 或者人工选择
[Y,I] = max(abs(curv));
I = 3;
mInv = mInvSequence(I+1,:);
%% 绘制对数域的L曲线
figure;
plot(log(dataMisfit),log(modelMisfit),'k-*');
hold on;
plot(log(dataMisfit(I+1)),log(modelMisfit(I+1)),'r-o');
xlabel('Log of Data misfit term');
ylabel('Log of Model misfit term');
%% 为绘制反演结构准备数据
VSIni = mIni(1:numLayer); HIni = mIni(numLayer+1:end);
VSInv = mInv(1:numLayer);HInv = mInv(numLayer+1:end);
pvPredicted = calcmulti(f,VSInv,HInv,vpdvs.*VSInv,den);

maxDepth = 60;
[VSIniSta,depthIniSta] = calcStairsData(VSIni,HIni, maxDepth);
[VSInvSta,depthInvSta] = calcStairsData(VSInv,HInv, maxDepth);
%% 绘制频散曲线拟合和反演的1D横波速度剖面
figure;
subplot(1,2,1);
plot(f,pvObserved,'r','marker','.','markersize',6,'linestyle','none');
hold on;
plot(f,pvPredicted,'k-');
legend('Observed','Predicted');
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');
title('Comparison of the observed and predicted values');
subplot(1,2,2);
h = plot(VSIniSta,depthIniSta,'b-',VSInvSta,depthInvSta,'k-');
set(gca,'YDir','reverse'); xlim([150, 700]);
legend(h,'Initial','Inverted');
xlabel('S-wave velocity (m/s)');ylabel('Depth (m)');
title('Initial and inverted result');