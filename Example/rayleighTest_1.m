 %% 基于预处理最速下降法的瑞雷波频散曲线反演：四层含软弱夹层模型测试
 %+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 %            四层含软弱夹层模型参数
 % 层序号    VS     H       VP       den
 %   1      400    4      692.0     2000
 %   2      120    4      207.6     2000
 %   3      500    32     865.0     2000
 %   4      600   半空间  1038.0    2000
 %
 %
 % 反演中使用的观测数据由对四层含软弱夹层模型正演得到，只取基阶，频带为1-100 Hz,
 % 频率间隔为1 Hz.
 % 采用薄层策略，只更新地层横波速度，层厚度为 1 m.
 %+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

clear;   % 清空matlab的工作区
%% 合成四层含软弱夹层模型的观测数据
f = 1:100;                       % 频率，Hz
m = [400 120 500 600 4 4 32];    % 模型，前4个参数为横波速度，后3个参数为层厚度
VS = m(1:4);                     % 模型的横波速度，m/s
H = m(5:7);                      % 模型的层厚度，m
VP = 1.73*VS;                    % 模型的纵波速度，m/s
den = [2000 2000 2000 2000];     % 模型的密度，kg/m^3
pv = calcmulti(f,VS,H,VP,den);   % 正演出的多阶面波相速度，m/s
pvObserved = pv(:,1);            % 只取基阶作反演，因此取第一列
Wd = ones(length(f),1);          % 数据加权矩阵，此处认为所有数据点的权重一样，即都是1

% 观测值的阶次序号向量，1代表只有基阶,[1,2]则是包含基阶和一阶高阶，观测值中的
% 每一列都代表一个阶次的面波相速度，观测值中的零元素代表该数据点处没有观测值
modeSN = [1]; 

% 以1 m薄层离散横波速度
VS_grid = [400*ones(1,4),120*ones(1,4),500*ones(1,32),600];
%% 设置初始模型
% 初始速度模型，由【横波速度向量、厚度向量组成】
% 以简单的半空间模型为初始模型
numLayer = 41;
mIni = [350*ones(1,numLayer),1*ones(1,numLayer-1)];  

% 横波速度与厚度的搜索范围
Bound = zeros(81,2);
Bound(1:4,1) = 300;
Bound(1:4,2) = 600;
Bound(5:8,1) = 50;
Bound(5:8,2) = 600;
Bound(9:41,1) = 50;
Bound(9:41,2) = 1000;

% 反演中层厚度固定为1 m, 即将上界和下界都设置为1 m.
Bound(42:81,1) = 1;   
Bound(42:81,2) = 1;

DEN = 2000*ones(1,numLayer);                 % 反演过程中的密度，kg/m^3
vpdvs = 1.73*ones(1,numLayer);               % 反演过程中每层的纵横波速比

% 反演过程中的参考模型，如果没有先验的参考模型，默认为0，只对反演模型做最光滑约束
mReference = zeros(1,2*numLayer-1);         
%% 设置反演参数
% 预处理最速下降法中的平衡因子，不要取1和0！！！
theta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1.1 10]; % 比较合理的设置
% 按Armijo准则测试出比较合适的百分比扰动步长
alpha = [0.01 0.02 0.04 0.08 0.16 0.32];
maxIter = 50;        % 设置最大的迭代次数为50
misfitType = "L2";   % 目标函数类型选为基于L2范数类型，也可以选择L1型
epsilon = 0.0; % 0.0意味着只拟合数据，不做模型约束，目标函数没有模型惩罚项

stateId = "Fast";    % "Fast"为快速计算模式，"Common"为普通计算模式
%% 设置反演算法需要使用的结构体
ObsInfo.f = f;
ObsInfo.pv = pvObserved;
ObsInfo.Wd = Wd;
ObsInfo.modeSN = modeSN;
ModelInfo.Bound = Bound;
ModelInfo.Ini = mIni;
ModelInfo.den = DEN;
ModelInfo.vpdvs = vpdvs;
ModelInfo.mReference = mReference;
InvParameterInfo.theta = theta;
InvParameterInfo.alpha = alpha;
InvParameterInfo.maxIter = maxIter;
InvParameterInfo.misfitType = misfitType;
InvParameterInfo.epsilon = epsilon;

InvParameterInfo.stateId = stateId; 
%% 反演
% 快速计算模式反演横波速度结构
fprintf("The fast computing mode has begin:\n");
tic;
[mInv_Fast,fValueSequence_Fast,~,~,~] = rayleighDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo);
t0_Fast = toc;

% 普通计算模式反演横波速度结构
stateId = "Common";
InvParameterInfo.stateId = stateId;
fprintf("The common computing mode has begin:\n");
tic;
[mInv_Common,fValueSequence_Common,~,~,~] = rayleighDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo);
t0_Common = toc;
%% 评价反演精度
err_Fast = sum(abs((VS_grid-mInv_Fast(1:numLayer))./VS_grid))/numLayer;
t0_Fast_EachIter = t0_Fast/(length(fValueSequence_Fast)-1);
err_Common = sum(abs((VS_grid-mInv_Common(1:numLayer))./VS_grid))/numLayer;
t0_Common_EachIter = t0_Common/(length(fValueSequence_Common)-1);

fprintf("The total computing time for fast computing mode is %f s, the computing time for each iteration is %f s.\n",t0_Fast,t0_Fast_EachIter);
fprintf("The relative error between the inverted result and true model of fast computing mode is %f.\n",err_Fast);
fprintf("The total computing time for common computing mode is %f s, the computing time for each iteration is %f s.\n",t0_Common,t0_Common_EachIter);
fprintf("The relative error between the inverted result and true model of common computing mode is %f.\n",err_Common);
%% 为绘图准备数据
VSIni = mIni(1:numLayer);HIni = mIni(numLayer+1:end);
VSInv_Fast = mInv_Fast(1:numLayer);
HInv_Fast = mInv_Fast(numLayer+1:end);
pvPredicted_Fast = calcmulti(f,VSInv_Fast,HInv_Fast,vpdvs.*VSInv_Fast,DEN);
VSInv_Common = mInv_Common(1:numLayer);
HInv_Common = mInv_Common(numLayer+1:end);
pvPredicted_Common = calcmulti(f,VSInv_Common,HInv_Common,vpdvs.*VSInv_Common,DEN);

maxDepth = 60;
[VSSta, depthSta] = calcStairsData(VS,H, maxDepth);

[VSIniSta,depthIniSta] = calcStairsData(VSIni,HIni, maxDepth);
[VSInv_FastSta,depthInv_FastSta] = calcStairsData(VSInv_Fast,HInv_Fast, maxDepth);
[VSInv_CommonSta,depthInv_CommonSta] = calcStairsData(VSInv_Common,HInv_Common, maxDepth);
%% 绘制频散曲线拟合和重建的1D横波速度剖面

% 绘制快速计算模式下的频散曲线拟合和反演结果
figure;
subplot(1,2,1);
plot(f,pvObserved,'r','marker','.','markersize',6,'linestyle','none');
hold on;
plot(f,pvPredicted_Fast,'k-');
legend('Observed','Predicted');
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');
title('Comparison of the observed and predicted values');
subplot(1,2,2);
h = plot(VSSta,depthSta,'r-',VSIniSta,depthIniSta,'b-',VSInv_FastSta,depthInv_FastSta,'k-');
set(gca,'YDir','reverse'); xlim([50 650]);
legend(h,'True','Initial','Inverted');
xlabel('S-wave velocity (m/s)');ylabel('Depth (m)');
title('True and inverted result');

% 绘制普通计算模式下的频散曲线拟合和反演结果
figure;
subplot(1,2,1);
plot(f,pvObserved,'r','marker','.','markersize',6,'linestyle','none');
hold on;
plot(f,pvPredicted_Common,'k-');
legend('Observed','Predicted');
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');
title('Comparison of the observed and predicted values');
subplot(1,2,2);
h = plot(VSSta,depthSta,'r-',VSIniSta,depthIniSta,'b-',VSInv_CommonSta,depthInv_CommonSta,'k-');
set(gca,'YDir','reverse'); xlim([50 650]);
legend(h,'True','Initial','Inverted');
xlabel('S-wave velocity (m/s)');ylabel('Depth (m)');
title('True and inverted result');

