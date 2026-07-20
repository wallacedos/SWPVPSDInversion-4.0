 %% 基于预处理最速下降法的瑞雷波多阶频散曲线反演：实测数据测试
 %+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 %
 % 反演中使用的观测数据属于沈白高铁线路基勘探项目中的一个测点，该测点以偏移距为7
 % m，道间距为1 m的24道线性观测系统采集主动源瑞雷波数据，使用F-J变换提取多道面波
 % 记录的频散能量图，然后从频散能量图中拾取得到两条频散曲线，频散曲线的阶次被先验
 % 地设置为基阶和二阶高阶。
 %
 %+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 
 clear;   % 清空matlab的工作区
 
 %% 加载测量的瑞雷波多阶面波相速度
 currDir = cd; filePath = strcat(cd,"\Data\RoadBedSWDC\");
 fileName = strcat(filePath,"rayleighObsData.mat");
 load(fileName);       % 加载观测数据
 [M,N] = size(pvObserved);
 Wd = zeros(M,N);
 for j=1:N
     [Index,~] = find(pvObserved(:,j)~=0);
     Wd(Index,j) = 1;                   % 设置数据加权矩阵，默认为1
 end
 %% 设置初始模型
 numLayer = 10;                                    % 地层设为10层
 [Index,~] = find(pvObserved(:,1)~=0);
 
 lambdaMax = pvObserved(Index(1),1)/f(Index(1));   % 最大波长
 alp = 0.63;
 dh = alp*lambdaMax/(numLayer-1);
 depth = 0:dh:(numLayer-1)*dh+0.001*dh;
 
 % 设置为1利用面波走时构建初始横波速度，设置为0利用面波波长构建初始横波速度
 mId = 1;  
 VSIni = construct1DIni(f(Index),pvObserved(Index,1),depth,alp,mId);
 HIni = depth(2:end)-depth(1:end-1);
 
%  pvBaseFEnd = pvObserved(Index(end),1);
%  pvBaseFStart = pvObserved(Index(1),1);
%  lambdaMax = pvObserved(Index(1),1)/f(Index(1));   % 最大波长
%  
%  VSIni = linspace(pvBaseFEnd/0.88,pvBaseFStart/0.88,numLayer);
%  HIni = 0.5*lambdaMax/(numLayer-1)*ones(1,numLayer-1);

 mIni = [VSIni,HIni];   % 初始模型
 den = 2000*ones(1,numLayer);
 vpdvs = 2.45*ones(1,numLayer); % 反演过程中每层的纵横波速比
 % 横波速度与厚度的搜索范围
 dBound = [50*ones(numLayer,1);0.5*ones(numLayer-1,1)];
 uBound = [2500*ones(numLayer,1); 30*ones(numLayer-1,1)];
 Bound = [dBound, uBound];
 
 % 反演过程中的参考模型，无先验信息时，设为0，只对反演模型做最光滑约束
mReference = zeros(1,2*numLayer-1);
 %% 设置反演参数
 % 预处理最速下降法中的平衡因子，不要取1和0！！！
 theta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1.1 10]; % 比较合理的设置
 % 速度扰动的候选步长，厚度扰动的步长则根据速度扰动步长自适应确定
 alpha = [0.005 0.01 0.02 0.04 0.08];
 
 maxIter = 100;       % 设置最大的迭代次数为80
 misfitType = "L2";   % 目标函数类型选为基于L2范数类型，也可以选择L1型
 stateId = "Fast";    % 以快速模式计算
 
 % 迭代中的反演模型与参考模型最光滑约束的超参数epsilon，通过试触测试确定
epsilonSequence = [1e-3 1e-2 1e-1 1 10];
 %% 设置反演算法需要使用的结构体
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

 tic;
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
 fprintf("The total computing time is %f s.\n",t0);
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
% 也可以人工选择，需要兼顾数据拟合.
[Y,I] = max(abs(curv));
mInv = mInvSequence(I+1,:);
%% 绘制对数域的L曲线
figure;
plot(log(dataMisfit),log(modelMisfit),'k-*');
hold on;
plot(log(dataMisfit(I+1)),log(modelMisfit(I+1)),'r-o');
xlabel('Log of Data misfit term');
ylabel('Log of Model misfit term');
 %% 为绘制频散曲线拟合和反演结果准备数据
 VSInv = mInv(1:numLayer);HInv = mInv(numLayer+1:end);
 maxModeSN = max(modeSN);
 pvPredicted = calcmulti(f,VSInv,HInv,vpdvs.*VSInv,den,maxModeSN);
 
maxDepth = 25;
[VSIniSta,depthIniSta] = calcStairsData(VSIni,HIni, maxDepth);
[VSInvSta,depthInvSta] = calcStairsData(VSInv,HInv, maxDepth);
 
 modeNum = length(modeSN);
 %% 绘制频散曲线拟合和重建的1D横波速度剖面
 figure;
 subplot(1,2,1);
 for j=1:modeNum-1
     [Index,~] = find(pvObserved(:,j)~=0);
     plot(f(Index),pvObserved(Index,j),'r','marker','.','markersize',6,'linestyle','none');
     hold on;
     [Index,~] = find(pvPredicted(:,modeSN(j))~=0);
     plot(f(Index),pvPredicted(Index,modeSN(j)),'k-');
 end
 j = modeNum;
 [Index,~] = find(pvObserved(:,j)~=0);
 plot(f(Index),pvObserved(Index,j),'r','marker','.','markersize',6,'linestyle','none');
 hold on;
 [Index,~] = find(pvPredicted(:,modeSN(j))~=0);
 plot(f(Index),pvPredicted(Index,modeSN(j)),'k-');
 legend('Observed','Predicted');
 xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)');
 title('Comparison of the observed and predicted values');
 subplot(1,2,2);
 h = plot(VSIniSta,depthIniSta,'b-',VSInvSta,depthInvSta,'k-');
 set(gca,'YDir','reverse'); xlim([200 1500]);
 legend(h,'Initial','Inverted');
 xlabel('S-wave velocity (m/s)');ylabel('Depth (m)');
 title('Inverted result');
 
 