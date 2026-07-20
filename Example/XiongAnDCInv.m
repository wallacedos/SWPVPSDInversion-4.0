%% 基于预处理最速下降法反演雄安测区的瑞雷波多阶频散曲线

clear;     % 清空matlab工作区

currDir = cd;
filePath = strcat(currDir,"\Data\XiongAnSWDC\"); % 频散文件所在路径
df = 0.05;

ptNum = 177;    % 测点数共177
numLayer = 41;  % 地层数设为41

maxIter = 80;        % 最大迭代次数设为80
misfitType = "L2";   % 目标函数类型选为基于L2范数类型，也可以选择L1型
% 预处理最速下降法中的平衡因子，不要取1和0！！！
theta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1.1 10]; % 比较合理的设置
% 模型扰动的百分比候选步长
alpha = [0.005 0.01 0.02 0.04 0.08];

fValueSequence = zeros(2*ptNum,maxIter+1);
fValueDataSequence = zeros(2*ptNum,maxIter+1);
fValueModelSequence = zeros(2*ptNum,maxIter+1);
mInvSequence = zeros(2*ptNum,2*numLayer-1);
t0 = zeros(ptNum,1);

den = 2000*ones(1,numLayer);   % 反演过程中的密度向量
vpdvs = 2.1*ones(1,numLayer);  % 反演过程中每层的纵横波速比
vpdvs(end) = 2.0;

alp = 0.5;                    % 假设经验勘探深度为0.5倍的最大波长

% 反演过程中的参考模型，无先验信息时，设为0，只对反演模型做最光滑约束
mReference = zeros(1,2*numLayer-1);

% 超参数epsilon，被固定为0.2
epsilon = 0.2;

stateId = "Fast"; % 以快速计算模式反演
%% 反演所有测点的频散曲线
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
% 每个测点先作基阶频散曲线反演，再以基阶反演结果为初始模型作多阶频散曲线反演，
% 测点的横向距离为200 m.
%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
for i=1:ptNum
    numStr = num2str(i);
    fileName = strcat(filePath,'chan49Merging(',numStr,').DC');
    fid = fopen(fileName,'r');
    A = textscan(fid,'%s %d\n');
    C = textscan(fid,'%s %f %f %f %d');
    fclose(fid);
    l = length(C{1,1});
    pvRaw = zeros(l,3);
    pvRaw(:,1) = C{1,2}; pvRaw(:,2) = C{1,3}; pvRaw(:,3) = C{1,5};
    modeSN = unique(pvRaw(:,3));
    modeTotalNum = length(modeSN);
    minf = min(pvRaw(:,1));maxf = max(pvRaw(:,1));
    s1 = ceil(minf/df); s2 = floor(maxf/df);
    pv = zeros(s2-s1+1,modeTotalNum);
    f = s1*df:df:s2*df;
    for j=1:modeTotalNum
        Ind = find(pvRaw(:,3)==modeSN(j));
        fq = ceil(pvRaw(Ind(1),1)/df)*df:df:floor(pvRaw(Ind(end),1)/df)*df;
        Vq = interp1(pvRaw(Ind,1),pvRaw(Ind,2),fq,'linear');
        fInd = floor(fq/df+0.1)-s1+1;
        pv(fInd,j) = Vq;
    end
    modeSN = modeSN+1;
    
    Ind = find(pv(:,1)~=0);
    pvObserved = pv(Ind,1);
    lambdaFundamentalMax = pv(Ind(1),1)/f(Ind(1));
    
    
    dh = alp*lambdaFundamentalMax/(numLayer-1);
    depth = 0:dh:(numLayer-1)*dh+0.001*dh;
    
    % 设置为1利用面波走时构建初始横波速度，设置为0利用面波波长构建初始横波速度
    mId = 1;
    VSIni = construct1DIni(f(Ind),pvObserved,depth,alp,mId);
    HIni = depth(2:end)-depth(1:end-1);
    mIni = [VSIni,HIni];
    
%     mIni = [(pvObserved(1)+pvObserved(end))/0.88/2*ones(1,numLayer),0.5*lambdaFundamentalMax/(numLayer-1)*ones(1,numLayer-1)];
   
    vsUBound = max(max(pv))/0.88*2.5;
    vsDBound = 40;
    Bound = [vsDBound*ones(numLayer,1),vsUBound*ones(numLayer,1);0.005*lambdaFundamentalMax*ones(numLayer-1,1),2*lambdaFundamentalMax*ones(numLayer-1,1)];
    Wd = ones(length(Ind),1);
    %% 设置反演算法需要使用的结构体
    ObsInfo.f = f(Ind);
    ObsInfo.pv = pvObserved;
    ObsInfo.Wd = Wd;
    ObsInfo.modeSN = modeSN(1);
    ModelInfo.Bound = Bound;
    ModelInfo.Ini = mIni;
    ModelInfo.den = den;
    ModelInfo.vpdvs = vpdvs;
    ModelInfo.mReference = mReference;
    InvParameterInfo.theta = theta;
    InvParameterInfo.alpha = alpha;
    InvParameterInfo.maxIter = maxIter;
    InvParameterInfo.misfitType = misfitType;
    InvParameterInfo.epsilon = epsilon;
    InvParameterInfo.stateId = stateId;
    
    tic
    % 第一阶段：作基阶频散曲线反演
    fprintf("The %d-th survey point, first stage...\n",i);
    [mInv,fValue,fValueData,fValueModel,~] = rayleighDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo);
    
    mInvSequence(2*i-1,:) = mInv;
    ll = length(fValue); fValueSequence(2*i-1,1:ll) = fValue;
    fValueDataSequence(2*i-1,1:ll) = fValueData;
    fValueModelSequence(2*i-1,1:ll) = fValueModel;
    [M,N] = size(pv);
    Wd = zeros(M,N);
    for j=1:N
        [Index,~] = find(pv(:,j)~=0);
        Wd(Index,j) = 1;                   % 设置数据加权矩阵，默认为1
    end
    
    ObsInfo.f = f;
    ObsInfo.pv = pv;
    ObsInfo.Wd = Wd;
    ObsInfo.modeSN = modeSN;
    ModelInfo.Ini = mInv;
    
    % 第二阶段：作多阶频散曲线反演
    fprintf("The %d-th survey point, second stage...\n",i);
    [mInv,fValue,fValueData,fValueModel,~] = rayleighDCPSDIter(ObsInfo,ModelInfo,InvParameterInfo);
    mInvSequence(2*i,:) = mInv;
    ll = length(fValue); fValueSequence(2*i,1:ll) = fValue;
    fValueDataSequence(2*i-1,1:ll) = fValueData;
    fValueModelSequence(2*i-1,1:ll) = fValueModel;
    t0(i) = toc;
end
%% 对多阶频散曲线作普通插值网格化
maxDepth = 800;
dh = 1;        % 插值后的深度方向的网格长度，m
depthGridNum = floor(maxDepth/dh)+1;
VSTomo = zeros(depthGridNum,ptNum);

DetectDepth = zeros(1,ptNum);
% 先做一维插值
for i=1:ptNum
    DepthInv_i = zeros(1,numLayer);
    VSInv_i = mInvSequence(2*i,1:numLayer);
    HInv_i = mInvSequence(2*i,numLayer+1:end);
    for j=2:numLayer
       DepthInv_i(j) = DepthInv_i(j-1)+HInv_i(j-1);
    end
    [VSInv_i_interp,DepthInv_i_interp] = vsInterp(VSInv_i,DepthInv_i,dh,maxDepth);
    DetectDepth(i) = floor(DepthInv_i(end)/dh)*dh;
    VSTomo(:,i) = VSInv_i_interp;
end
dx = 200;     % 测点间的距离，点距，单位 m
x = 0:dx:(ptNum-1)*dx+0.001*dx;
[xGrid,depthGrid] = meshgrid(x,DepthInv_i_interp);

dx_i = 5;    % 插值后的横向网格长度，m
x_i = 0:dx_i:(ptNum-1)*dx+0.001*dx_i;
[xGrid_i,depthGrid_i] = meshgrid(x_i,DepthInv_i_interp);

VSTomo_i = interp2(xGrid,depthGrid,VSTomo,xGrid_i,depthGrid_i,"cubic");
VSTomo_i_g = imgaussfilt(VSTomo_i,50);  % 高斯平滑后的横波速度

% 白化方法为后向格式白化
IndFirst = 1;
for i=2:ptNum
    Ind = find(x_i<x(i));
    Ind_Depth = floor(DetectDepth(i)/dh+1e-6)+1;
    VSTomo_i(Ind_Depth+1:end,IndFirst:Ind(end)) = nan;
    VSTomo_i_g(Ind_Depth+1:end,IndFirst:Ind(end)) = nan;
    IndFirst = Ind(end)+1;
end
VSTomo_i(Ind_Depth+1:end,end) = nan;
VSTomo_i_g(Ind_Depth+1:end,end) = nan;

% 绘制网格化后的横波速度拟2D剖面
figure;
pcolor(xGrid_i./1000,depthGrid_i./1000,VSTomo_i_g);shading interp;
set(gca,'Ydir','reverse');c = colorbar;
xlabel(c,'S-wave velocity (m/s)');
xlabel('Distance (km)');ylabel('Depth (km)');
title('Inversion result of multimodal surface wave');
colormap(jet); caxis([300 1000]);
pbaspect([5,1,1]);


% %% 对多阶频散曲线反演的结果作克里金网格化
% % 上方普通插值的代码块注释后，可执行克里金网格化，计算时间较长
% X = zeros(ptNum*numLayer,1); % 横向位置
% Y = zeros(ptNum*numLayer,1); % 垂直位置，深度
% Z = zeros(ptNum*numLayer,1); % 反演的横波速度
% 
% dx = 200; % 测点间的距离，点距，单位 m
% for i=1:ptNum
%     for j=1:numLayer
%         X((i-1)*numLayer+j) = (i-1)*200;
%         Y((i-1)*numLayer+j) = sum(mInvSequence(2*i,numLayer+1:numLayer+j-1));
%         Z((i-1)*numLayer+j) = mInvSequence(2*i,j);
%     end
% end
% 
% [VSTomo,gridX,gridY] = kriging(X,Y,Z);
% 
% % 绘制网格化后的横波速度拟2D剖面
% figure;
% pcolor(gridX,gridY,VSTomo');shading interp;
% set(gca,'Ydir','reverse');colorbar;
% xlabel('Position (m)');ylabel('Depth (m)');
% title('Inversion result of multimodal surface wave');

