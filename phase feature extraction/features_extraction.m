close all,clear all,clc;

%% 细胞相位特征提取
folder_path = '/Users/taylor/Downloads/Data/cell_data/Breast_cancer_subtype/HER2_SKBR3/phase_example'; %修改文件夹名称
files_list = dir(fullfile(folder_path,'*.mat'));
num_files = numel(files_list); % 获取文件夹中符合要求的文件个数
for t =1:num_files
    filename = files_list(t).name; % 获取文件名
    filepath = fullfile(files_list(t).folder, filename); % 获取文件路径
%%%%%%%%%%%%%%%%开始编写处理图片的代码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %细胞图像读取
    I=importdata(filename);
    se=strel('disk',5);
    I=imclose(I,se); %%闭合（先膨胀后腐蚀），填充物体内细小空洞，连接邻近物体，平滑边界
    I=imfill(I); %%填充灰度图像中的由边界包围的空洞区域
      %figure,imshow(I);
      %colormap(parula);
    [CM]=maxLianTongYu(I); %寻找最大连通域；Cell mask function
      %figure,imshow(CM);colormap(parula);
    %se=strel('disk',5);
    %CM=imclose(CM,se); %%闭合（先膨胀后腐蚀），填充物体内细小空洞，连接邻近物体，平滑边界
    %CM=imfill(CM,'holes'); %%填充二值图像中的由边界包围的空洞区域
    cell_unwrapped_phase=CM.*I; %最大连通域作为细胞二值掩膜
      %figure,imshow(cell_unwrapped_phase);colormap(parula);
      %figure,mesh(cell_unwrapped_phase);colormap(parula);

    %细胞特征提取
    % (1)体特征-Bulk features: cell size, cell mass, and the cell shape
    % 0)基本参数设定
    Np=length(find(CM==1)); %Np: 细胞(ROI)所占像素点数;the number of pixels included in the ROI
    Lp=2.5; %Lp: 像元尺寸，单位um
    M=40; %M：像面全息放大倍数，为40倍
    lambda=0.6328; %光波长，单位um
    [B,L]=bwboundaries(CM,'noholes'); %L是对连通区域进行标号的结果
    C=regionprops(L,'centroid'); %计算每个细胞（白色连通域）的中心点坐标
      %imshow(CM);
      %hold on
    centroids=cat(1, C.Centroid); %质心横纵坐标
    xc=centroids(1,1); %xc存放的是第i个细胞质心的横坐标
    yc=centroids(1,2); %yc存放的是第i个细胞质心的纵坐标
      %plot(xc, yc, '*', 'Color', 'r');
      %hold on
    % 1）形状测量值
    ShapeArea=(Lp^2*Np)/(M^2); %Area:细胞投影面积，单位um^2
    Per=regionprops(L,'Perimeter'); 
    ShapePerimeter=cat(1, Per.Perimeter)*Lp/M; %P:细胞周长，单位um
    Lmax=regionprops(L,"MajorAxisLength");
    ShapeLmajor=cat(1, Lmax.MajorAxisLength)*Lp/M; %Lmajor:长轴长度，单位um
    Lmin=regionprops(L,"MinorAxisLength");
    ShapeLminor=cat(1, Lmin.MinorAxisLength)*Lp/M; %Lminor:短轴长度，单位um
    ShapeElogation=ShapeLminor/ShapeLmajor; %Aspect_Ratio:细胞纵横比，单位1
    %Volume=4/3*pi*(Lminor/2)^2*(Lmajor/2); %volume: 细胞像素体积，单位um^3
    %Circularity=4*pi*Area/(P^2); %Circularity：细胞圆度，单位1
    Cir=regionprops(L,'Circularity'); 
    ShapeCircularity=cat(1, Cir.Circularity); %Circularity：细胞圆度，单位1
    Ecc=regionprops(L,'Eccentricity'); 
    ShapeEccentricity=cat(1, Ecc.Eccentricity); %Eccentricity：细胞偏心率，单位1
    %Ori=regionprops(L,'Orientation'); 
    %Orientation=cat(1, Ori.Orientation); %Orientation：x轴与长轴（该椭圆与区域具有相同的二阶矩）之间的角度，单位°
    x=B{1}(:,2); %x存储的是第i个细胞边缘所有点的横坐标
    y=B{1}(:,1); %y存储的是第i个细胞边缘所有点的纵坐标
    distance=sqrt((x-xc).^2 + (y-yc).^2)*Lp/M; %distance中存储的是第i个细胞边缘所有点到中心点的距离，单位um
    distmin=min(distance); %distmin中存储的是变量distance中最短的距离（细胞最短径），单位um
    distmax=max(distance); %distmax中存储的是变量distance中最长的距离（细胞最长径），单位um
    dislen=length(distance);
        for j=1:dislen
            if j==1
                smoo(j)=abs(distance(1)-((distance(2)+distance(dislen))/2));
            elseif j==dislen
                smoo(j)=abs(distance(dislen)-((distance(1)+distance(dislen-1))/2));
            else
                smoo(j)=abs(distance(j)-((distance(j-1)+distance(j+1))/2));
            end
        end
    ShapeSmoothness=sum(smoo); %Smoothness：细胞光滑度，单位um
    % 2）像素值测量值
    Phase=cell_unwrapped_phase*lambda/(2*pi); %以OPD光延迟图像表示"相位"，单位um
      %figure,imshow(Phase);
      %colormap(parula);
    %Phase_Mean=mean(mean(Phase));%细胞相位均值，单位um
    %Phase_Max=max(max(Phase));%细胞相位最大值，单位um
    alpha=0.2; %折射增量，近似为0.18-0.21mL/g
    BoundingBox=regionprops(L,"BoundingBox");
    BoundingBox=cat(1, BoundingBox.BoundingBox);
    ShapeVolume=0; %Volume:细胞相位体积，单位um^3;以OPD表示相位因此并非为细胞实际体积
    for i=ceil(BoundingBox(1,1)):(floor(BoundingBox(1,1))+BoundingBox(1,3))
        for j=ceil(BoundingBox(1,2)):(floor(BoundingBox(1,2))+BoundingBox(1,4))
            ShapeVolume=ShapeVolume+Phase(i,j)*Lp^2/(M^2); %细胞相位体积，单位um^3（fL）
        end
    end
    [Phase_gradientx,Phase_gradienty]=gradient(Phase);
    ShapeSurface=0+ShapeArea; %细胞相位表面积，单位um^2;以OPD表示相位因此并非为细胞实际表面积
    for i=ceil(BoundingBox(1,1)):(floor(BoundingBox(1,1))+BoundingBox(1,3))
        for j=ceil(BoundingBox(1,2)):(floor(BoundingBox(1,2))+BoundingBox(1,4))
            ShapeSurface=ShapeSurface+(1+Phase_gradientx(i,j)^2+Phase_gradienty(i,j)^2)^(1/2)*Lp^2/(M^2); %细胞相位表面积，单位um^2
        end
    end
    ShapeDryMass=ShapeVolume/alpha; %Dry Mass: 细胞干重，单位pg
    Phasemap=imcrop(Phase,BoundingBox); %用最小邻接矩阵裁剪图像，使特征提取更精准，单位um
      %figure,imshow(Phasemap);colormap(parula);
    [m,n]=size(Phasemap);
    [B1,L1]=bwboundaries(Phasemap,'noholes'); %L是对连通区域进行标号的结果
    C1=regionprops(L1,'centroid'); %计算每个细胞（白色连通域）的中心点坐标
    centroids1=cat(1, C1.Centroid); %质心横纵坐标
    xcen1=centroids1(1,1); %细胞相位图像质心的横坐标
    ycen1=centroids1(1,2); %细胞相位图像质心的纵坐标
    Cell_height_map=zeros(m,n);
    for i=1:m
        for j=1:n
            Cell_height_map(i,j)=((ShapeLminor+ShapeLmajor)^2/4-((i-xcen1)^2+(j-ycen1)^2)*(Lp^2)/(M)^2)^(1/2); %细胞高度，单位um
        end
    end
      %figure,mesh(Cell_height_map);colormap(parula);
    Dry_mass_density_map=(Phasemap./Cell_height_map)/alpha; %干重密度分布图像，单位g/mL
      %figure,imshow(Dry_mass_density_map);colormap(parula);
    ShapeDryMassDensity=mean2(Dry_mass_density_map);%Dry Mass Density:干重密度，单位g/mL;以OPD表示相位因此并非为细胞实际干重密度
    ShapeSurface_to_Volume_Ratio=ShapeSurface/ShapeVolume; %细胞相位体积与表面积比，单位um-1
    %Phase_Surface_to_Dry_Mass=DM/S; %细胞干重与表面积比，单位g/m^2
    %Phase_Area_to_Volume=V/Area; %细胞体积与投影面积比，单位um
    ShapeSphericity=(36*pi*ShapeVolume^2)^(1/3)/ShapeSurface; %Sphericity：细胞相位球度，单位1
    
    % (2)一阶统计量/直方图特征/灰度特征-First Order Features
    %CMmap=imcrop(CM,BoundingBox); %用最小邻接矩阵裁剪初始二值掩膜图像
      %figure,imshow(CMmap);colormap(parula);
    %[B0,L0]=bwboundaries(CMmap,'noholes'); %L是对连通区域进行标号的结果
    %C0=regionprops(L0,'centroid'); %计算每个细胞（白色连通域）的中心点坐标
    %centroids0=cat(1, C0.Centroid); %质心横纵坐标
    %xcen0=centroids0(1,1); %细胞二值掩膜质心的横坐标
    %ycen0=centroids0(1,2); %细胞二值掩膜质心的纵坐标
    
    %Phasemap的特征
    Phasemap=imcrop(Phase,BoundingBox); %用最小邻接矩阵裁剪图像，使特征提取更精准，单位um
      %figure,imshow(Phasemap);colormap(parula);
    PhaseCM=imcrop(CM,BoundingBox); 
      %figure,imshow(PhaseCM);colormap(parula);
    pixelnum=sum(sum(PhaseCM==1));%该变量存放ROI内的像素数
    PhasemapROI=Phasemap(PhaseCM==1);%该变量存放对应于ROI区域的原始图像中的灰度值
    [count,level]=imhist(PhasemapROI);%count变量存放每个像素值的数量，level变量存放图像的灰度级
      %figure,imhist(PhasemapROI);
    FoFMean=mean(PhasemapROI);%相位均值
    FoFMedian=median(PhasemapROI);%相位中值
    FoFRange=range(PhasemapROI,"all"); %相位极差
    FoFMinimum=min(PhasemapROI);%相位最小值
    FoFMaximum=max(PhasemapROI);%相位最大值
    FoFVariation=var(PhasemapROI,0,'all');%相位方差
    FoFSkewness=skewness(PhasemapROI,0,'all'); %相位偏度
    FoFKurtosis=kurtosis(PhasemapROI,0,'all'); %相位峰度
    FoFEntropy=entropy(PhasemapROI); %相位熵
    possibility=count/pixelnum; %该变量存放每个灰度值的概率
    FoFUniformity=sum(possibility.^2);%相位均匀度
    FoFPercentile10=prctile(PhasemapROI,10);%相位10分位点
    FoFPercentile90=prctile(PhasemapROI,90);%相位90分位点
    FoFInterquartileRange=prctile(PhasemapROI,75)-prctile(PhasemapROI,25);%相位四分位距离
    minu=PhasemapROI-sum(PhasemapROI)/pixelnum;%该变量存放各像素值的离差
    FoFMAD=sum(abs(minu))/pixelnum;%相位平均绝对误差 Mean Absolute Deviation
    PhasemapROI1090=PhasemapROI(FoFPercentile10<PhasemapROI<FoFPercentile90);
    pixelnum1090=length(PhasemapROI1090);
    rminu=PhasemapROI1090-sum(PhasemapROI1090)/pixelnum1090;%该变量存放各像素值的离差
    FoFrMAD=sum(abs(rminu))/pixelnum1090;%相位鲁棒平均绝对误差 robust Mean Absolute Deviation
    FoFRMS=sqrt(sum(PhasemapROI.^2)/pixelnum); %相位均方根 Root Mean Squared
    FoFEnergy=sum(PhasemapROI.^2); %相位能量
    %[B1,L1]=bwboundaries(Phasemap,'noholes'); %L是对连通区域进行标号的结果
    %C1=regionprops(L1,'centroid'); %计算每个细胞（白色连通域）的中心点坐标
    %centroids1=cat(1, C1.Centroid); %质心横纵坐标
    %xcen1=centroids1(1,1); %细胞相位图像质心的横坐标
    %ycen1=centroids1(1,2); %细胞相位图像质心的纵坐标
    %Phase_Centroid_Displacement=((xcen1-xcen0)^2+(ycen1-ycen0)^2)^(1/2)*Lp/M; %细胞相位质心偏移距离，单位um
    
    % (3)纹理特征-Texture Features
    % 1)灰度共生矩阵 (Grey Level Co-occurrence Matrix, GLCM）
    %获取GLCM
    offsets=[0 1; -1 1;-1 0;-1 -1];
    [glcms,SI]=graycomatrix(Phasemap,'NumLevels',16,'G',[],'Offset',offsets);
      %imshow(rescale(SI));
    %计算灰度共生矩阵四个方向的特征值均值
    [out1]=GLCM_Features1(glcms(:,:,1),0);
    [out2]=GLCM_Features1(glcms(:,:,2),0);
    [out3]=GLCM_Features1(glcms(:,:,3),0);
    [out4]=GLCM_Features1(glcms(:,:,4),0);
    GLCMAutoc=(out1.autoc+out2.autoc+out3.autoc+out4.autoc)/4; % Autocorrelation
    GLCMContr=(out1.contr+out2.contr+out3.contr+out4.contr)/4; % Contrast
    GLCMCorre=(out1.corrm+out2.corrm+out3.corrm+out4.corrm)/4; % Correlation 
    GLCMCprom=(out1.cprom+out2.cprom+out3.cprom+out4.cprom)/4; % Cluster Prominence
    GLCMCshad=(out1.cshad+out2.cshad+out3.cshad+out4.cshad)/4; % Cluster Shade
    GLCMDissi=(out1.dissi+out2.dissi+out3.dissi+out4.dissi)/4; % Dissimilarity=Difference Average
    GLCMEnerg=(out1.energ+out2.energ+out3.energ+out4.energ)/4; % Energy
    GLCMEntro=(out1.entro+out2.entro+out3.entro+out4.entro)/4; % Entropy
    GLCMHomom=(out1.homom+out2.homom+out3.homom+out4.homom)/4; % Homogeneity=Inverse difference (INV)
    GLCMMaxpr=(out1.maxpr+out2.maxpr+out3.maxpr+out4.maxpr)/4; % Maximum probability
    GLCMSosvh=(out1.sosvh+out2.sosvh+out3.sosvh+out4.sosvh)/4; % Sum of sqaures: Variance
    GLCMSavgh=(out1.savgh+out2.savgh+out3.savgh+out4.savgh)/4; % Sum average
    GLCMSvarh=(out1.svarh+out2.svarh+out3.svarh+out4.svarh)/4; % Sum variance
    GLCMSenth=(out1.senth+out2.senth+out3.senth+out4.senth)/4; % Sum entropy
    GLCMDvarh=(out1.dvarh+out2.dvarh+out3.dvarh+out4.dvarh)/4; % Difference variance 
    GLCMDenth=(out1.denth+out2.denth+out3.denth+out4.denth)/4; % Difference entropy
    GLCMInf1h=(out1.inf1h+out2.inf1h+out3.inf1h+out4.inf1h)/4; % Information measure of correlation1
    GLCMInf2h=(out1.inf2h+out2.inf2h+out3.inf2h+out4.inf2h)/4; % Informaiton measure of correlation2 
    GLCMIndnc=(out1.indnc+out2.indnc+out3.indnc+out4.indnc)/4; % Inverse difference normalized (INN)
    GLCMIdmnc=(out1.idmnc+out2.idmnc+out3.idmnc+out4.idmnc)/4; % Inverse difference moment normalized
    %stats=graycoprops(glcms,'all');
    %Con=[stats.Contrast];%对比度
    %H=[stats.Homogeneity];%同质性
    %Cor=[stats.Correlation];%相关性
    %En=[stats.Energy];%能量
    
    % 2)灰度梯度共生矩阵 (Gray Gradient Co-occurrence Matrix，GGCM）
    [GGCMout]=GrayGradinet(Phasemap);
    GGCMsmgrd=GGCMout.smgrd; %Small Grads Dominance
    GGCMbigrd=GGCMout.bigrd; %Big Grads Dominance
    GGCMgrasy=GGCMout.grasy; %Gray Asymmetry
    GGCMgdasy=GGCMout.gdasy; %Grads Asymmetry
    GGCMenerg=GGCMout.energ; %Energy
    GGCMgrmea=GGCMout.grmea; %Gray Mean
    GGCMgdmea=GGCMout.gdmea; %Grads Mean
    GGCMgrvar=GGCMout.grvar; %Gray Variance
    GGCMgdvar=GGCMout.gdvar; %Grads Variance
    GGCMcorre=GGCMout.corre; %Correlation 
    GGCMgrent=GGCMout.grent; %Gray Entropy
    GGCMgdent=GGCMout.gdent; %Grads Entropy
    GGCMentro=GGCMout.entro; %Entropy
    GGCMinert=GGCMout.inert; %Inertia
    GGCMhomog=GGCMout.homog; %Homogeneity
    
    % 3)灰度游程矩阵（Gray Level Run Length Matrix，GLRLM)
    [GLRLMout]=GLRLM_Features(Phasemap);
    GLRLMsre=GLRLMout.sre; % short run emphasis
    GLRLMlre=GLRLMout.lre; % long run emphasis
    GLRLMgln=GLRLMout.gln; % gray level non-uniformity
    GLRLMglnn=GLRLMout.glnn; % Gray Level Non-Uniformity Normalized (GLNN)
    GLRLMrp=GLRLMout.rp; % run percentage
    GLRLMrln=GLRLMout.rln; % run length non-uniformity
    GLRLMrlnn=GLRLMout.rlnn; % run length non-uniformity normalized
    GLRLMrv=GLRLMout.rv; % run variance
    GLRLMglv=GLRLMout.glv; % run variance
    GLRLMlgre=GLRLMout.lgre; % low gray level run emphasis
    GLRLMhgre=GLRLMout.hgre; % high gray level run emphasis
    GLRLMsrlgle=GLRLMout.srlgle; % run variance
    GLRLMsrhgle=GLRLMout.srhgle; % run variance
    GLRLMlrlgle=GLRLMout.lrlgle; % low gray level run emphasis
    GLRLMlrhgle=GLRLMout.lrhgle; % high gray level run emphasis

    
  filename=['/Users/taylor/Downloads/Data/cell_data/Breast_cancer_subtype/HER2_SKBR3/feature_example/cell_feature_',num2str(t),'.mat'];  %保存 %修改文件夹名称 
  save(filename,'ShapeArea','ShapePerimeter','ShapeLmajor','ShapeLminor','ShapeElogation','ShapeCircularity','ShapeEccentricity','ShapeSmoothness','ShapeVolume', ...
      'ShapeSurface','ShapeDryMass','ShapeDryMassDensity','ShapeSurface_to_Volume_Ratio','ShapeSphericity', ...
      'FoFMean','FoFMedian','FoFRange','FoFMinimum','FoFMaximum','FoFVariation','FoFSkewness','FoFKurtosis', ...
      'FoFEntropy','FoFUniformity','FoFPercentile10','FoFPercentile90','FoFInterquartileRange','FoFMAD','FoFrMAD','FoFRMS','FoFEnergy', ...
      'GLCMAutoc','GLCMContr','GLCMCorre','GLCMCprom','GLCMCshad','GLCMDissi','GLCMEnerg','GLCMEntro','GLCMHomom','GLCMMaxpr',...
      'GLCMSosvh','GLCMSavgh','GLCMSvarh','GLCMSenth','GLCMDvarh','GLCMDenth','GLCMInf1h','GLCMInf2h','GLCMIndnc','GLCMIdmnc',...
      'GGCMsmgrd','GGCMbigrd','GGCMgrasy','GGCMgdasy','GGCMenerg','GGCMgrmea','GGCMgdmea','GGCMgrvar','GGCMgdvar','GGCMcorre','GGCMgrent','GGCMgdent','GGCMentro','GGCMinert','GGCMhomog',...
      'GLRLMsre','GLRLMlre','GLRLMgln','GLRLMglnn','GLRLMrp','GLRLMrln','GLRLMrlnn','GLRLMrv','GLRLMglv','GLRLMlgre','GLRLMhgre','GLRLMsrlgle','GLRLMsrhgle','GLRLMlrlgle','GLRLMlrhgle');
end

%% 4.细胞相位图特征数据集汇总
close all,clear all,clc;

% 创建初始数据集
ShapeArea=[];
ShapePerimeter=[];
ShapeLmajor=[];
ShapeLminor=[];
ShapeElogation=[];
ShapeCircularity=[];
ShapeEccentricity=[];
ShapeSmoothness=[];
ShapeVolume=[];
ShapeSurface=[];
ShapeDryMass=[];
ShapeDryMassDensity=[];
ShapeSurface_to_Volume_Ratio=[];
ShapeSphericity=[];

FoFEnergy=[];
FoFEntropy=[];
FoFInterquartileRange=[];
FoFKurtosis=[];
FoFMAD=[];
FoFMaximum=[];
FoFMean=[];
FoFMedian=[];
FoFMinimum=[];
FoFPercentile10=[];
FoFPercentile90=[];
FoFRange=[];
FoFrMAD=[];
FoFRMS=[];
FoFSkewness=[];
FoFUniformity=[];
FoFVariation=[];

GLCMAutoc=[];
GLCMContr=[];
GLCMCorre=[];
GLCMCprom=[];
GLCMCshad=[];
GLCMDissi=[];
GLCMEnerg=[];
GLCMEntro=[];
GLCMHomom=[];
GLCMMaxpr=[];
GLCMSosvh=[];
GLCMSavgh=[];
GLCMSvarh=[];
GLCMSenth=[];
GLCMDvarh=[];
GLCMDenth=[];
GLCMInf1h=[];
GLCMInf2h=[];
GLCMIndnc=[];
GLCMIdmnc=[];

GGCMsmgrd=[];
GGCMbigrd=[];
GGCMgrasy=[];
GGCMgdasy=[];
GGCMenerg=[];
GGCMgrmea=[];
GGCMgdmea=[];
GGCMgrvar=[];
GGCMgdvar=[];
GGCMcorre=[];
GGCMgrent=[];
GGCMgdent=[];
GGCMentro=[];
GGCMinert=[];
GGCMhomog=[];

GLRLMsre=[];
GLRLMlre=[];
GLRLMgln=[];
GLRLMglnn=[];
GLRLMrp=[];
GLRLMrln=[];
GLRLMrlnn=[];
GLRLMrv=[];
GLRLMglv=[];
GLRLMlgre=[];
GLRLMhgre=[];
GLRLMsrlgle=[];
GLRLMsrhgle=[];
GLRLMlrlgle=[];
GLRLMlrhgle=[];

%% 遍历特征
folder_path='/Users/taylor/Downloads/Data/cell_data/Breast_cancer_subtype/HER2_SKBR3/feature_example'; %修改文件夹名称
List=dir(fullfile(folder_path,'*.mat'));
k=numel(List); % 获取文件夹中符合要求的文件个数
for i=1:k
    file_name{i}=List(i).name;
    temp=importdata(file_name{i});
        ShapeArea(i)=getfield(temp,'ShapeArea');
        ShapePerimeter(i)=getfield(temp,'ShapePerimeter');
        ShapeLmajor(i)=getfield(temp,'ShapeLmajor');
        ShapeLminor(i)=getfield(temp,'ShapeLminor');
        ShapeElogation(i)=getfield(temp,'ShapeElogation');
        ShapeCircularity(i)=getfield(temp,'ShapeCircularity');
        ShapeEccentricity(i)=getfield(temp,'ShapeEccentricity');
        ShapeSmoothness(i)=getfield(temp,'ShapeSmoothness');
        ShapeVolume(i)=getfield(temp,'ShapeVolume');
        ShapeSurface(i)=getfield(temp,'ShapeSurface');
        ShapeDryMass(i)=getfield(temp,'ShapeDryMass');
        ShapeDryMassDensity(i)=getfield(temp,'ShapeDryMassDensity');
        ShapeSurface_to_Volume_Ratio(i)=getfield(temp,'ShapeSurface_to_Volume_Ratio');
        ShapeSphericity(i)=getfield(temp,'ShapeSphericity');

        FoFMean(i)=getfield(temp,'FoFMean');
        FoFMedian(i)=getfield(temp,'FoFMedian');
        FoFRange(i)=getfield(temp,'FoFRange');
        FoFMinimum(i)=getfield(temp,'FoFMinimum');
        FoFMaximum(i)=getfield(temp,'FoFMaximum');
        FoFVariation(i)=getfield(temp,'FoFVariation');
        FoFSkewness(i)=getfield(temp,'FoFSkewness');
        FoFKurtosis(i)=getfield(temp,'FoFKurtosis');
        FoFEntropy(i)=getfield(temp,'FoFEntropy');
        FoFUniformity(i)=getfield(temp,'FoFUniformity');
        FoFPercentile10(i)=getfield(temp,'FoFPercentile10');
        FoFPercentile90(i)=getfield(temp,'FoFPercentile90');
        FoFInterquartileRange(i)=getfield(temp,'FoFInterquartileRange');
        FoFMAD(i)=getfield(temp,'FoFMAD');
        FoFrMAD(i)=getfield(temp,'FoFrMAD');
        FoFRMS(i)=getfield(temp,'FoFRMS');
        FoFEnergy(i)=getfield(temp,'FoFEnergy');
       
        GLCMAutoc(i)=getfield(temp,'GLCMAutoc');
        GLCMContr(i)=getfield(temp,'GLCMContr');
        GLCMCorre(i)=getfield(temp,'GLCMCorre');
        GLCMCprom(i)=getfield(temp,'GLCMCprom');
        GLCMCshad(i)=getfield(temp,'GLCMCshad');
        GLCMDissi(i)=getfield(temp,'GLCMDissi');
        GLCMEnerg(i)=getfield(temp,'GLCMEnerg');
        GLCMEntro(i)=getfield(temp,'GLCMEntro');
        GLCMHomom(i)=getfield(temp,'GLCMHomom');
        GLCMMaxpr(i)=getfield(temp,'GLCMMaxpr');
        GLCMSosvh(i)=getfield(temp,'GLCMSosvh');
        GLCMSavgh(i)=getfield(temp,'GLCMSavgh');
        GLCMSvarh(i)=getfield(temp,'GLCMSvarh');
        GLCMSenth(i)=getfield(temp,'GLCMSenth');
        GLCMDvarh(i)=getfield(temp,'GLCMDvarh');
        GLCMDenth(i)=getfield(temp,'GLCMDenth');
        GLCMInf1h(i)=getfield(temp,'GLCMInf1h');
        GLCMInf2h(i)=getfield(temp,'GLCMInf2h');
        GLCMIndnc(i)=getfield(temp,'GLCMIndnc');
        GLCMIdmnc(i)=getfield(temp,'GLCMIdmnc');

        GGCMsmgrd(i)=getfield(temp,'GGCMsmgrd');
        GGCMbigrd(i)=getfield(temp,'GGCMbigrd');
        GGCMgrasy(i)=getfield(temp,'GGCMgrasy');
        GGCMgdasy(i)=getfield(temp,'GGCMgdasy');
        GGCMenerg(i)=getfield(temp,'GGCMenerg');
        GGCMgrmea(i)=getfield(temp,'GGCMgrmea');
        GGCMgdmea(i)=getfield(temp,'GGCMgdmea');
        GGCMgrvar(i)=getfield(temp,'GGCMgrvar');
        GGCMgdvar(i)=getfield(temp,'GGCMgdvar');
        GGCMcorre(i)=getfield(temp,'GGCMcorre');
        GGCMgrent(i)=getfield(temp,'GGCMgrent');
        GGCMgdent(i)=getfield(temp,'GGCMgdent');
        GGCMentro(i)=getfield(temp,'GGCMentro');
        GGCMinert(i)=getfield(temp,'GGCMinert');
        GGCMhomog(i)=getfield(temp,'GGCMhomog');

        GLRLMsre(i)=getfield(temp,'GLRLMsre');
        GLRLMlre(i)=getfield(temp,'GLRLMlre');
        GLRLMgln(i)=getfield(temp,'GLRLMgln');
        GLRLMglnn(i)=getfield(temp,'GLRLMglnn');
        GLRLMrp(i)=getfield(temp,'GLRLMrp');
        GLRLMrln(i)=getfield(temp,'GLRLMrln');
        GLRLMrlnn(i)=getfield(temp,'GLRLMrlnn');
        GLRLMrv(i)=getfield(temp,'GLRLMrv');
        GLRLMglv(i)=getfield(temp,'GLRLMglv');
        GLRLMlgre(i)=getfield(temp,'GLRLMlgre');
        GLRLMhgre(i)=getfield(temp,'GLRLMhgre');
        GLRLMsrlgle(i)=getfield(temp,'GLRLMsrlgle');
        GLRLMsrhgle(i)=getfield(temp,'GLRLMsrhgle');
        GLRLMlrlgle(i)=getfield(temp,'GLRLMlrlgle');
        GLRLMlrhgle(i)=getfield(temp,'GLRLMlrhgle');

end

%保存特征
filename=['/Users/taylor/Downloads/Data/cell_data/Breast_cancer_subtype/HER2_SKBR3/feature_example/feature_example','.mat'];  %保存 %修改文件夹名称 
save(filename,'ShapeArea','ShapePerimeter','ShapeLmajor','ShapeLminor','ShapeElogation','ShapeCircularity','ShapeEccentricity','ShapeSmoothness','ShapeVolume', ...
      'ShapeSurface','ShapeDryMass','ShapeDryMassDensity','ShapeSurface_to_Volume_Ratio','ShapeSphericity', ...
      'FoFMean','FoFMedian','FoFRange','FoFMinimum','FoFMaximum','FoFVariation','FoFSkewness','FoFKurtosis', ...
      'FoFEntropy','FoFUniformity','FoFPercentile10','FoFPercentile90','FoFInterquartileRange','FoFMAD','FoFrMAD','FoFRMS','FoFEnergy', ...
      'GLCMAutoc','GLCMContr','GLCMCorre','GLCMCprom','GLCMCshad','GLCMDissi','GLCMEnerg','GLCMEntro','GLCMHomom','GLCMMaxpr',...
      'GLCMSosvh','GLCMSavgh','GLCMSvarh','GLCMSenth','GLCMDvarh','GLCMDenth','GLCMInf1h','GLCMInf2h','GLCMIndnc','GLCMIdmnc',...
      'GGCMsmgrd','GGCMbigrd','GGCMgrasy','GGCMgdasy','GGCMenerg','GGCMgrmea','GGCMgdmea','GGCMgrvar','GGCMgdvar','GGCMcorre','GGCMgrent','GGCMgdent','GGCMentro','GGCMinert','GGCMhomog',...
      'GLRLMsre','GLRLMlre','GLRLMgln','GLRLMglnn','GLRLMrp','GLRLMrln','GLRLMrlnn','GLRLMrv','GLRLMglv','GLRLMlgre','GLRLMhgre','GLRLMsrlgle','GLRLMsrhgle','GLRLMlrlgle','GLRLMlrhgle','-v7.3');