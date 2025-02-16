function [GGCMout]=GrayGradinet(IMG)
% 灰度梯度共生矩阵 H
% 归一化灰度梯度矩阵 H_basic
GGCMout.smgrd=zeros(1,1); % 小梯度优势 T1 Small Grads Dominance
GGCMout.bigrd=zeros(1,1); % 大梯度优势 T2 Big Grads Dominance
GGCMout.grasy=zeros(1,1); % 灰度分布的不均匀性 T3 Gray Asymmetry
GGCMout.gdasy=zeros(1,1); % 梯度分布的不均匀性 T4 Grads Asymmetry
GGCMout.energ=zeros(1,1); % 能量 T5 Energy
GGCMout.grmea=zeros(1,1); % 灰度平均 T6 Gray Mean
GGCMout.gdmea=zeros(1,1); % 梯度平均 T7 Grads Mean
GGCMout.grvar=zeros(1,1); % 灰度均方差 T8 Gray Variance
GGCMout.gdvar=zeros(1,1); % 梯度均方差 T9 Grads Variance
GGCMout.corre=zeros(1,1); % 相关 T10 Correlation 
GGCMout.grent=zeros(1,1); % 灰度熵 T11 Gray Entropy
GGCMout.gdent=zeros(1,1); % 梯度熵 T12 Grads Entropy
GGCMout.entro=zeros(1,1); % 混合熵 T13 Entropy
GGCMout.inert=zeros(1,1); % 惯性 T14 Inertia
GGCMout.homog=zeros(1,1); % 逆差矩 T15 Homogeneity

%图像输入
IMG=im2uint8(IMG);
  %figure,imshow(IMG);
gray=256;
[R,C]=size(IMG);
%采用平方求和计算梯度矩阵
GM=zeros(R-1,C-1);
for i=1:R-1
    for j=1:C-1
        n_GM=(IMG(i,j+1)-IMG(i,j))^2+(IMG(i+1,j)-IMG(i,j))^2;
        GM(i,j)=sqrt(double(n_GM));
    end
end
%figure,imshow(GM);

%找出最大值最小值        
n_min=min(GM(:));
n_max=max(GM(:));
%把梯度图象灰度级离散化
%设置新的灰度级为new_gray
new_gray=256;
%新的梯度矩阵为new_GM
new_GM=zeros(R-1,C-1);
new_GM=uint8((GM-n_min)/(n_max-n_min)*(new_gray-1));
%figure,imshow(new_GM);

%计算灰度梯度共生矩阵
%梯度矩阵比轨度矩阵维数少1，忽略灰度矩阵最外围
H=zeros(gray,new_gray);
for i=1:R-1
    for j=1:C-1
        H(IMG(i,j)+1,new_GM(i,j)+1)= H(IMG(i,j)+1,new_GM(i,j)+1)+1;
    end
end

%归一化灰度梯度矩阵 H_basic
total=i*j;
H_basic=H/total;
%小梯度优势 T1
TT=sum(H);
T1=0;
for j=1:new_gray
    T1=T1+TT(1,j)/j^2;
end
T1=T1/total;
%计算大梯度优势 T2
T2=0;
for j=1:new_gray
    T2=T2+TT(1,j)*(j-1);
end
T2=T2/total;
%计算灰度分布的不均匀性 T3
T3=0;
TT1=sum(H');
for j=1:gray
    T3=T3+TT1(1,j)^2;
end
T3=T3/total;
%计算梯度分布的不均匀性 T4
T4=0;
for j=1:new_gray
    T4=T4+TT(1,j)^2;
end
T4=T4/total;
%计算能量 T5
T5=0;
for i=1:gray
    for j=1:new_gray
        T5=T5+H_basic(i,j)^2;
    end
end
%计算灰度平均 T6
TT2=sum((H_basic)');
T6=0;
for j=1:gray
    T6=T6+(j-1)*TT2(1,j);
end
%计算梯度平均 T7
T7=0;
TT3=sum(H_basic);
for j=1:new_gray
    T7=T7+(j-1)*TT3(1,j);
end
%计算灰度均方差 T8
T8=0;
for j=1:gray
    T8=T8+(j-1-T6)^2*TT2(1,j);
end
T8=sqrt(T8);
%计算梯度均方差 T9
T9=0;
for j=1:new_gray
    T9=T9+(j-1-T7)^2*TT3(1,j);
end
T9=sqrt(T9);
% 计算相关 T10
T10=0;
for i=1:gray
    for j=1:new_gray
        T10=T10+(i-1-T6)*(j-1-T7)*H_basic(i,j);
    end
end
%计算灰度熵 T11
T11=0;
for j=1:gray
    T11=T11+TT2(1,j)*log10(TT2(1,j)+eps);
end
T11=-T11;
%计算梯度熵 T12
T12=0;
for j=1:new_gray
    T12=T12+TT3(1,j)*log10(TT3(1,j)+eps);
end
T12=-T12;
%计算混合熵 T13
T13=0;
for i=1:gray
    for j=1:new_gray
        T13=T13+H_basic(i,j)*log10(H_basic(i,j)+eps);
    end
end
T13=-T13;
%计算惯性 T14
T14=0;
for i=1:gray
    for j=1:new_gray
        T14=T14+(i-j)^2*H_basic(i,j);
    end
end
%计算逆差矩 T15
T15=0;
for i=1:gray
    for j=1:new_gray
        T15=T15+H_basic(i,j)/(1+(i-j)^2);
    end
end

%x=1:50:750;

%输出特征
GGCMout.smgrd=T1;
GGCMout.bigrd=T2;
GGCMout.grasy=T3;
GGCMout.gdasy=T4;
GGCMout.energ=T5;
GGCMout.grmea=T6;
GGCMout.gdmea=T7;
GGCMout.grvar=T8;
GGCMout.gdvar=T9;
GGCMout.corre=T10;
GGCMout.grent=T11;
GGCMout.gdent=T12;
GGCMout.entro=T13;
GGCMout.inert=T14;
GGCMout.homog=T15;

% if num>2
%     plot(x,OUT,'-');
%     hold on;
% else
%     plot(x,OUT,'-*r');
%     hold on;
% end
