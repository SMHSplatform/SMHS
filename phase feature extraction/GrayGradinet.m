function [GGCMout]=GrayGradinet(IMG)
% �Ҷ��ݶȹ������� H
% ��һ���Ҷ��ݶȾ��� H_basic
GGCMout.smgrd=zeros(1,1); % С�ݶ����� T1 Small Grads Dominance
GGCMout.bigrd=zeros(1,1); % ���ݶ����� T2 Big Grads Dominance
GGCMout.grasy=zeros(1,1); % �Ҷȷֲ��Ĳ������� T3 Gray Asymmetry
GGCMout.gdasy=zeros(1,1); % �ݶȷֲ��Ĳ������� T4 Grads Asymmetry
GGCMout.energ=zeros(1,1); % ���� T5 Energy
GGCMout.grmea=zeros(1,1); % �Ҷ�ƽ�� T6 Gray Mean
GGCMout.gdmea=zeros(1,1); % �ݶ�ƽ�� T7 Grads Mean
GGCMout.grvar=zeros(1,1); % �ҶȾ����� T8 Gray Variance
GGCMout.gdvar=zeros(1,1); % �ݶȾ����� T9 Grads Variance
GGCMout.corre=zeros(1,1); % ��� T10 Correlation 
GGCMout.grent=zeros(1,1); % �Ҷ��� T11 Gray Entropy
GGCMout.gdent=zeros(1,1); % �ݶ��� T12 Grads Entropy
GGCMout.entro=zeros(1,1); % ����� T13 Entropy
GGCMout.inert=zeros(1,1); % ���� T14 Inertia
GGCMout.homog=zeros(1,1); % ���� T15 Homogeneity

%ͼ������
IMG=im2uint8(IMG);
  %figure,imshow(IMG);
gray=256;
[R,C]=size(IMG);
%����ƽ����ͼ����ݶȾ���
GM=zeros(R-1,C-1);
for i=1:R-1
    for j=1:C-1
        n_GM=(IMG(i,j+1)-IMG(i,j))^2+(IMG(i+1,j)-IMG(i,j))^2;
        GM(i,j)=sqrt(double(n_GM));
    end
end
%figure,imshow(GM);

%�ҳ����ֵ��Сֵ        
n_min=min(GM(:));
n_max=max(GM(:));
%���ݶ�ͼ��Ҷȼ���ɢ��
%�����µĻҶȼ�Ϊnew_gray
new_gray=256;
%�µ��ݶȾ���Ϊnew_GM
new_GM=zeros(R-1,C-1);
new_GM=uint8((GM-n_min)/(n_max-n_min)*(new_gray-1));
%figure,imshow(new_GM);

%����Ҷ��ݶȹ�������
%�ݶȾ���ȹ�Ⱦ���ά����1�����ԻҶȾ�������Χ
H=zeros(gray,new_gray);
for i=1:R-1
    for j=1:C-1
        H(IMG(i,j)+1,new_GM(i,j)+1)= H(IMG(i,j)+1,new_GM(i,j)+1)+1;
    end
end

%��һ���Ҷ��ݶȾ��� H_basic
total=i*j;
H_basic=H/total;
%С�ݶ����� T1
TT=sum(H);
T1=0;
for j=1:new_gray
    T1=T1+TT(1,j)/j^2;
end
T1=T1/total;
%������ݶ����� T2
T2=0;
for j=1:new_gray
    T2=T2+TT(1,j)*(j-1);
end
T2=T2/total;
%����Ҷȷֲ��Ĳ������� T3
T3=0;
TT1=sum(H');
for j=1:gray
    T3=T3+TT1(1,j)^2;
end
T3=T3/total;
%�����ݶȷֲ��Ĳ������� T4
T4=0;
for j=1:new_gray
    T4=T4+TT(1,j)^2;
end
T4=T4/total;
%�������� T5
T5=0;
for i=1:gray
    for j=1:new_gray
        T5=T5+H_basic(i,j)^2;
    end
end
%����Ҷ�ƽ�� T6
TT2=sum((H_basic)');
T6=0;
for j=1:gray
    T6=T6+(j-1)*TT2(1,j);
end
%�����ݶ�ƽ�� T7
T7=0;
TT3=sum(H_basic);
for j=1:new_gray
    T7=T7+(j-1)*TT3(1,j);
end
%����ҶȾ����� T8
T8=0;
for j=1:gray
    T8=T8+(j-1-T6)^2*TT2(1,j);
end
T8=sqrt(T8);
%�����ݶȾ����� T9
T9=0;
for j=1:new_gray
    T9=T9+(j-1-T7)^2*TT3(1,j);
end
T9=sqrt(T9);
% ������� T10
T10=0;
for i=1:gray
    for j=1:new_gray
        T10=T10+(i-1-T6)*(j-1-T7)*H_basic(i,j);
    end
end
%����Ҷ��� T11
T11=0;
for j=1:gray
    T11=T11+TT2(1,j)*log10(TT2(1,j)+eps);
end
T11=-T11;
%�����ݶ��� T12
T12=0;
for j=1:new_gray
    T12=T12+TT3(1,j)*log10(TT3(1,j)+eps);
end
T12=-T12;
%�������� T13
T13=0;
for i=1:gray
    for j=1:new_gray
        T13=T13+H_basic(i,j)*log10(H_basic(i,j)+eps);
    end
end
T13=-T13;
%������� T14
T14=0;
for i=1:gray
    for j=1:new_gray
        T14=T14+(i-j)^2*H_basic(i,j);
    end
end
%�������� T15
T15=0;
for i=1:gray
    for j=1:new_gray
        T15=T15+H_basic(i,j)/(1+(i-j)^2);
    end
end

%x=1:50:750;

%�������
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
