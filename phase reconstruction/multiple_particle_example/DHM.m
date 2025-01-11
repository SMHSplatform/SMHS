%%
clear;clc;close all
 
%% 读取图像

I=imread('sample.jpg');%读所拍的细胞全息图
holography=im2double(I(:,:,1));%选取单通道图像
I2=imread('background.jpg');%读所拍的背景全息图
bg=im2double(I2(:,:,1));%选取单通道图像
figure,imshow(holography);%显示所拍的细胞全息图

%% 傅里叶变换与居中
FFTbg0=fftshift(fft2(bg));% 背景图 傅里叶变换与居中
FFTholography0=fftshift(fft2(holography));% 细胞图 傅里叶变换与居中
figure,imshow(abs(FFTholography0),[0,max(max(abs(FFTholography0)))/500]);%显示细胞图频谱
figure,imshow(abs(FFTbg0),[0,max(max(abs(FFTbg0)))/500]);%显示背景图频谱
%[xc,yc]=ginput(1)        % 获取+1级频谱中心坐标

%% 绘像面滤波窗 这一大块程序就是为了在图上画一个小框框
xc=2464;yc=2158;                  % 要正确读取+1级频谱的中心坐标位置。修改数值
dx=190;%选取所需的频谱范围
dy=190;
hold on;                 
xt=[xc-dx xc-dx];
yt=[yc-dy yc+dy];     
plot(xt,yt,'r');
xt=[xc+dx xc+dx];
yt=[yc-dy yc+dy];
plot(xt,yt,'r');
xt=[xc-dx xc+dx];
yt=[yc-dy yc-dy];
plot(xt,yt,'r');
xt=[xc-dx xc+dx];
yt=[yc+dy yc+dy];
plot(xt,yt,'r');  
hold off;

%N=2456;
N=5120;
Ufi=zeros(N,N);% 预设一个窗口，准备将提取的成分移到图像中间
Ufibg=zeros(N,N);
if xc-dx<1% 这一大块程序就是为了边缘检测防止越界
    xb0=1; 
else
    xb0=xc-dx; 
end     % xc是滤波窗中心坐标，dx为滤波窗半宽度
if yc-dy<1
    yb0=1; 
else
    yb0=yc-dy;
end
if     xc+dx>N
    xb1=N;
else
    xb1=xc+dx;
end
if yc+dy>N
    yb1=N; 
else
    yb1=yc+dy;
end
Ufi(yb0-yc+N/2:yb1-yc+N/2,xb0-xc+N/2:xb1-xc+N/2)=FFTholography0(yb0:yb1,xb0:xb1);%像面滤波；提取Uf中左上角复振幅成分
Ufibg(yb0-yc+N/2:yb1-yc+N/2,xb0-xc+N/2:xb1-xc+N/2)=FFTbg0(yb0:yb1,xb0:xb1);%背景图一样的操作
Uf=Ufi;
Ufbg=Ufibg;
figure,imshow(log(Uf),[])
%% 首先获取振幅信息-(光强信息)
I2=ifft2(ifftshift(Uf));%傅里叶逆变换
I3=abs(I2);
I4=mat2gray(I3);
figure,imshow(I4);title("强度-振幅信息")

%% 对背景图的操作
I2bg=ifft2(ifftshift(Ufbg));%傅里叶逆变换
%I5=abs(I2bg);
%figure,imshow(mat2gray(I3-I5));title("强度-振幅信息")

%% 获取包裹相位信息
a=real(I2); %提取实部
b=imag(I2); %提取虚部
phase=atan2(imag(I2),real(I2)); %反正切函数   
figure;imshow(phase,[]);title("包裹相位图");
%figure;mesh(phase);

abg=real(I2bg);
bbg=imag(I2bg);
phase_bg=atan2(imag(I2bg),real(I2bg));   
figure;imshow(phase_bg,[]);title("背景包裹相位图");
%figure;mesh(phase_bg);

figure;imshow(phase-phase_bg,[]);title("包裹相位图（减去背景）");
phase_wrapped=phase-phase_bg;
%figure;mesh(phase_wrapped);

%% 最小二乘解包裹运算
[M,N]=size(phase_wrapped);                 %计算二维包裹相位的大小（行、列数）
dx=zeros(M,N);dy=zeros(M,N);   %预设包裹相位沿x方向和y方向的梯度
m=1:M-1; 
dx(m,:)=phase_wrapped(m+1,:)-phase_wrapped(m,:);       %计算包裹相位沿x方向的梯度
dx=dx-pi*round(dx/pi);         %去除梯度中的跳跃
n=1:N-1;
dy(:,n)=phase_wrapped(:,n+1)-phase_wrapped(:,n);       %计算包裹相位沿y方向的梯度
dy=dy-pi*round(dy/pi);         %去除梯度中的跳跃
p=zeros(M,N);p1=zeros(M,N);p2=zeros(M,N); %为计算ρnm作准备
m=2:M;
p1(m,:)=dx(m,:)-dx(m-1,:);     %计算Δgxnm-Δgx(n-1)m
n=2:N;
p2(:,n)=dy(:,n)-dy(:,n-1);     %计算Δgynm–Δgyn(m-1)
p=p1+p2;                       %计算ρnm
p(1,1)=dx(1,1)+dy(1,1);        %计算ρnm
n=2:N;
p(1,n)=dx(1,n)+dy(1,n)-dy(1,n-1);%赋值Neumann边界条件
m=2:M;
p(m,1)=dx(m,1)-dx(m-1,1)+dy(m,1);
pp=dct2(p)+eps;                %计算ρnm的DCT
fi=zeros(M,N);
for m=1:M                      %计算Φnm在DCT域的精确解
   for n=1:N  
      fi(m,n)=pp(m,n)/(2*cos(pi*(m-1)/M)+2*cos(pi*(n-1)/N)-4+eps);
   end
end
fi(1,1)=pp(1,1);                %赋值DCT域的Φ11
phs=idct2(fi);                  %用iDCT计算解包裹相位在空域中的值

figure,imshow(-phs,[]);title("解包裹相位") ;          %显示解包裹相位
%colormap (jet);
%figure,mesh(phs)                %显示解包裹相位
%colormap (jet);
%% 二次曲面拟合去除系统畸变相位
curve_phase = fitt(phs);            % 畸变项
%figure,imshow(curve_phase,[])        % 显示畸变相位
%figure,mesh(curve_phase) % 显示畸变相位
%colormap (jet);

phase_unwrap_no_curve = -(phs - curve_phase);   % 去除畸变的相位(被测物体的纯相位)
figure,imshow(phase_unwrap_no_curve,[]);title('去除畸变后的相位')
colormap (jet);
%figure,mesh(phase_unwrap_no_curve)        % 显示去畸变相位
%colormap (jet);

phase_unwrap_no_curve(find(phase_unwrap_no_curve<0))=0; %减相位负值
figure,imshow(phase_unwrap_no_curve,[]);colormap (jet);
%figure,mesh(phase_unwrap_no_curve);
%colormap (jet);

%T2=graythresh(phase_unwrap_no_curve);
%BW2=im2bw(phase_unwrap_no_curve,T2);%Otus阈值进行分割
%figure;imshow(BW2),title('Otus阈值进行分割');

%figure,mesh(phase_unwrap_no_curve) ;
%axis on;
%saveas(gcf,['Reconstruction_mesh',num2str(1)],'jpg');
%colormap (jet);
%saveas(gcf,['Reconstruction_mesh_color',num2str(1)], 'jpg');
%figure,imshow(phase_unwrap_no_curve,[]);
%axis on;
%set(gca,'xtick',0:100:3000); 
%set(gca,'ytick',0:100:3000); 
%saveas(gcf, ['Reconstruction_imshow',num2str(1)], 'jpg');
%colormap (jet);
%saveas(gcf,['Reconstruction_imshow_color',num2str(1)], 'jpg');