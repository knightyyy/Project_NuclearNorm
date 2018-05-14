%% 基于核范数的稀疏重构算法，二维仿真(均匀降采样)
clc;
clear all;
close all;
%% 基本参数设置
c=3e8;%光速
f0=26.5e9;%初始频率
fcutoff=40e9;%截止频率
B=fcutoff-f0;%带宽
fc=(f0+fcutoff)/2;%中心频点
lambda=c/fc;%波长
TRZ=0.3;%目标平面与观测平面之间的距离
%% 合成孔径设置
D=0.3;%合成孔径长度
cofx=1;%方位维欠采样倍数
deltaX=lambda/2*cofx;%方位维采样间隔，粗略设置为半波长
Nx=floor(D/deltaX);%确定总共的采样位置
if mod(Nx,2)==0
    Nx=Nx;
    
else
    Nx=Nx-1;
end
TRX_pos_origin=(-D/2)+(0:Nx-1)*deltaX;%确定方位维采样位置
N_random=55;
Numberarray_Nu=sort([1 randperm(Nx-3,N_random-2) Nx-2]);%随机产生位置
Nx=length(Numberarray_Nu);
Nc=67;%表示要插值结果的数目
TRX_pos=TRX_pos_origin(Numberarray_Nu(1:end));
%% 扫频设置
Nf=40;%设扫频频点个数为80,满足距离维奈奎斯特采样定理
deltaF=B/Nf;%每个频点之间的距离
f=f0+(0:Nf-1)*deltaF;%扫频频点
%% 将天线方向角与合成孔径辐射角度进行比较
thetaX_annt=30*pi/180;%天线方向角
thetaX_span=2*asin(D/sqrt(TRZ^2+D^2));%目标方位向距离和探测器方位向确定的展开角度
thetaX=min(thetaX_annt,thetaX_span);%确定波数域展开角度
%% 分辨率分析
rhox=c*sqrt(TRZ^2+(D/2)^2)/(2*f0*D/2);%方位维分辨率
rhoz=c/2/B;%距离维分辨率
%% 目标点参数设置
Ptar=[-8,4,1;-6,0,1;-4,-4,1;
      -2,0,1;0,4,1;2,0,1;
      4,-4,1;6,0,1;8,4,1]*diag([rhox rhoz 1]);
Ntar=length(Ptar(:,1));%计算目标点的个数
%% 回波信号构建
ECHO=zeros(Nx,Nf);
j=sqrt(-1);
for index_X=1:Nx
    s=zeros(1,Nf);
    for index_tar=1:Ntar
        x_tar=Ptar(index_tar,1);
        z_tar=Ptar(index_tar,2);
        A=Ptar(index_tar,3);
        R=sqrt((TRX_pos(index_X)-x_tar)^2+(TRZ-z_tar)^2);
        s=s+A*exp(-j*2*pi*f*2*R/c);
    end
    ECHO(index_X,:)=s;
end
% figure,
% imagesc(real(ECHO));
% figure,
% surf(real(ECHO));
%% 确定波数域频率
K=2*pi*f/c;%波数域频率
Kx=zeros(Nx,Nf);
Kz=zeros(Nx,Nf);
for index_f=1:Nf
    k=K(index_f);
    kxmax=k*sin(thetaX/2);
    kxorign=linspace(-kxmax,kxmax,Nc);
    Kx(:,index_f)=kxorign(Numberarray_Nu(1:end));%将方位维波数域均匀平分
    Kz(:,index_f)=sqrt(4*k^2-Kx(:,index_f).^2);%确定距离维波数域频率
end
%% 对方位维进行傅里叶变换

 S=zeros(Nc,Nf);
for index_f=1:Nf
    echo=ECHO(:,index_f);
    kx=Kx(:,index_f);
    s_ftx(:,index_f)=exp(-j*kx*TRX_pos)*echo;

%% 利用核范数进行处理

nsplikes=Nx;%表示目标点的个数

kk=0:Nc-1;
T=Nx;
N_s=(Nc+1)/2;
interval=2;
% temp = 1:Nx+1;
temp=Numberarray_Nu;
k1 = temp(1:T)-1;
%% 利用核范数进行处理

cvx_solver sdpt3
cvx_begin sdp quiet
    variable Y(N_s,N_s) hermitian
    variable Z(N_s,N_s) hermitian
    variable u(Nc,1) complex 
    minimize 0.5*trace(Y)+0.5*trace(Z)
    subject to
        [Y, hankel(u(1:N_s),u(N_s:end)); 
         hankel(u(1:N_s),u(N_s:end))', Z] >=0
        u(k1+1)==s_ftx(:,index_f);    
cvx_end

 S(:,index_f)=u;

end
% xgrid=(-Nc/2:Nc/2)*rhox;
% zgrid=(-Nf/2:Nf/2)*rhoz;
%% 设定目标平面参数
Lx=0.5;
Lz=0.3;
xgrid=linspace(-Lx/2,Lx/2,Nc);
zgrid=linspace(-Lz/2,Lz/2,Nf);
%% 确定插值之后的波数域频率
for index_f=1:Nf
    k=K(index_f);
    kxmax=k*sin(thetaX/2);
    Kxvec(:,index_f)=linspace(-kxmax,kxmax,Nc);%将方位维波数域均匀平分
    Kzvec(:,index_f)=sqrt(4*k^2-Kxvec(:,index_f).^2);%确定距离维波数域频率
end
%% 进行相位补偿处理
Phasecomp=Kzvec*TRZ;
s_comp_Nu=S.*exp(j*Phasecomp);
% [M_comp,N_comp]=size(s_comp_Nu);
% nn=2;%二阶
% d0=50;
% m_comp=fix(M_comp/2);
% n_comp=fix(N_comp/2);
% for index_M=1:M_comp
%     for index_N=1:N_comp
%         d_comp=sqrt((index_M-m_comp)^2+(index_N-n_comp)^2);
%         if d_comp==0
%             h=0;
%         else
%             h_comp=1/(1+(d_comp/d0)^(2*nn));
% %             h_comp=exp(-(d_comp*d_comp)/(2*d0*d0));
%         end
%         s_comp_Nu_non(index_M,index_N)=h_comp* s_comp_Nu(index_M,index_N);
%     end
% end
%% 进行二维傅里叶变换
for index_f=1:Nf    
    kxvec=Kxvec(:,index_f);
    temp=s_comp_Nu(:,index_f);
    s_iftx(:,index_f)=exp(j*xgrid.'*kxvec.')*temp/Nc;
end
for index_X=1:Nc
    kzvec=Kzvec(index_X,:);
    temp=s_iftx(index_X,:);
    s_iftz(index_X,:)=temp*exp(j.*kzvec.'*(zgrid))./Nf;
end

figure
%imagesc(xgrid*100,zgrid*100,imrotate(abs((s_iftz).'),-90));%降采样x
imagesc(xgrid*100,zgrid*100,imrotate(abs(s_iftz),90));%降采样f
set(gca, 'YDir', 'normal');
hold on
plot(Ptar(:,1)*100,Ptar(:,2)*100,'ro','MarkerSize',8,'LineWidth',1.5);
xlabel('Azimuth(cm)','Fontname','Times New Roman','FontSize',14);
ylabel('Range(cm)','Fontname','Times New Roman','FontSize',14);
set(gca, 'YDir', 'normal');



figure
%imagesc(xgrid*100,zgrid*100,imrotate(abs((s_iftz).'),-90));%降采样x
imagesc(xgrid*100,zgrid*100,imrotate(abs(s_iftz),90)); %降采样f
set(gca, 'YDir', 'normal');
hold on
plot(Ptar(:,1)*100,Ptar(:,2)*100,'ro','MarkerSize',8,'LineWidth',1.5);
xlabel('Azimuth(cm)','Fontname','Times New Roman','FontSize',14);
ylabel('Range(cm)','Fontname','Times New Roman','FontSize',14);
colormap(gray);
set(gca, 'YDir', 'normal');

for index_X=1:Nc
    for index_f=1:Nf 
        if(abs(s_iftz(index_X,index_f))>1.9)
            s_iftz(index_X,index_f)=s_iftz(index_X,index_f);
        else  s_iftz(index_X,index_f)=0;
        end
    end
end
%  

% Y=imread('C:\Users\303_thz\Desktop\33.png');
% M1= double(rgb2gray(Y))/255;
% J = grayslice(M1,30);%原灰度图灰度分成16层，J是索引图
% figure,imshow(J,jet(20));%figure,imshow(J,hot(16));

%频域高斯低通滤波和频域锐化的Matlab实现
I=double(s_iftz); 
nn=2;
[M,N]=size(I);  
m=(M+1)/2;  
n=(N+1)/2;  
d0=60;                      %截止频率  
F=fftshift(fft2(I));        %空域转频域，平移中心  
for i = 1:M  
   for j = 1:N  
        d = sqrt((i-m)^2+(j-n)^2);  
       if(d==0)  
            h=0;  
       else  
            h = 1/(1+(d0/d)/(2*nn));  
       end  
        result(i,j) = h*F(i,j);  
   end  
  end  
result = ifftshift(result);  
J2 = ifft2(result);  
% J3 = uint8(real(J2));  
figure
imagesc(xgrid*100,zgrid*100,imrotate(abs((J2)), 90));
hold on
plot(Ptar(:,1)*100,Ptar(:,2)*100,'ro','MarkerSize',8,'LineWidth',1.5);
xlabel('Azimuth(cm)','Fontname','Times New Roman','FontSize',14);   %字体和字号
ylabel('Range(cm)','Fontname','Times New Roman','FontSize',14);
colormap(gray);
set(gca, 'YDir', 'normal');

