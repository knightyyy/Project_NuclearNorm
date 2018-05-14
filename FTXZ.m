%% 基于核范数的稀疏重构算法，二维
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
% cofx=1;%方位维欠采样倍数
% deltaX=lambda/2*cofx;%方位维采样间隔，粗略设置为半波长
% Nx=floor(D/deltaX);%确定总共的采样位置
% if mod(Nx,2)==0
%     Nx=Nx;
% else
%     Nx=Nx-1;
% end
Nx=67;
deltaX=(D/Nx);
TRX_pos=(-D/2)+(0:Nx-1)*deltaX;%确定方位维采样位置
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
    Kx(:,index_f)=linspace(-kxmax,kxmax,Nx);%将方位维波数域均匀平分
    Kz(:,index_f)=sqrt(4*k^2-Kx(:,index_f).^2);%确定距离维波数域频率
end
%% 对方位维进行傅里叶变换
for index_f=1:Nf
    echo=ECHO(:,index_f);
    kx=Kx(:,index_f);
    s_ftx(:,index_f)=exp(-j*kx*TRX_pos)*echo;
end
%% 进行相位补偿处理
Phasecomp=Kz*TRZ;
s_comp_FT=s_ftx.*exp(j*Phasecomp);
%% 目标平面绘制
% cofgrid=1;
% xgrid=(-Nx/2:Nx/2)*rhox*cofgrid;
% zgrid=(-Nf/2:Nf/2)*rhoz*cofgrid;
Lx=0.5;
Lz=0.3;
xgrid=linspace(-Lx/2,Lx/2,Nx);
zgrid=linspace(-Lz/2,Lz/2,Nf);
%% 进行傅里叶逆变换
% [M_comp,N_comp]=size(s_comp_FT);
% nn=2;%二阶
% d0=60;
% m_comp=fix(M_comp/2);     %一半位置取整
% n_comp=fix(N_comp/2);
% for index_M=1:M_comp
%     for index_N=1:N_comp
%         d_comp=sqrt((index_M-m_comp)^2+(index_N-n_comp)^2);   %
%         if d_comp==0
%             h=0;
%         else
%             h_comp=1/(1+(d_comp/d0)^(2*nn));    %
% %             h_comp=exp(-(d_comp*d_comp)/(2*d0*d0));
%         end
%         s_comp_FT(index_M,index_N)=h_comp*s_comp_FT(index_M,index_N);
%     end
% end
for index_f=1:Nf
    kx=Kx(:,index_f);
    temp=s_comp_FT(:,index_f);
    s_iftx(:,index_f)=exp(j*xgrid.'*kx.')*temp/Nx;
end
for index_X=1:Nx
    kz=Kz(index_X,:);
    temp=s_iftx(index_X,:);
    s_iftz(index_X,:)=temp*exp(1j.*kz.'*(zgrid))./Nf;
end
figure
imagesc(xgrid*100,zgrid*100,imrotate(abs((s_iftz)), 90));
hold on
plot(Ptar(:,1)*100,Ptar(:,2)*100,'ro','MarkerSize',8,'LineWidth',1.5);
xlabel('Azimuth(cm)','Fontname','Times New Roman','FontSize',14);
ylabel('Range(cm)','Fontname','Times New Roman','FontSize',14);
set(gca, 'YDir', 'normal');
save('D:\yyy\核范数成像\实部虚部仿真1\s_comp_FT.mat','s_comp_FT');
