%% ���ں˷�����ϡ���ع��㷨����ά����
clc;
clear all;
close all;
%% ������������
c=3e8;%����
f0=26.5e9;%��ʼƵ��
fcutoff=40e9;%��ֹƵ��
B=fcutoff-f0;%����
fc=(f0+fcutoff)/2;%����Ƶ��
lambda=c/fc;%����
TRZ=0.3;%Ŀ��ƽ����۲�ƽ��֮��ľ���
%% �ϳɿ׾�����
D=0.3;%�ϳɿ׾�����
% cofx=1;%��λάǷ��������
% deltaX=lambda/2*cofx;%��λά�����������������Ϊ�벨��
% Nx=floor(D/deltaX);%ȷ���ܹ��Ĳ���λ��
% if mod(Nx,2)==0
%     Nx=Nx;
% else
%     Nx=Nx-1;
% end
Nx=67;
deltaX=(D/Nx);
TRX_pos=(-D/2)+(0:Nx-1)*deltaX;%ȷ����λά����λ��
%% ɨƵ����
Nf=40;%��ɨƵƵ�����Ϊ80,�������ά�ο�˹�ز�������
deltaF=B/Nf;%ÿ��Ƶ��֮��ľ���
f=f0+(0:Nf-1)*deltaF;%ɨƵƵ��
%% �����߷������ϳɿ׾�����ǶȽ��бȽ�
thetaX_annt=30*pi/180;%���߷����
thetaX_span=2*asin(D/sqrt(TRZ^2+D^2));%Ŀ�귽λ������̽������λ��ȷ����չ���Ƕ�
thetaX=min(thetaX_annt,thetaX_span);%ȷ��������չ���Ƕ�
%% �ֱ��ʷ���
rhox=c*sqrt(TRZ^2+(D/2)^2)/(2*f0*D/2);%��λά�ֱ���
rhoz=c/2/B;%����ά�ֱ���
%% Ŀ����������
Ptar=[-8,4,1;-6,0,1;-4,-4,1;
      -2,0,1;0,4,1;2,0,1;
      4,-4,1;6,0,1;8,4,1]*diag([rhox rhoz 1]);
Ntar=length(Ptar(:,1));%����Ŀ���ĸ���
%% �ز��źŹ���
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
%% ȷ��������Ƶ��
K=2*pi*f/c;%������Ƶ��
Kx=zeros(Nx,Nf);
Kz=zeros(Nx,Nf);
for index_f=1:Nf
    k=K(index_f);
    kxmax=k*sin(thetaX/2);
    Kx(:,index_f)=linspace(-kxmax,kxmax,Nx);%����λά���������ƽ��
    Kz(:,index_f)=sqrt(4*k^2-Kx(:,index_f).^2);%ȷ������ά������Ƶ��
end
%% �Է�λά���и���Ҷ�任
for index_f=1:Nf
    echo=ECHO(:,index_f);
    kx=Kx(:,index_f);
    s_ftx(:,index_f)=exp(-j*kx*TRX_pos)*echo;
end
%% ������λ��������
Phasecomp=Kz*TRZ;
s_comp_FT=s_ftx.*exp(j*Phasecomp);
%% Ŀ��ƽ�����
% cofgrid=1;
% xgrid=(-Nx/2:Nx/2)*rhox*cofgrid;
% zgrid=(-Nf/2:Nf/2)*rhoz*cofgrid;
Lx=0.5;
Lz=0.3;
xgrid=linspace(-Lx/2,Lx/2,Nx);
zgrid=linspace(-Lz/2,Lz/2,Nf);
%% ���и���Ҷ��任
% [M_comp,N_comp]=size(s_comp_FT);
% nn=2;%����
% d0=60;
% m_comp=fix(M_comp/2);     %һ��λ��ȡ��
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
save('D:\yyy\�˷�������\ʵ���鲿����1\s_comp_FT.mat','s_comp_FT');