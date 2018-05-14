%% ���ں˷�����ϡ���ع��㷨����ά����(���Ƚ�����)
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
cofx=1;%��λάǷ��������
deltaX=lambda/2*cofx;%��λά�����������������Ϊ�벨��
Nx=floor(D/deltaX);%ȷ���ܹ��Ĳ���λ��
if mod(Nx,2)==0
    Nx=Nx;
    
else
    Nx=Nx-1;
end
TRX_pos_origin=(-D/2)+(0:Nx-1)*deltaX;%ȷ����λά����λ��
N_random=55;
Numberarray_Nu=sort([1 randperm(Nx-3,N_random-2) Nx-2]);%�������λ��
Nx=length(Numberarray_Nu);
Nc=67;%��ʾҪ��ֵ�������Ŀ
TRX_pos=TRX_pos_origin(Numberarray_Nu(1:end));
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
    kxorign=linspace(-kxmax,kxmax,Nc);
    Kx(:,index_f)=kxorign(Numberarray_Nu(1:end));%����λά���������ƽ��
    Kz(:,index_f)=sqrt(4*k^2-Kx(:,index_f).^2);%ȷ������ά������Ƶ��
end
%% �Է�λά���и���Ҷ�任

 S=zeros(Nc,Nf);
for index_f=1:Nf
    echo=ECHO(:,index_f);
    kx=Kx(:,index_f);
    s_ftx(:,index_f)=exp(-j*kx*TRX_pos)*echo;

%% ���ú˷������д���

nsplikes=Nx;%��ʾĿ���ĸ���

kk=0:Nc-1;
T=Nx;
N_s=(Nc+1)/2;
interval=2;
% temp = 1:Nx+1;
temp=Numberarray_Nu;
k1 = temp(1:T)-1;
%% ���ú˷������д���

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
%% �趨Ŀ��ƽ�����
Lx=0.5;
Lz=0.3;
xgrid=linspace(-Lx/2,Lx/2,Nc);
zgrid=linspace(-Lz/2,Lz/2,Nf);
%% ȷ����ֵ֮��Ĳ�����Ƶ��
for index_f=1:Nf
    k=K(index_f);
    kxmax=k*sin(thetaX/2);
    Kxvec(:,index_f)=linspace(-kxmax,kxmax,Nc);%����λά���������ƽ��
    Kzvec(:,index_f)=sqrt(4*k^2-Kxvec(:,index_f).^2);%ȷ������ά������Ƶ��
end
%% ������λ��������
Phasecomp=Kzvec*TRZ;
s_comp_Nu=S.*exp(j*Phasecomp);
% [M_comp,N_comp]=size(s_comp_Nu);
% nn=2;%����
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
%% ���ж�ά����Ҷ�任
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
%imagesc(xgrid*100,zgrid*100,imrotate(abs((s_iftz).'),-90));%������x
imagesc(xgrid*100,zgrid*100,imrotate(abs(s_iftz),90));%������f
set(gca, 'YDir', 'normal');
hold on
plot(Ptar(:,1)*100,Ptar(:,2)*100,'ro','MarkerSize',8,'LineWidth',1.5);
xlabel('Azimuth(cm)','Fontname','Times New Roman','FontSize',14);
ylabel('Range(cm)','Fontname','Times New Roman','FontSize',14);
set(gca, 'YDir', 'normal');



figure
%imagesc(xgrid*100,zgrid*100,imrotate(abs((s_iftz).'),-90));%������x
imagesc(xgrid*100,zgrid*100,imrotate(abs(s_iftz),90)); %������f
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
% J = grayslice(M1,30);%ԭ�Ҷ�ͼ�Ҷȷֳ�16�㣬J������ͼ
% figure,imshow(J,jet(20));%figure,imshow(J,hot(16));

%Ƶ���˹��ͨ�˲���Ƶ���񻯵�Matlabʵ��
I=double(s_iftz); 
nn=2;
[M,N]=size(I);  
m=(M+1)/2;  
n=(N+1)/2;  
d0=60;                      %��ֹƵ��  
F=fftshift(fft2(I));        %����תƵ��ƽ������  
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
xlabel('Azimuth(cm)','Fontname','Times New Roman','FontSize',14);   %������ֺ�
ylabel('Range(cm)','Fontname','Times New Roman','FontSize',14);
colormap(gray);
set(gca, 'YDir', 'normal');

