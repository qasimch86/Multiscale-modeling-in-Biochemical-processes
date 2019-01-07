clear all
tic
global Alpha Vc Sc K d_t t figure1
Alpha=1;
Vc=5*10^-16;
Sc=3*10^-10;
K=[ 1./3000   1./30000   1.e-4  1e-9  1e4  1e-3]; % 6 kinetic constants
% Numerical discretization for Number Density Function
N=1001;
tmax=40000;
M=N;
d_t(1:N)=tmax/(M-1);
t=0:d_t(1):tmax;
%% Initial concentration of all proteins
c_init=[3.24*10^-10;0;0;0;0;0];% 'CD3','CD3i','CD25','IL 2','CD25-IL2'
%% Growth Function;
d=growth(c_init(1));
conc_Fig(d(:,:,:));
c1(1:M,1:M)=d(1,:,:);
[c1min,t1min]=min(min(c1));
[c1max,t1max]=max(max(c1));
%% Interpolation of c3 as x3, c4 as x4 and c5 as x5
for k=1:M
   G(k,:)=-K(1).*c1(k,:);
end
%%  Reducing Problem
% discretization for c
Lc=round((N-1)/5);
dC1=abs(c1min-c1max)/(Lc-1);
C1=c1min:dC1:c1max;
% Discretization of t for interpolated data
Lt=M;
dt=tmax/(Lt-1);
t1=0:dt:tmax;
t0=2;
%% Interpolation of new points
G1=zeros(Lc,Lt);
G1(1,t0-1)=G(t0-1,1);
for j=t0:Lt
%     cv=M-j+1;
    cmin=min(c1(j,1:j));
    cmax=max(c1(j,1:j));
    c=cmin:dC1:cmax;
    k=length(c);
    G1(1:k,j)=interp1(c1(j,j:-1:1),G(j,j:-1:1),c(1:k),'cubic','extrap');
end
for j=1:Lc
    G2(Lc-j+1,:)=G1(j,:);
end
%% Given condition
k=0;
n1=zeros(Lc,Lt);
B=@(t1) 3*10^8*exp(-t1/7200);
n1=zeros(Lc,Lt);
for i=1:Lc
   n1(i,1)=0;
end
for j=1:Lt
        n1(Lc,j)=-B(t1(j))/G2(Lc,j);
end
    
%% Numerical Method
for j=1:Lt-1
    for i=2:Lc-1
        v=dt/dC1;
% Flux at 3 points
        F0=G2(i-1,j)*n1(i-1,j);
        F1=G2(i,j)*n1(i,j);
        F2=G2(i+1,j)*n1(i+1,j);
% mean growth function
    if i==2
        if not(n1(i,j)==n1(i-1,j))
            G3(i-1,j)=v*(F1-F0)/(n1(i,j)-n1(i-1,j));
        else
            G3(i-1,j)=v*G2(i-1,j);
        end
        r(i-1)=1;%(n1(i,j)-n1(i-1,j))/(n1(i,j)-n1(i-1,j));
        PHI(i-1,j)=max(0,max(min(2*r(i-1),1),min(r(i-1),2)));%max(0,min(1,r(i-1)));%(abs(r(i-1))+r(i-1))/(1+abs(r(i-1)));   
        if G3(i-1,j)>0
            F(i-1,j)=F0+0.5*(1-G3(i-1,j))*(F1-F0)*PHI(i-1,j);
        elseif G3(i-1,j)<0
            F(i-1,j)=F1-0.5*(1+G3(i-1,j))*(F1-F0)*PHI(i-1,j);
        else
            F(i-1,j)=0;%(F0+F1-G3(i-1,j)*(F1-F0)*PHI(i-1,j))/2;
        end
    end
    
    if not(n1(i+1,j)==n1(i,j))
        G3(i,j)=v*(F2-F1)/(n1(i+1,j)-n1(i,j));
    else
        G3(i,j)=v*G2(i,j);
    end
    k=i-sign(G3(i,j));
    if k==Lc
        k=k-1;
    end
    r(i)=(n1(k+1,j)-n1(k,j))/(n1(i+1,j)-n1(i,j));
    PHI(i,j)=max(0,max(min(2*r(i),1),min(r(i),2)));%max(0,min(1,r(i)));%(abs(r(i))+r(i))/(1+abs(r(i)));%
% LWF Method with Flux limiter
    if G3(i,j)>0
        F(i,j)=F1+0.5*(1-G3(i,j))*(F2-F1)*PHI(i,j);
    elseif G3(i,j)<0
        F(i,j)=F2-0.5*(1+G3(i,j))*(F2-F1)*PHI(i,j);
    else
        F(i,j)=0;%(F1+F2-G3(i,j)*(F2-F1)*PHI(i,j))/2;
    end
    if not(G2(i+1,j)==0)
        n1(i,j+1)=n1(i,j)-v*(F(i,j)-F(i-1,j));  
    end
    end
end
%% Match the c3 values and make new nA.
% figure
% % hold on
% for j=1:5:Lt-t0+1
% %     C5=sort(c5(1:M,k),'descend');
%     plot(C1(1:Lc),n1(1:Lc,j));
%     pause(0.05);
% end
for k=1:4
  tvl=tmax*[0.3 0.6 0.8 1];
  len_t(k)=length(0:d_t:tvl(k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figure
figure1 = figure;

% Create axes
axes('Parent',figure1,'YTickLabel','','YTick',zeros(1,0),...
    'YColor',[0.8 0.8 0.8],...
    'XTickLabel','',...
    'XTick',zeros(1,0),...
    'XColor',[0.8 0.8 0.8],...
    'Position',[0.099 0.08269 0.8595 0.8145],...
    'CLim',[0 1]);

% Create title
title(['Population Density of T-cells for CD3 Protein',sprintf('\n'),'Method of Characteristics'],...
    'FontSize',12);

% Create xlabel
xlabel('State variable "c_1"','FontSize',16,'Color',[0 0 0]);

% Create ylabel
ylabel('Population Density','FontSize',16,'Color',[0 0 0]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
j=10;
k=0;
m=length(C1(1:Lc));%%In concentration variable, c5, values are repeating after m;
PBE_Figure(C1(i:j:m+k),n1(i:j:m+k,len_t(1)),n1(i:j:m+k,len_t(2)),n1(i:j:m+k,len_t(3)),n1(i:j:m+k,len_t(4)),len_t);
toc
