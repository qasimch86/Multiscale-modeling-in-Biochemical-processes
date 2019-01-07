clear all
tic
global Alpha Vc Sc K d_t t figure1
Alpha=1;
Vc=5*10^-16;
Sc=3*10^-10;
K=[ 1./3000   1./30000   1.e-4  3.e-11  1.e4]; % 6 kinetic constants
% Numerical discretization for Number Density Function
N=4001;
tmax=24000;
M=N;
d_t(1:N)=tmax/(M-1);
t=0:d_t(1):tmax;
%% Initial concentration of all proteins
c_init=[3.24*10^-10;0;0;0;0;0];% 'CD3','CD3i','CD25','IL 2','CD25-IL2'
%% Growth Function;
d=growth(c_init(1));
conc_Fig(d(:,:,:));
c2(1:M,1:M)=d(2,:,:);
c3(1:M,1:M)=d(3,:,:);
c4(1:M)=d(4,:,1);
c3min=min(min(c3));
c3max=max(max(c3));
for k=1:M
   G(k,:)=K(3).*Vc./Sc.*c2(k,:)-K(5).*c3(k,:).*c4(k);
end
[cm,tm]=max(c3(:,1));
%%  Reducing Problem
% discretization for c
Lc=round((N-1)/15);
dC3=abs(c3min-c3max)/(Lc-1);
C3=c3min:dC3:c3max;
D3(1)=C3(1);
for i=2:Lc
    D3(i)=(C3(i)+C3(i-1))/2;
    dD3(i-1)=abs(D3(i)-D3(i-1));
end
% Discretization of t for interpolated data
Lt=N;
dt=tmax/(Lt-1);
t1=0:dt:tmax;
t0=3;
%% Interpolation of new points
k=1;
G1(1,1)=G(1,1);
for j=t0:Lt
    cm=max(c3(j,:));
    c=0:dC3:cm;
    k1(j)=length(c);
G1(1:k1(j),j)=interp1(c3(j,1:j-1),G(j,1:j-1),c(1:k1(j)),'cubic','extrap');
c=0;
for i=1:Lc-1
        if c(i)+dD3(i)>cm
            break;
        end
        c(i+1)=c(i)+dD3(i);
    end
    k2(j)=length(c);
G2(1:k2(j),j)=interp1(c3(j,1:j-1),G(j,1:j-1),c(1:k2(j)),'cubic','extrap');
end
for i=2:Lc
    G3(i,:)=(G1(i-1,:)+G1(i,:))/2;
end
G3(1,:)=G1(1,:);
%% Given condition
k=0;
n1=zeros(Lc,Lt);
B=@(t1) 3*10^8*exp(-t1/7200);
% figure
% hold on
double(n1);
n1=zeros(Lc,Lt);
for i=1:Lc
   n1(i,1)=0;
end
for j=t0:Lt
%     if not(G2(1,j)==0)
        n1(1,j)=B(t1(j))/G1(1,j);
%     end
end
% G1(Lc,:)=0;
% for j=1:Lt
%% Numerical Method
for j=1:M-1
    S(j)=n1(1,j);
    for i=2:Lc-1
        v=dt/dC3;
        F0=G2(i-1,j)*n1(i-1,j);
        F1=G2(i,j)*n1(i,j);
        F2=G2(i+1,j)*n1(i+1,j);
        if G2(i-1,j)>0
            n1(i,j+1)=n1(i,j)-v*(F2-F0)/2+0.5*v^2*(G3(i,j)*(F2-F1)-G3(i-1,j)*(F1-F0));
        elseif G2(i,j)<0
            n1(i,j+1)=n1(i,j)-v*(F2-F0)/2+0.5*v^2*(G3(i+1,j)*(F2-F1)-G3(i,j)*(F1-F0));
        end
    S(j)=S(j)+n1(i,j);
    end
end
%% Match the c3 values and make new nA.
figure
hold on
for j=1:10:M
%     C5=sort(c5(1:M,k),'descend');
    plot(C3(1:Lc-1),n1(1:Lc-1,j));
    pause(0.05);
end
for k=1:4
  tvl=tmax*[0.125 0.375 0.625 1];
  len_t(k)=length(0:d_t:tvl(k));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
title(['Population Density of T-cells for CD25 Protein',sprintf('\n'),'Lax-Wendroff Scheme'],...
    'FontSize',12);

% Create xlabel
xlabel('State variable "c_5"','FontSize',16,'Color',[0 0 0]);

% Create ylabel
ylabel('Population Density','FontSize',16,'Color',[0 0 0]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
j=2;
k=0;
m=length(C3(1:Lc));%%In concentration variable, c5, values are repeating after m;
PBE_Figure(C3(i:j:m+k),n1(i:j:m+k,len_t(1)),n1(i:j:m+k,len_t(2)),n1(i:j:m+k,len_t(3)),n1(i:j:m+k,len_t(4)),len_t);

toc