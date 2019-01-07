clear all
tic
global Alpha Vc Sc K d_t t figure1
Alpha=1;
Vc=5*10^-16;
Sc=3*10^-10;
K=[ 1./3000   1./30000   1.e-4  3.e-11  1.e4]; % 6 kinetic constants
% Numerical discretization for Number Density Function
N=2001;
tmax=20000;
M=N;
d_t(1:N)=tmax/(M-1);
t=0:d_t(1):tmax;

%% Initial concentration of all proteins
c_init=[3.24*10^-10;0;0;0;0;0];% 'CD3','CD3i','CD25','IL 2','CD25-IL2'
%% Growth Function;
d=growth(c_init(1));
conc_Fig(d(:,:,:));
c3(1:M,1:M)=d(3,:,:);
c3=c3';
[c3min,t3min]=min(c3(1,:));
[c3max,t3max]=max(c3(1,:));
d_C=zeros(N,M);
% for j=1:M
    for i=2:N
        d_C(i-1,:)=(c3(i-1,:)-c3(i,:));
    end
% end
Lc=round(N/1);
d_c=(c3max-c3min)/(Lc-1);
C3=c3min:d_c:c3max;
t0=4;
%% Derivatives - Speed function
G=zeros(N,M);
for i=1:N-1
% for j=1:M
G(i,:)=d_C(i,:)./d_t(i);
% end
end
%% Analytical Solution
B1=@(t) 3*10^8*exp(-t/7200);
% B2=@(t) 3*10^5*exp(-(t-t_thld)/7200);
n1=zeros(N,M);
% n2=zeros(N,M);
% n12=zeros(N,M);
n=zeros(N,M);
freq=round(M/100);
% figure
% hold on
for j=1:freq:M
    cm=zeros(M,2);
     for i=1:j-1
        cm(i,1)=min(c3(i,j),c3(i+1,j));
        cm(i,2)=max(c3(i,j),c3(i+1,j));
        taum=0.5*(t(i)+t(i+1));
    if not(G(i,j)==0)
%         if (or(c_thld<=c3(i,j),G(i,j)>0))
           n1(j-i,j)=B1(taum)/abs(G(i,j));
%            if G(i,j)>0
%                plot(t(j),c3(i,j),'b');
%            else
%                plot(t(j),c3(i,j),'r');
%            end
%         end
    end
%         if t(i)>t_thld
%             n2(j-i,j)=B2(taum)/abs(G(i,j));
%         end
%         n12(j-i,j)=n1(j-i,j)+n2(j-i,j);
    end
for jt=1:M
    n(jt,j)=0;
    for jtau=1:j-1
        if and(C3(jt)>=cm(jtau,1), C3(jt)<cm(jtau,2))
           n(jt,j)=n(jt,j)+n1(j-jtau,j);
        end
    end
end
end
%% Match the c3 values and make new nA. 
figure
hold on
for k=1:freq:M
    plot(C3(3:N),n(3:N,k));
    pause(0.01);
end
for k=1:4
  tvl=tmax*[0.2 0.5 0.8 1];
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
title(['Population Density of T-cells for CD25-IL2 Protein',sprintf('\n'),'Analytical Method'],...
    'FontSize',12);

% Create xlabel
xlabel('State variable "c_3"','FontSize',16,'Color',[0 0 0]);

% Create ylabel
ylabel('Population Density','FontSize',16,'Color',[0 0 0]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
j=1;
k=0;
m=length(C3(1:Lc));%%In concentration variable, c3, values are repeating after m;
PBE_Figure(C3(i:j:m+k),n(i:j:m+k,len_t(1)),n(i:j:m+k,len_t(2)),n(i:j:m+k,len_t(3)),n(i:j:m+k,len_t(4)),len_t);

toc