clear all
tic
global Alpha Vc Sc K d_t t figure1
Alpha=1;
Vc=5*10^-16;
Sc=3*10^-10;
K=[ 1./3000   1./30000   1.e-4  3.e-11  1.e4]; % 6 kinetic constants
% Numerical discretization for Number Density Function
N=2001;
tmax=24000;
M=N;
d_t(1:N)=tmax/(M-1);
t=0:d_t(1):tmax;

%% Initial concentration of all proteins
c_init=[3.24*10^-10;0;0;0;0;0];% 'CD3','CD3i','CD25','IL 2','CD25-IL2'
%% Growth Function;
d=growth(c_init(1));
conc_Fig(d(:,:,:));
c1(1:M,1:M)=d(1,:,:);
% c1=c1';
% [c1min,t1min]=min(min(c1));
% [c1max,t1max]=max(max(c1));
% for j=1:M
%     for i=2:N
%         d_C(i-1,j)=abs(c1(i,j)-c1(i-1,j));
%     end
% end
Lc=round(N);
C1=sort(c_init(1).*exp(-K(1).*t));%sort(c1(:,1));
t0=2;
d_C=diff(c1(:,1));
%% Derivatives - Speed function
for i=1:N-1
    for j=1:M
        G(i)=d_C(i)/d_t(j);
    end
end

% for j=1:M
%     G1(1:N-1,j)=G(N-1:-1:1,j);
% end
%% Analytical Solution
k=0;
n1=zeros(N,M);
B=@(t) 3*10^8*exp(-t/7200);
n=zeros(N,M);
freq=round(M/N);
for j=1:freq:M
    for i=1:j-1
%         if not(G1(i)==0)
           n(Lc-i,j)=B(t(j-i))/abs(G(i));
%         end
    end
end
%% Match the c3 values and make new nA. 
figure
% hold on
for k=1:freq*10:M-1
    plot(C1(1:Lc),n(:,k));
    pause(0.01);
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
title(['Population Density of T-cells for CD3 Protein',sprintf('\n'),'Analytical Method'],...
    'FontSize',12);

% Create xlabel
xlabel('State variable "c_1"','FontSize',16,'Color',[0 0 0]);

% Create ylabel
ylabel('Population Density','FontSize',16,'Color',[0 0 0]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
j=1;
k=0;
m=length(C1(1:Lc));%%In concentration variable, c1, values are repeating after m;
PBE_Figure(C1(i:j:m+k),n(i:j:m+k,len_t(1)),n(i:j:m+k,len_t(2)),n(i:j:m+k,len_t(3)),n(i:j:m+k,len_t(4)),len_t);

toc