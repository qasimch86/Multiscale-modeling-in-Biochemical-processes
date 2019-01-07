clear all
tic
global Alpha Vc Sc K d_t t figure1
Alpha=1;
Vc=5*10^-16;
Sc=3*10^-10;
K=[ 1./3000   1./30000   1.e-4  3e-11  1e4]; % 6 kinetic constants
% Numerical discretization for Number Density Function
N=2001;
tmax=100000;
M=N;
d_t(1:N)=tmax/(M-1);
t=0:d_t(1):tmax;
%% Initial concentration of all proteins
c_init=[3.24*10^-10;0;0;0;0;0];% 'CD3','CD3i','CD25','IL 2','CD25-IL2'
% c=growth(c_init);
%% Finiding concentration for the CD3
d=growth(c_init(1));
conc_Fig(d(:,:,:));

c1(1:M,1:M)=d(1,:,:);
[c1min,t1min]=min(min(c1));
[c1max,t1max]=max(max(c1));
% c1=c1';
% c1max=max(c(1,:));
%%  Reducing Problem
% discretization for c
Lc=round((N));
% for i=1:N-1
dC1=abs(diff(c1(:,1)));
% end
C1=sort(c1(:,1));
Lt=N;
%% Exact Solution
h1=-1/K(1)*log(c1max);
for j=1:Lt
        x1(j)= 3.24*10^-10*exp(-K(1)*t(j));
        x1(j)=x1(j)*10^60;
        x1(j)=fix(x1(j));
        x1(j)=x1(j)/10^60;
        
    for i=1:Lc
%         %Analytical Method
%          if t(j)>t(N-i+1) && t(i)>=0     
%             nA(i,j)=-B(t(j)-t(N-i+1))/(G1(i));
%          end
%         % Method of Characteristics
         if (C1(i)<x1(j))
            n1(i,j)=0;
        elseif (C1(i)>=x1(j))
            n1(i,j)=3e+8*exp(-t(j)/7200-h1/7200-1/K(1)/7200*log(C1(i)))/(K(1)*C1(i)^(1));
         end
    end
end
% figure
% % hold on
% for j=1:10:Lt
% %     C5=sort(c5(1:M,k),'descend');
%     plot(C1(1:Lc),n1(1:Lc,j));
%     pause(0.05);
% end
%Plot at different times
for k=1:4
  tvl=tmax*[0.2 0.5 0.8 1];
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
j=25;
k=0;
m=length(C1(1:Lc));%%In concentration variable, c5, values are repeating after m;
PBE_Figure(C1(i:j:m+k),n1(i:j:m+k,len_t(1)),n1(i:j:m+k,len_t(2)),n1(i:j:m+k,len_t(3)),n1(i:j:m+k,len_t(4)),len_t);
toc