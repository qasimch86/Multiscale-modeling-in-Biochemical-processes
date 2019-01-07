% Main program
% solving the set of integro differential equations
% C : matrix of concentrations C(i,n,k): ith protein, time n, activation k
% we may keep only 2 values fo k: previous and next time step
% kinetic constants
% close all

% Kinetic constants
function C=growth(c0)
global Alpha Vc Sc K d_t t
%K(5)=K(5)*0.01;
% Sc= 3e-10;
% Vc= 5e-16;
% r1=1/300;%2;%k1
% r2=1/30000;%0.4;%k2
% r3=1/10000;%2;%k3
% r4=10^-9;%0.8;%k4
% r5=10^6;%3;%k5
% %  r5p=0.03;%3;%k5p
% r6=10^-3;%2;%k6 
% time
% tmax=20000;
time=t;
npas=length(time)
% concentration
% C1ini=3.24e-12;   % initial C1
C=zeros(6,npas,npas);
b=funcb(time);

% conditions initiales
C(1,1,1:npas)=c0;
 
% boucle sur le temps
temps=0;
n=1;
Time=1;
% while temps < tmax-1e-10
%     temps=temps+delt;
for n=2:npas
%     n=n+1;
    delt=d_t(n-1);%(tmax/(npas-1));
    C(1,n,n:npas)=c0;
%           x(1:6,:)=C(1:6,n,:);
      for k=1:n-1
%           RK Method (npas=1250 for tmax=100000)
%         DO something so that it runs RK somewhere and send value here
%           k1 = rhs(x(:,k),C(4,n-1,1));
%           k2 = rhs(x(:,k) + delt*k1./2,C(4,n-1,1));
%           k3 = rhs(x(:,k) + delt*k2./2,C(4,n-1,1));
%           k4 = rhs(x(:,k) + delt*k3,C(4,n-1,1));
%           y(:,k)=x(:,k)+ delt*(k1 + 2*k2 + 2*k3 + k4)/6;
%           C(1:3,n,k)=y(1:3,k); 
%           C(5:6,n,k)=y(5:6,k); 

%           Euler Method 
          C(1,n,k)=C(1,n-1,k)-delt*K(1)*C(1,n-1,k);
          C(2,n,k)=C(2,n-1,k)+delt*K(1)*(Sc/Vc)*C(1,n-1,k)-delt*K(2)*C(2,n-1,k);
          C(4,n,1)=K(4)*time(n-1);
          C(3,n,k)=C(3,n-1,k)+delt*K(3)*(Vc/Sc)*C(2,n-1,k)-delt*K(5)*C(3,n-1,k)*C(4,n,1);          
     end 
%     d4(i+1)=d4(i)+d_t/2*(fj_11-fj_21*d4(i) + fj_12-fj_22*(d4(i)+d_t*(fj_11-fj_21*d4(i))));
          sum1=0;   %trapezes integration
          sum2=0;
          if n>1       
              sum1=delt*(0.5*C(2,n-1,1)*b(1) + C(2,n-1,n-1)*b(n-1))  ;
              sum2=delt*(0.5*C(3,n-1,1)*b(1) + C(3,n-1,n-1)*b(n-1))  ;
             if n>2
                 Cc2=reshape(C(2,n-1,1:n-2),1,n-2);
                 Cc3=reshape(C(3,n-1,1:n-2),1,n-2);
                 sum1=sum1+delt*sum(Cc2.*b(1:n-2));
                 sum2=sum2+delt*sum(Cc3.*b(1:n-2));
             end
          end
Time=Time+1;      
end
     C(1,n-1,n:npas)=c0;

% Plots
% figure
% nfreq=floor(npas/10);
% for nupro=1:6
%  subplot(2,3,nupro) 
%   for nuc=1:nfreq:npas
%     plot(time,C(nupro,1:length(time),nuc))
%     hold on
%   end
%     axis tight;
%     grid
%     hold off
% end