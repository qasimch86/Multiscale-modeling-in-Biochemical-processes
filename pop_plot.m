%% Graph for Population density function
function pop_plot(len_t,c1,n,y,nA)
global t
hold on
i=10;
j=1;
k=-4; %% k<=0
m=length(c1);%%In concentration variable, c1, values are repeating after m;
plot(c1(i:j:m+k),n(i:j:m+k,len_t(1)),'g-','MarkerSize',4);
% plot(c1(i:j:m+k),y(i:j:m+k,len_t(1)),'g-o','MarkerSize',4);
% plot(c1(i:j:m+k),nA(i:j:m+k,len_t(1)),'g-o','MarkerSize',4);

plot(c1(i:j:m+k),n(i:j:m+k,len_t(2)),'b-','MarkerSize',4);
% plot(c1(i:j:m+k),y(i:j:m+k,len_t(2)),'b-*','MarkerSize',4);
% plot(c1(i:j:m+k),nA(i:j:m+k,len_t(2)),'b-*','MarkerSize',4);

plot(c1(i:j:m+k),n(i:j:m+k,len_t(3)),'m-','MarkerSize',4);
% plot(c1(i:j:m+k),y(i:j:m+k,len_t(3)),'m->','MarkerSize',4);
% plot(c1(i:j:m+k),nA(i:j:m+k,len_t(3)),'m->','MarkerSize',4);

plot(c1(i:j:m+k),n(i:j:m+k,len_t(4)),'k-','MarkerSize',4);
% plot(c1(i:j:m+k),y(i:j:m+k,len_t(4)),'k-s','MarkerSize',4);
% plot(c1(i:j:m+k),nA(i:j:m+k,len_t(4)),'k-s','MarkerSize',4);

% plot(c1(i:j:m+k),n(i:j:m+k,len_t(5)),'r-d','MarkerSize',4);
hold off
legend(strcat('Approx. Sol. t= ',num2str(t(len_t(1)))),strcat('Approx. Sol. t= ',num2str(t(len_t(2)))),strcat('Approx. Sol. t= ',num2str(t(len_t(3)))),strcat('Approx. Sol. t= ',num2str(t(len_t(4)))));%,strcat('Approx. Sol. t= ',num2str(t(len_t(5)))));
% legend(strcat('Approx. Sol. t=',num2str(len_t(1))','","Approx. Sol. t=',num2str(len_t(2)),'","Approx. Sol. t=',num2str(len_t(3)),'","Approx. Sol. t=',num2str(len_t(4)),'","Approx. Sol. t=',num2str(len_t(5)),'");
title('Population Density of T-cells for CD3 Protein','FontSize',12)
ylabel('Population Density','FontSize',16)
xlabel('State variable "c_2"','FontSize',16)
grid
axis tight
% axis ([0 max(c1(i:j:m+k)) 0 1e14]);