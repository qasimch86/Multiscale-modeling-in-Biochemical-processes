%% Graph for concentration
function conc_Fig(C)
global t
time=t;
npas=length(time);
nvar=length(C(:,1,1));
figure
nfreq=floor(npas/100);
for nupro=1:nvar
 subplot(ceil(nvar/2),ceil(nvar/3),nupro) 
  for nuc=1:nfreq:npas
    plot(time,C(nupro,1:length(time),nuc))
    hold on
  end
    axis tight;
    grid
    hold off
end
