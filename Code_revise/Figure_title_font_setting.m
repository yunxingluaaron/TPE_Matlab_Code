
%% Temperature plotting
ylim([0 220]);
ly=ylabel('Temperature Change(K)');
lx=xlabel('Radius(m)');
set(gca, 'linewidth', 1.1)
set(gca, 'FontSize', 15)
set(lx,'FontSize',24)
set(ly,'FontSize',24)
box on;


%% Temperature plotting
ylim([0 45]);
%ly=ylabel('Temperature Change(K)');
ly=ylabel('\Delta T  (K)');
lx=xlabel('Radius(m)');
set(gca, 'linewidth', 1.1)
set(gca, 'FontSize', 15)
set(lx,'FontSize',24)
set(ly,'FontSize',24)
box on;



%% Pore pressure plotting
%ylim([0 3e7]);
%ylim([0 2.5e7]);
%ylim([0 6000]);
ylim([0 12e7]);
%ly=ylabel('\Delta p  (Pa)');
ly=ylabel('p  (Pa)');
lx=xlabel('Radius(m)');
set(gca, 'linewidth', 1.1)
set(gca, 'FontSize', 15)
set(lx,'FontSize',24)
set(ly,'FontSize',24)
box on;


%% Effective stress plotting
ylim([-5e7 1e7]);
ly=ylabel('\sigma^\prime_{r}   (Pa)');
lx=xlabel('Radius(m)');
set(gca, 'linewidth', 1.1)
set(gca, 'FontSize', 15)
set(lx,'FontSize',24)
set(ly,'FontSize',24)
box on;


x1 = [0 -5e7; 1 -5e7; 1 0e7;0 0];
y1 = [0 0; 1 0; 1 1e7;0 1e7];
fill(x1,'g')
fill(y1,'r')


%% Effective stress plotting
ylim([0 3.5e8]);
ylim([0 3e7]);
ylim([0 5200]);
xlim([1e-1 1e4]);
lx=xlabel('Time(days)');
ylim([-0.5e8 3e8]);

ylim([0 220]);


ly=ylabel('');
lx=xlabel('');

set(gca,'TickDir','out');
ly=ylabel('\sigma^\prime_{r}   (Pa)');
lx=xlabel('Radius(m)');