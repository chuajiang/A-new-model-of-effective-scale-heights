clear all;
close all;

%%%%matched with GNSS TEC
load('timeD.mat');
load('gTEC.mat');
load('timeY.mat');

load('TECdata_10.mat');
Vtec=TECdata_10.Vtec;
gTEC=TECdata_10.gTEC;
t=TECdata_10.t;
iv=TECdata_10.iv;
ig=TECdata_10.ig;
HFin=TECdata_10.HFin;
startTime = t(1);
endTime = t(end);

%%%%Figure 2
figure;
set(gcf,'unit', 'normalized', 'position',[0.1,0.1,0.9,0.5]);
plot(t, Vtec(iv), 'o', 'MarkerSize', 2, 'MarkerFaceColor','b');
hold on;
plot(t, gTEC(ig), 'r', 'LineWidth', 2);
hold on;
plot(t, Vtec(iv)-gTEC(ig), 'o', 'MarkerSize', 2, 'MarkerFaceColor','g');
xlim([startTime endTime]);
set(gca, 'xtick', startTime : datenum(0, 0, 0, 24, 0, 0) : endTime);
grid on;
datetick('x', 'dd');
% legend('\beta-Chapman', 'GNSS', 'Chapman-GNSS');
legend('Epstein', 'GNSS', 'Epstein-GNSS');
legend('boxoff');

xlabel('Day');
ylabel('VTEC/TECU');
title('Octorber, 2021');
set(gca,'FontSize',18, 'FontWeight','bold');

%%%%Figure 3
HFin(find(HFin==0)) = NaN;
HFin = HFin(iv);
figure;
set(gcf,'unit', 'normalized', 'position',[0.1,0.1,0.9,0.5]);
plot(t, HFin, 'o', 'MarkerSize', 2);
xlim([startTime endTime]);
set(gca, 'xtick', startTime : datenum(0, 0, 0, 24, 0, 0) : endTime);
x_label = 1001 : 1 : 1031;
x_label = [x_label 1101];
set(gca, 'xticklabel', x_label);
grid on;
datetick('x', 'dd');

xlabel('Day');
ylabel('Scale Height/km');
ylim([50 250]);
% yticks([40 80 120 160 200]);
title('Octorber, 2021');
set(gca,'FontSize',18, 'FontWeight','bold');



%%%Figure 4 monthly meadian value
% clear all;
% close all;
Htdata_2021=[];
for imon=1:12
    datafile=sprintf('Htdata_%02d.mat',imon);
    load(datafile);
   
end

Htdata_2021=[Htdata_01.month_Ht_data;Htdata_02.month_Ht_data;Htdata_03.month_Ht_data;Htdata_04.month_Ht_data;
           Htdata_05.month_Ht_data;Htdata_06.month_Ht_data;Htdata_07.month_Ht_data; Htdata_08.month_Ht_data;
           Htdata_09.month_Ht_data; Htdata_10.month_Ht_data; Htdata_11.month_Ht_data; Htdata_12.month_Ht_data];
F = Htdata_2021;
F(1, :) = fillmissing(F(1, :), 'makima', 2);
F(end, :) = fillmissing(F(end, :), 'makima', 2);
F(:, 1) = fillmissing(F(:, 1), 'makima', 2);
F(:, end) = fillmissing(F(:, end), 'makima', 2);
F = fillmissing(F, 'makima', 2);
F(F < 40) = NaN;
F(F > 250) = NaN;
F = fillmissing(F, 'makima', 1);
F(F < 40) = NaN;
F(F > 250) = NaN;
F = fillmissing(F, 'makima', 2);
F(F > 250) = NaN;
F = fillmissing(F, 'makima', 1);

Htdata_2021_makima=F;
save('Htdata_2021_makima.mat','Htdata_2021_makima');

figure;
subplot(3,4,1);
plot(0:23,Htdata_01.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_01.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
%xlabel('Local Time (hour)');
%ylabel('Effetive Scale Heights (km)');
title('Jan');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,2);
plot(0:23,Htdata_02.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_02.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
%xlabel('Local Time (hour)');
%ylabel('Effetive Scale Heights (km)');
title('Feb');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,3);
plot(0:23,Htdata_03.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_03.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
%xlabel('Local Time (hour)');
%ylabel('Effetive Scale Heights (km)');
title('Mar');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,4);
plot(0:23,Htdata_04.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_04.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
%xlabel('Local Time (hour)');
%ylabel('Effetive Scale Heights (km)');
title('Apr');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,5);
plot(0:23,Htdata_05.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_05.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
ylabel('Effetive Scale Heights (km)');
title('May');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,6);
plot(0:23,Htdata_06.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_06.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);

title('Jun');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,7);
plot(0:23,Htdata_07.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_07.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);

title('Jul');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,8);
plot(0:23,Htdata_08.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_08.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);

title('Aug');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,9);
plot(0:23,Htdata_09.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_09.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
xlabel('Local Time (hour)');

title('Sept');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,10);
plot(0:23,Htdata_10.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_10.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
xlabel('Local Time (hour)');

title('Oct');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,11);
plot(0:23,Htdata_11.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_11.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
xlabel('Local Time (hour)');

title('Nov');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;

subplot(3,4,12);
plot(0:23,Htdata_12.month_Ht_data','Color',[0.831372559070587 0.815686285495758 0.7843137383461]);hold on;
month_median_value=nanmedian(Htdata_12.month_Ht_data,1);
plot(0:23,month_median_value,'r','linewidth',2);
axis([0 23 50 250]);
set(gca, 'xtick', 0: 4: 23);
set(gca, 'ytick', 50: 50: 250);
set(gca, 'XTickLabelRotation',0);
xlabel('Local Time (hour)');

title('Dec');
set(gca,'FontSize',18, 'FontWeight','bold');
grid on;


%%%%Figures 5-7  EOF model
load('Htdata_2021_makima.mat');
%%0.
tmp = Htdata_2021_makima';
[row, col] = size(tmp);
for i = 1 : row
    Tmp(i, :) = tmp(i, :) - mean(tmp(i, :));
end

%%1.
datatmp=tmp;%(1,:);
data_cov=datatmp*datatmp';%cov(datatmp');
%%2.
[eof,eigenvalue]=eig(data_cov);
ei = sum(eigenvalue);
k = ei/sum(ei);


a=data_cov*eof;
b=eof*eigenvalue;
%%3.PC
PC=eof'*datatmp;%dot(eof,datatmp_avg);
startTime = datenum('00:00:00');
endTime = datenum('23:59:59');

me = mean(tmp, 2);

%%4.
eof_tmp=eof;
eoftem=eof_tmp(:,20:24)*PC(20:24,:);%(1,:);

%%%%%%Figure 5
figure;
set(gcf,'unit', 'normalized', 'position',[0.1,0.1,0.9,0.4]);
xData = 0:23;%startTime : datenum(0, 0, 0, 1, 0, 0) : endTime;
yyaxis left
plot(xData, eof_tmp(:,24), 'b-o', 'MarkerSize', 6,'linewidth',1.5);
ylabel('E1');
ylim([0.15 0.3]);

hold on;
yyaxis right
plot(xData, mean(tmp, 2),  'r-^', 'MarkerSize', 6,'linewidth',1.5)

ylabel('Scale Height (km)');
xlabel('Local Time (hour)');
set(gca, 'xtick', 0 : 2 : 23);
legend( 'E1','Scale Height');
legend('boxoff');
grid on;
set(gca,'FontSize',18, 'FontWeight','bold');
axis([0 23 100 180]);
title('The First Eigen Function');

%%%%%%Figure 6
figure;
set(gcf,'unit', 'normalized', 'position',[0.1,0.1,0.9,0.4]);
plot(xData, eof_tmp(:,23), 'r-o', 'MarkerSize', 6,'linewidth',1.5);
hold on;
plot(xData, eof_tmp(:,22),  'g-^', 'MarkerSize', 6,'linewidth',1.5)
hold on;
plot(xData, eof_tmp(:,21),  'b-p', 'MarkerSize', 6,'linewidth',1.5)
hold on;
plot(xData, eof_tmp(:,20),  'k-p', 'MarkerSize', 6,'linewidth',1.5)

legend('E2', 'E3', 'E4', 'E5');
legend('boxoff');
grid on;
ylabel('Base Functions');
xlabel('Local Time (hour)');
set(gca, 'xtick', 0 : 2 : 23);
set(gca,'FontSize',18, 'FontWeight','bold');
set(gca, 'xtick', 0 : 2 : 23);
set(gca, 'ytick', -0.5 : 0.25 : 0.5);
axis([0 23 -0.7 0.5]);
title('The Second, Third, Fourth and Fifth Eigen Functions');


%%%%%%Figure 7
newH = reshape(eoftem, 1, 24 * 365);
oldH = reshape(tmp, 1, 24 * 365);
x = oldH;
y = newH;
R = corrcoef(x, y);
Htdata_2021_makima_EOF=eoftem;
save('Htdata_2021_makima_EOF.mat','Htdata_2021_makima_EOF');

p = polyfit(x,y,1);
yfit = polyval(p,x);
% nyfit(nyfit > 200) = 200;
figure;
plot(x, y, 'o', 'MarkerSize', 2);
hold on;
plot(x, yfit, 'r', 'LineWidth', 2);
xlim([40 250]);
ylim([40 250]);
% xlabel('原始标高/km');
% ylabel('重构标高/km');
xlabel('Best-fit Scale Height (km)');
ylabel('Reconstructed Scale Height (km)');
set(gca,'FontSize',18, 'FontWeight','bold');

yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;

%%%%Figure 8 reconstructed TEC in October for 5min
load('ReCON_TECdata_10.mat');
Vtec=ReCON_TECdata_10.Vtec;
gTEC=ReCON_TECdata_10.gTEC;
t=ReCON_TECdata_10.t;
iv=ReCON_TECdata_10.iv;
ig=ReCON_TECdata_10.ig;
HFin=ReCON_TECdata_10.HFin;
deleHC = Vtec(iv)-gTEC(ig);
%%%%%%%%
figure;
plot(t, Vtec(iv), 'o', 'MarkerSize', 2, 'MarkerFaceColor','b');
hold on;
plot(t, gTEC(ig), 'r', 'LineWidth', 2);
hold on;
plot(t, Vtec(iv)-gTEC(ig), 'o', 'MarkerSize', 2, 'MarkerFaceColor','g');
xlim([startTime endTime]);
set(gca, 'xtick', startTime : datenum(0, 0, 0, 24, 0, 0) : endTime);
set(gcf,'unit', 'normalized', 'position',[0.1,0.1,0.9,0.4]);
% set(gca, 'xticklabel', x_label);
datetick('x', 'dd');
grid on;
legend('Epstein', 'GNSS', 'Epstein-GNSS');
legend('boxoff');
% xlabel('鏃堕棿/鏃ユ湡');
xlabel('Day');
ylabel('VTEC/TECU');
title('Octorber, 2021');
set(gca,'FontSize',18, 'FontWeight','bold');

%%%%%Figure 9 reconstructed TEC in 2021
%%%%the reconstructed data includes the density profiles, so it is large.
%%%%If interesting, please access it by contact me (chuajiang@whu.edu.cn)
% load('reconstructed_data.mat');
% time=reconstructed_data.time;
% t=reconstructed_data.t;
% iv=reconstructed_data.iv;
% ig=reconstructed_data.ig;
% Vtec=reconstructed_data.Vtec;
% f2VTEC=reconstructed_data.f2VTEC;
% iono_profiles=reconstructed_data.iono_profiles;
% 
% deleHC = Vtec(iv)-gTEC(ig);
% %%%Figure 9 histogram of error
% figure;
% edges = [-10 : 0.5 : 10];
% histogram(deleHC, edges);
% xlim([-10, 10]);
% % ylim([0, 800]);
% grid on;
% xlabel("\DeltaVTEC/TECU");
% ylabel("Number");
% title('Reconstructed TEC compared with GNSS TEC in 2021');
% set(gca,'FontSize',18, 'FontWeight','bold');
% 
% 
% length(deleHC(deleHC >= -2.5 & deleHC <= 2.5)) / length(deleHC)

%%%%Figure 10 Swarm B at 530 km compared with topside profiles
% load('SwarmB_03.mat');
% load('SwarmB_08.mat');
% load('SwarmB_12.mat');
% 
% figure;
% swarm_ne=SwarmB_03.ne;
% swarm_time=SwarmB_03.time+datenum(0,0,0,8,0,0);
% tindex=1;
% [minvalue,ipos]=min(abs(swarm_time(tindex)-reconstructed_data.time));%minvalue< 1hour
% subplot(2,2,1);
% plot(reconstructed_data.iono_profiles(ipos,:),90:0.1:90+0.1*(length(reconstructed_data.iono_profiles(ipos,:))-1),'k','linewidth',2); hold on;
% scatter(swarm_ne(tindex)*1e6,530,'*','b','linewidth',8,'Marker','pentagram');
% strtitle=sprintf('%s BJT',datestr(swarm_time(tindex)));
% title(strtitle);
% axis([0 10e11 90 800]);
% xlabel('Density (m^-^3)');
% ylabel('Altitude (km)');
% set(gca,'FontSize',14, 'FontWeight','bold');
% 
% 
% swarm_ne=SwarmB_08.ne;
% swarm_time=SwarmB_08.time+datenum(0,0,0,8,0,0);
% tindex=1;
% [minvalue,ipos]=min(abs(swarm_time(tindex)-reconstructed_data.time));%minvalue< 1hour
% subplot(2,2,2);
% plot(reconstructed_data.iono_profiles(ipos,:),90:0.1:90+0.1*(length(reconstructed_data.iono_profiles(ipos,:))-1),'k','linewidth',2); hold on;
% scatter(swarm_ne(tindex)*1e6,530,'*','b','linewidth',8,'Marker','pentagram');
% strtitle=sprintf('%s BJT',datestr(swarm_time(tindex)));
% title(strtitle);
% axis([0 4e11 90 800]);
% xlabel('Density (m^-^3)');
% ylabel('Altitude (km)');
% set(gca,'FontSize',14, 'FontWeight','bold');
% tindex=4;
% [minvalue,ipos]=min(abs(swarm_time(tindex)-reconstructed_data.time));%minvalue< 1hour
% subplot(2,2,3);
% plot(reconstructed_data.iono_profiles(ipos,:),90:0.1:90+0.1*(length(reconstructed_data.iono_profiles(ipos,:))-1),'k','linewidth',2); hold on;
% scatter(swarm_ne(tindex)*1e6,530,'*','b','linewidth',8,'Marker','pentagram');
% strtitle=sprintf('%s BJT',datestr(swarm_time(tindex)));
% title(strtitle);
% axis([0 4e11 90 800]);
% xlabel('Density (m^-^3)');
% ylabel('Altitude (km)');
% set(gca,'FontSize',14, 'FontWeight','bold');
% 
% swarm_ne=SwarmB_12.ne;
% swarm_time=SwarmB_12.time+datenum(0,0,0,8,0,0);
% tindex=2;
% [minvalue,ipos]=min(abs(swarm_time(tindex)-reconstructed_data.time));%minvalue< 1hour
% subplot(2,2,4);
% plot(reconstructed_data.iono_profiles(ipos,:),90:0.1:90+0.1*(length(reconstructed_data.iono_profiles(ipos,:))-1),'k','linewidth',2); hold on;
% scatter(swarm_ne(tindex)*1e6,530,'*','b','linewidth',8,'Marker','pentagram');
% strtitle=sprintf('%s BJT',datestr(swarm_time(tindex)));
% title(strtitle);
% axis([0 10e11 90 800]);
% xlabel('Density (m^-^3)');
% ylabel('Altitude (km)');
% set(gca,'FontSize',14, 'FontWeight','bold');


