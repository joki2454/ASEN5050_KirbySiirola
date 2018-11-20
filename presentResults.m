%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/17/2018
%
% Purpose:  This function presents various results for the project in the
% form of plots.
%
% Inputs:
%   nominal - struct containing all data from the nominal (DVpJ = [0 0 0]'
%             km/s) case
%   DVpJ    - SET.PRESENT.numSteps x SET.PRESENT.numSteps x 3 matrix of
%             DeltaV vectors in J2000 performed at perijove to attain a
%             desired sma and inc (each corresponding to the first and second
%             dimension of the matrix, respectively), km/s
%   SET     - struct of settings, initial conditions, and options
%
function [] = presentResults(nominal,DVpJ,SET)
%% Allocation
Jcolor = [255,165,0]./256;

%% Calculate sma and inc vectors and produce meshgrid
sma = linspace(SET.RANGES.sma(1),SET.RANGES.sma(2),SET.PRESENT.numSteps); % sma
inc = linspace(SET.RANGES.inc(1),SET.RANGES.inc(2),SET.PRESENT.numSteps); % deg
[SMA,INC] = meshgrid(sma,inc);
SMA = SMA';
INC = INC';

%% Calculate transfer parameters over entire range
cm = 1; % mesh counter
for i = 1:size(DVpJ,1)
  for j = 1:size(DVpJ,2)
    [KE_PARK(i,j),TOF(i,j),T0(i,j),JSOI1(i,j),PJ(i,j),JSOI2(i,j),SSOI(i,j),PS(i,j),optim_badness] = ...
          transferSequence(squeeze(DVpJ(i,j,:)),SET);
    % Make sure all solutions actually enter Saturn's SOI, otherwise there
    %   is a serious problem
    if optim_badness
      error(['One of the solutions for DVpJ produces an orbit which does not ',...
              'intersect Saturn''s sphere-of-influence, which is not handled by ',...
              'the code as-is.  This most likely occurred because the upper range ',...
              'for target semi-major axes (SET.RANGES.sma(2) in projectInitialize.m) ',...
              'is too large, try reducing it to something ',...
              'significantly less than Saturn''s sphere-of-influence'])
    end
  end
end

%% Calculate all |DeltaV|s and residuals
for i = 1:length(sma)
  for j = 1:length(inc)
    res(i,j)    = sqrt(((KE_PARK(i,j).a-sma(i))/sma(i))^2 + ((KE_PARK(i,j).i-inc(j))/inc(j))^2); % %
    Dv_pJ(i,j)  = norm(squeeze(DVpJ(i,j,:)));                    % km/s
    Dv_pS(i,j)  = norm(PS(i,j).DV);                              % km/s
    Dv_tot(i,j) = norm(squeeze(DVpJ(i,j,:))) + norm(PS(i,j).DV); % km/s
  end
end

%% Calculate trajectories for varied parking orbit conditions
sma_ind = round(linspace(1,length(inc),2));
inc_ind = round(linspace(1,length(inc),5));
for i = 1:length(sma_ind)
  for j = 1:length(inc_ind)
    [~,R_hc(i,j).data,R_jc(i,j).data,R_sc(i,j).data,P_EPHEM(i,j)] = transferStates(TOF(sma_ind(i),inc_ind(j)),...
      T0(sma_ind(i),inc_ind(j)),JSOI1(sma_ind(i),inc_ind(j)),PJ(sma_ind(i),inc_ind(j)),...
      JSOI2(sma_ind(i),inc_ind(j)),SSOI(sma_ind(i),inc_ind(j)),...
      PS(sma_ind(i),inc_ind(j)),KE_PARK(sma_ind(i),inc_ind(j)),SET);
  end
end


%% |DeltaV| at perijove
figure('name','Dv_pJ')
hold on
grid on
grid minor
surf(SMA,INC,Dv_pJ,'edgealpha',0.2)
plot3(nominal.KE_PARK.a,nominal.KE_PARK.i,norm(nominal.DVpJ),'r.','markersize',15)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{pJ} (km/s)')
title('|\DeltaV|_{perijove} vs. SMA and INC')
hold off

%% |DeltaV| at perisaturnium
lims = [min(min(Dv_pS)) Dv_pS(end,end)];
figure('name','Dv_pS','units','normalized','position',[0.125 0.4 0.75 0.5])
subplot(1,2,1)
hold on
grid on
grid minor
surf(SMA,INC,Dv_pS,'edgealpha',0.2)
plot3(nominal.KE_PARK.a,nominal.KE_PARK.i,norm(nominal.PS.DV),'r.','markersize',15)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{pS} (km/s)')
caxis(lims)
title('|\DeltaV|_{perisaturnium} vs. SMA and INC')
hold off

subplot(1,2,2)
hold on
grid on
grid minor
surf(SMA,INC,Dv_pS,'edgealpha',0.2)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{pS} (km/s)')
plot3(nominal.KE_PARK.a,nominal.KE_PARK.i,norm(nominal.PS.DV),'r.','markersize',15)
zlim(lims)
caxis(lims)
title('Zoomed |\DeltaV|_{perisaturnium} vs. SMA and INC')
hold off

%% |DeltaV| total
lims = [min(min(Dv_tot)) Dv_tot(end,end)];
figure('name','Dv_tot','units','normalized','position',[0.125 0.4 0.75 0.5])
subplot(1,2,1)
hold on
grid on
grid minor
surf(SMA,INC,Dv_tot,'edgealpha',0.2)
plot3(nominal.KE_PARK.a,nominal.KE_PARK.i,norm(nominal.PS.DV),'r.','markersize',15)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{total} (km/s)')
caxis(lims)
title('|\DeltaV|_{tot} vs. SMA and INC')
hold off

subplot(1,2,2)
hold on
grid on
grid minor
surf(SMA,INC,Dv_tot,'edgealpha',0.2)
plot3(nominal.KE_PARK.a,nominal.KE_PARK.i,norm(nominal.PS.DV),'r.','markersize',15)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{total} (km/s)')
zlim(lims)
caxis(lims)
title('Zoomed |\DeltaV|_{tot} vs. SMA and INC')
hold off

%% Sum of squared residuals
figure('name','residuals')
hold on
grid on
grid minor
surf(SMA,INC,res,'edgealpha',0.2)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('(\deltaSMA^2+\deltaINC^2)^{1/2}')
title('Sum of Squared SMA/INC Residuals vs. SMA and INC')
hold off

%% RAAN
figure('name','RAAN')
hold on
grid on
grid minor
for i = 1:size(KE_PARK,2)
  raan(:,i) = extractfield(KE_PARK(:,i),'raan');
end
surf(SMA,INC,raan,'edgealpha',0.2)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('RAAN (deg)')
title('RAAN vs. SMA and INC')
hold off

%% DVpJ Solution Space
% Make color data for DVpJ
for i = 1:length(sma)
  for j = 1:length(inc)
    CDataInc(i,j) = inc(j);
    CDataSma(i,j) = sma(i);
  end
end

figure('name','DVpJ Soln Space','units','normalized','position',[0.125 0.125 0.75 0.75])
subplot(1,2,1)
hold on
grid on
grid minor
h = surf(DVpJ(:,:,1),DVpJ(:,:,2),DVpJ(:,:,3),CDataSma,'edgealpha',0.2);
h.EdgeAlpha = 0.3;
xlabel('Vx (km/s)')
ylabel('Vy (km/s)')
zlabel('Vz (km/s)')
axis equal
title('$\Delta\vec{V}_{perijove}$ for Various Final SMAs','interpreter','latex')
c = colorbar('location','southoutside');
set(get(c,'Label'),'String','Final SMA (deg)');
hold off

subplot(1,2,2)
hold on
grid on
grid minor
h = surf(DVpJ(:,:,1),DVpJ(:,:,2),DVpJ(:,:,3),CDataInc,'edgealpha',0.2);
h.EdgeAlpha = 0.3;
xlabel('Vx (km/s)')
ylabel('Vy (km/s)')
zlabel('Vz (km/s)')
axis equal
title('$\Delta\vec{V}_{perijove}$ for Various Final INCs','interpreter','latex')
c = colorbar('location','southoutside');
set(get(c,'Label'),'String','Final INC (deg)');
hold off


%% Nominal (DVpJ = 0) heliocentric trajectory
[~,nominal.R_hc,nominal.R_jc,nominal.R_sc,nominal.P_EPHEM] = transferStates(nominal.TOF,nominal.T0,nominal.JSOI1,...
  nominal.PJ,nominal.JSOI2,nominal.SSOI,nominal.PS,nominal.KE_PARK,SET);
figure('name','3D Heliocentric','units','normalized','position',[0.125 0.125,0.75 0.75])
hold on
grid on
c = get(gca,'colororder');
% plot trajectory
plot3(nominal.R_hc(1,:),nominal.R_hc(2,:),nominal.R_hc(3,:),'color',c(1,:))
% plot jupiter (not to scale) and jupiter trajectory
plotBody3D(nominal.P_EPHEM.RJ_pJ,150*SET.CONST.RJ,Jcolor);
plot3(nominal.P_EPHEM.RJ_hc(1,:),nominal.P_EPHEM.RJ_hc(2,:),nominal.P_EPHEM.RJ_hc(3,:),'color',c(2,:))
% plot saturn (not to scale) and saturn trajectory
plotBody3D(nominal.P_EPHEM.RS_pS,150*SET.CONST.RS,[1 1 1]);
plot3(nominal.P_EPHEM.RS_hc(1,:),nominal.P_EPHEM.RS_hc(2,:),nominal.P_EPHEM.RS_hc(3,:),'color',c(3,:))
% plot sun (not to scale)
plotBody3D([0 0 0]',50*SET.CONST.RSun,[1 1 0]);
set(gca,'color',[0 0 0],'gridcolor',[0.9 0.9 0.9],'minorgridcolor',[0.9 0.9 0.9],'gridalpha',0.4,'minorgridalpha',0.4)
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('Spacecraft (|\DeltaV_{perijove}|=0 km/s) in Heliocentric J2000')
h = legend('Cassini Trajectory','Jupiter (150x Scale)',...
  'Jupiter Trajectory','Saturn (150x Scale)','Saturn Trajectory',...
  'Sun (50x Scale)','location','westoutside');
h.Color = [0.8 0.8 0.8];
hold off

%% Nominal jupiter-centric and saturn-centric trajectories
figure('name','Nominal Planet-Centric','units','normalized','position',[0.125 0.125 0.75 0.75])
subplot(1,2,1)
hold on
grid on
plot3(nominal.R_jc(1,:),nominal.R_jc(2,:),nominal.R_jc(3,:),'color',c(1,:))
plotBody3D([0 0 0]',10*SET.CONST.RJ,Jcolor);
[x,y,z] = sphere;
h = surf(SET.CONST.JSOI*x,SET.CONST.JSOI*y,SET.CONST.JSOI*z,ones(size(x,1),size(x,2),3));
h.FaceAlpha = 0.1;
h.EdgeAlpha = 0.3;
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
set(gca,'color',[0 0 0],'gridcolor',[0.9 0.9 0.9],'minorgridcolor',[0.9 0.9 0.9],'gridalpha',0.4,'minorgridalpha',0.4)
title({'Jupiter Flyby in Jupiter-Centric J2000','Nominal (|\DeltaV_{perijove}|=0 km/s) Case'})
h = legend('Cassini Trajectory','Jupiter (10x Scale)','Jupiter SOI','location','southoutside');
h.Color = [0.8 0.8 0.8];
hold off

subplot(1,2,2)
hold on
grid on
plot3(nominal.R_sc(1,:),nominal.R_sc(2,:),nominal.R_sc(3,:),'color',c(1,:))
plotBody3D([0 0 0]',10*SET.CONST.RS,[1 1 1]);
[x,y,z] = sphere;
h = surf(SET.CONST.SSOI*x,SET.CONST.SSOI*y,SET.CONST.SSOI*z,ones(size(x,1),size(x,2),3));
h.FaceAlpha = 0.1;
h.EdgeAlpha = 0.3;
axis equal
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
set(gca,'color',[0 0 0],'gridcolor',[0.9 0.9 0.9],'minorgridcolor',[0.9 0.9 0.9],'gridalpha',0.4,'minorgridalpha',0.4)
title({'Saturn Flyby in Saturn-Centric J2000','Nominal (|\DeltaV_{perijove}|=0 km/s) Case'})
h = legend('Cassini Trajectory','Saturn (10x Scale)','Saturn SOI','location','southoutside');
h.Color = [0.8 0.8 0.8];
hold off

%% Jupiter centric trajectories
figure('name','3D Jupiter-Centric','units','normalized','position',[0.125 0.125 0.75 0.75])
% plot Cassini trajectories
%plot3(nominal.R_sc(1,:),nominal.R_sc(2,:),nominal.R_sc(3,:),'color',c(1,:))
for i = 1:length(sma_ind)
  subplot(1,length(sma_ind),i)
  hold on
  grid on
  lc = 1; % legend counter
  for j = 1:length(inc_ind)
    p(lc) = plot3(R_jc(i,j).data(1,:),R_jc(i,j).data(2,:),R_jc(i,j).data(3,:));
    legend_str{lc} = ['Final INC = ',num2str(round(inc(inc_ind(j)),2)),' deg'];
    lc = lc + 1;
  end
  axis equal
  xlabel('x (km)')
  ylabel('y (km)')
  zlabel('z (km)')
  % plot jupiter (10*scale)
  p(lc) = plotBody3D([0 0 0]',10*SET.CONST.RJ,Jcolor);
  legend_str{lc} = 'Jupiter (10x Scale)';
  h = legend(p,legend_str,'location','southoutside');
  h.Color = [0.8 0.8 0.8];
  set(gca,'color',[0 0 0],'gridcolor',[0.9 0.9 0.9],'minorgridcolor',[0.9 0.9 0.9],'gridalpha',0.4,'minorgridalpha',0.4)
  title({'Jupiter Flyby in Jupiter-Centric J2000',['(Final SMA = ',num2str(round(sma(sma_ind(i)),2)),' km) ']})
  hold off
end




%% Saturn centric trajectories
clear p h legend_str
figure('name','3D Saturn-Centric','units','normalized','position',[0.125 0.125 0.75 0.75])
% plot Cassini trajectories
%plot3(nominal.R_sc(1,:),nominal.R_sc(2,:),nominal.R_sc(3,:),'color',c(1,:))
for i = 1:length(sma_ind)
  subplot(1,length(sma_ind),i)
  hold on
  grid on
  lc = 1; % legend counter
  for j = 1:length(inc_ind)
    p(lc) = plot3(R_sc(i,j).data(1,:),R_sc(i,j).data(2,:),R_sc(i,j).data(3,:));
    legend_str{lc} = ['Final INC = ',num2str(round(inc(inc_ind(j)),2)),' deg'];
    lc = lc + 1;
  end
  axis equal
  xlabel('x (km)')
  ylabel('y (km)')
  zlabel('z (km)')
  if i == 1
    % plot saturn (1*scale)
    p(lc) = plotBody3D([0 0 0]',1*SET.CONST.RS,[1 1 1]);
    legend_str{lc} = 'Saturn (1x Scale)';
    lims = [-1 1].*5*SET.RANGES.sma(1);
    h = legend(p,legend_str,'location','southoutside');
    h.Color = [0.8 0.8 0.8];
  else
    % plot saturn (10*scale)
    p(lc) = plotBody3D([0 0 0]',10*SET.CONST.RS,[1 1 1]);
    legend_str{lc} = 'Saturn (10x Scale)';
    lims = [-1 1].*2*SET.RANGES.sma(2);
    h = legend(p,legend_str,'location','southoutside');
    h.Color = [0.8 0.8 0.8];
  end
  xlim(lims)
  ylim(lims)
  zlim(lims)
  set(gca,'color',[0 0 0],'gridcolor',[0.9 0.9 0.9],'minorgridcolor',[0.9 0.9 0.9],'gridalpha',0.4,'minorgridalpha',0.4)
  title({'Saturn Arrival in Saturn-Centric J2000',['(Final SMA = ',num2str(round(sma(sma_ind(i)),2)),' km) ']})
  hold off
end





end