function [] = dvPlots(ID)

% variation with dz ----
% dvPlots('??2106032') %0.4 atm; center
% dvPlots('??2206032') %0.4 atm; center; seed2
% dvPlots('??2107132') %0.4 atm; UR; 
% dvPlots('??3107132') %0.6 atm; UR;
% dvPlots('??5107132') %1.0 atm; UR;
% dvPlots('??7107132') %1.4 atm; UR;
% dvPlots('??3117132') %0.6 atm; UR; M2 x decenter 0.5mm
% dvPlots('??3127132') %0.6 atm; UR; M2 x-tilt 0.01 deg
% dvPlots('??3117142') %0.6 atm; UR; M2 x decenter 0.5mm; 150s exp;
% dvPlots('??3102132') %0.6 atm; UR; r-band
% dvPlots('??3112132') %0.6 atm; UR; M2 x decenter 0.5mm; r-band
% dvPlots('??3122132') %0.6 atm; UR; M2 x tilt 0.01 deg; r-band
% dvPlots('??3117122') %0.6 atm; UR; M2 x decenter 0.5mm; 10s exp;
% dvPlots('??3117112') %0.6 atm; UR; M2 x decenter 0.5mm; 1s exp;

% variation with atm ---
% dvPlots('10?107132') % 1.0mm; design; 770nm; UR; 15s; default CCD
% dvPlots('20?107132') % 2.0mm; design; 770nm; UR; 15s; default CCD

% variation with system state --
% dvPlots('1031?7132') % 1.0mm; 770nm/ UR; 15s; default CCD
% dvPlots('2031?7132') % 2.0mm; 770nm/ UR; 15s; default CCD
% dvPlots('1031?2132') % 1.0mm; r-band / UR; 15s; default CCD
% dvPlots('2031?2132) % 2.0mm; r-band / UR; 15s; default CCD

% variation with exposure time ---
% dvPlots('1031171?2') % 1.0mm; M2 decenter; 770nm; UR; default CCD

znmax = 22;
if strcmp(ID(1:2),'??') %dz
    dv = [0.5 1.0 1.5 2.0 2.5];
    % dv = [0.5 1.0 1.5 2.0 ];
    % dv = [0.5 1.0 2.0 ];
    % dv = [1.0 2.0 ];
    ndv = length(dv);
    dvString=cell(ndv,1);
    for idv=1:length(dv),dvString{idv}=sprintf('%02d',dv(idv) *10);end
    xL = 'sensor offset (mm)';
elseif strcmp(ID(3),'?') % atm
    dv = [0 0.4 0.6 1.0 1.4];
    ndv = length(dv);
    dvString=cell(ndv,1);
    for idv=1:length(dv),dvString{idv}=sprintf('%d',int8(dv(idv)/0.2));end    
    xL = 'seeing (outer scale included) (arcsec)';
elseif strcmp(ID(5),'?') % optics
    dv = [0 1 2];
    ndv = length(dv);
    dvString=cell(ndv,1);
    for idv=1:length(dv),dvString{idv}=sprintf('%d',int8(dv(idv)));end    
    xL = 'design;M2 decenter; M2 tilt';
elseif strcmp(ID(8),'?') % exp time
    dv = [1 10 15 150];
    ndv = length(dv);
    dvString=cell(ndv,1);
    for idv=1:length(dv),dvString{idv}=sprintf('%d',int8(idv));end    
    xL = 'exposure time (second)';
end

znmax3 = znmax-3;
zMean = zeros(ndv, znmax3);
zDev = zeros(ndv, znmax3);
for idv = 1:ndv
    IDString = strrep(ID,repmat('?',1,length(dvString{idv})), dvString{idv});
    sumTXT = sprintf('output/wfs_%s_100_sum.txt', IDString);
    txtData = load(sumTXT);
    if idv == 1
        zTrue = txtData(1,:);
    end
    zMean(idv, :) = txtData(2,:);
    zDev(idv, :) = txtData(3,:);
end

figure(1);clf;
subplot(2,2,1);
dd = zMean-repmat(zTrue,ndv,1);
[~, idx] = sort(rms(dd),'descend');
plot(dv,dd(:,idx)','.-','linewidth',2); 
idx = idx+3;
legend(['Z' num2str(idx(1))],['Z' num2str(idx(2))],['Z' num2str(idx(3))],...
    ['Z' num2str(idx(4))],'Location','southoutside');
xlabel(xL);
ylabel('Measured - Truth (nm)');
title('(Measured - Truth) for each Zernike');
% y1=ylim;

subplot(4,2,5);
plot(dv, sqrt(sum((zMean-repmat(zTrue,ndv,1)).^2,2)),'-ro','linewidth',2);
xlabel(xL);
ylabel('Measured - Truth (nm)');
title('(Measured - Truth) RSSed');
grid on;
% ylim([0 y1(2)-y1(1)]);

subplot(2,2,2);
[~, idx] = sort(rms(zDev),'descend');
plot(dv,zDev(:,idx)','.-','linewidth',2);
idx = idx+3;
legend(['Z' num2str(idx(1))],['Z' num2str(idx(2))],['Z' num2str(idx(3))],...
    ['Z' num2str(idx(4))],'Location','southoutside');
xlabel(xL);
ylabel('STD (nm)');
title('STD for each Zernike');

subplot(4,2,6);
if strcmp(ID(8),'?')
    loglog(dv, sqrt(sum(zDev.^2,2)),'-ro','linewidth',2);    
else
    plot(dv, sqrt(sum(zDev.^2,2)),'-ro','linewidth',2);
end
xlabel(xL);
ylabel('STD (nm)');
title('STD RSSed');
grid on;

end
