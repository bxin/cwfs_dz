function [] = dzPlots(fixed)

% dzPlots('2106032') %0.4 atm; center
% dzPlots('2206032') %0.4 atm; center; seed2
% dzPlots('2107132') %0.4 atm; UR; 
% dzPlots('3107132') %0.6 atm; UR;
% dzPlots('5107132') %1.0 atm; UR;
% dzPlots('7107132') %1.4 atm; UR;
% dzPlots('3117132') %0.6 atm; UR; M2 x decenter 0.5mm
% dzPlots('3127132') %0.6 atm; UR; M2 x-tilt 0.01 deg
% dzPlots('3117142') %0.6 atm; UR; M2 x decenter 0.5mm; 150s exp;

% dzPlots('3102132') %0.6 atm; UR; r-band
% dzPlots('3112132') %0.6 atm; UR; M2 x decenter 0.5mm; r-band


znmax = 22;
dz = [0.5 1.0 1.5 2.0 2.5];
% dz = [0.5 1.0 1.5 2.0 ];
% dz = [0.5 1.0 2.0 ];
% dz = [1.0 2.0 ];

ndz = length(dz);
znmax3 = znmax-3;
zMean = zeros(ndz, znmax3);
zDev = zeros(ndz, znmax3);
for idz = 1:ndz
    sumTXT = sprintf('output/wfs_%02d%s_100_sum.txt', dz(idz)*10, fixed);
    txtData = load(sumTXT);
    if idz == 1
        zTrue = txtData(1,:);
    end
    zMean(idz, :) = txtData(2,:);
    zDev(idz, :) = txtData(3,:);
end

figure(1);clf;
subplot(2,2,1);
dd = zMean-repmat(zTrue,ndz,1);
[~, idx] = sort(rms(dd),'descend');
plot(dz,dd(:,idx)','.-','linewidth',2); 
idx = idx+3;
legend(['Z' num2str(idx(1))],['Z' num2str(idx(2))],['Z' num2str(idx(3))],...
    ['Z' num2str(idx(4))],'Location','southoutside');
xlabel('sensor offset (mm)');
ylabel('Measured - Truth (nm)');
title('(Measured - Truth) for each Zernike');
% y1=ylim;

subplot(4,2,5);
plot(dz, sqrt(sum((zMean-repmat(zTrue,ndz,1)).^2,2)),'-ro','linewidth',2);
xlabel('sensor offset (mm)');
ylabel('Measured - Truth (nm)');
title('(Measured - Truth) RSSed');
grid on;
% ylim([0 y1(2)-y1(1)]);

subplot(2,2,2);
[~, idx] = sort(rms(zDev),'descend');
plot(dz,zDev(:,idx)','.-','linewidth',2);
idx = idx+3;
legend(['Z' num2str(idx(1))],['Z' num2str(idx(2))],['Z' num2str(idx(3))],...
    ['Z' num2str(idx(4))],'Location','southoutside');
xlabel('sensor offset (mm)');
ylabel('STD (nm)');
title('STD for each Zernike');

subplot(4,2,6);
plot(dz, sqrt(sum(zDev.^2,2)),'-ro','linewidth',2);
xlabel('sensor offset (mm)');
ylabel('STD (nm)');
title('STD RSSed');
grid on;

end
