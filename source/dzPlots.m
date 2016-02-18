function [] = dzPlots(fixed)

% dzPlots('2106032') %0.4 atm; center
% dzPlots('2206032') %0.4 atm; center; seed2
% dzPlots('2107132') %0.4 atm; UR; 
% dzPlots('5107132') %1.0 atm; UR;

znmax = 22;
dz = [0.5 1.0 1.5 2.0 2.5];
dz = [0.5 1.0 1.5 2.0 ];
% dz = [0.5 1.0 2.0 ];

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
for iz=4:znmax
    plot(dz,zMean(:, iz-3)-zTrue(iz-3)); hold on;
end
xlabel('sensor offset (mm)');
ylabel('Measured - Truth (nm)');
title('(Measured - Truth) for each Zernike');

subplot(2,2,3);
plot(dz, sqrt(sum((zMean-repmat(zTrue,ndz,1)).^2,2)),'-ro');
xlabel('sensor offset (mm)');
ylabel('Measured - Truth (nm)');
title('(Measured - Truth) RSSed');
grid on;

subplot(2,2,2);
for iz=4:znmax
    plot(dz,zDev(:, iz-3)); hold on;
end
xlabel('sensor offset (mm)');
ylabel('STD (nm)');
title('STD for each Zernike');

subplot(2,2,4);
plot(dz, sqrt(sum(zDev.^2,2)),'-ro');
xlabel('sensor offset (mm)');
ylabel('STD (nm)');
title('STD RSSed');
grid on;

end
