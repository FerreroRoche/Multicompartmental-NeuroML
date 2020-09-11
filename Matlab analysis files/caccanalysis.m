
%// Define the x values
x = t;
xMat = repmat(t, 1, 21); %// For plot3

%// Define y values
y = linspace(0,21,21);
yMat = repmat(y, numel(x), 1); %//For plot3

zMat = [Baseline cacc05 cacc10 cacc15 cacc20 cacc25 cacc30 cacc35 cacc40 cacc45 cacc50 cacc55 cacc60 cacc65 cacc70 cacc75 cacc80 cacc85 cacc90 cacc95 cacc100];



figure(2)
plot3(xMat,yMat,zMat)
grid;
xlabel('t(s)'); ylabel('mV'); zlabel('');
view(40,40); 

figure(3)
waterfall(xMat,yMat,zMat)


%% 

for i = 1:length(zMat(1,:))

    %figure
    [pks, locs] = findpeaks(zMat(:,i),t);
    peakInterval = diff(locs);
    %hist(peakInterval)
    %grid on
    %xlabel('second Intervals')
    %ylabel('Frequency of Occurrence')
    %title('Histogram of Peak Intervals (seconds)')
    AverageDistance_Peaks = mean(diff(locs));
    Averagefrequency(i) = 1/AverageDistance_Peaks;
    
    
end

density = linspace(0,100,21);
plot(density(),Averagefrequency())
grid on
xlabel('CaCC Density (mS_per_cm2)')
ylabel('Oscillatory frequency (Hz)')
title('CaCC infulence on subthreshold oscilaltion frequency')
%%


h = surf(xMat,yMat,zMat);
set(h,'LineStyle','none')



%%


%% Cacc analysis why is it so powerfull???
figure(2)
fs = 40000;                    % Sampling frequency (samples per second)
dt = 1/fs;                   % seconds per sample
StopTime = 3;             % seconds
tsin = (0:dt:StopTime-dt)';     % seconds
F = 9.1591;                      % Sine wave frequency (hertz)
data = sin(2*pi*F*tsin);
plot(data)


%%

gcl = 6.5*10^-9; % conductance
Ecl = -45*10^-3; % volt
Clh = 0.37*10^-6; %%microMolar
caconc = 30*data;

for i = 1:length(Baseline)-1
    
    m = 1 / (1 + exp((Clh - caconc(i))/0.09));

    I_Cacc(i) = gcl * m *(Baseline(i) - Ecl);
end
figure(1)
plot(t,I_Cacc)
figure(2)
plot(t,Baseline)
