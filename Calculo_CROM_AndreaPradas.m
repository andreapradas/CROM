clear
close all
clc
%% Data reading from the file of Patients
Patient_ID = 10;

opts = delimitedTextImportOptions("NumVariables", 13);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

opts.VariableNames = ["Var1", "Time", "Device", "X", "Y", "Z", "qW", "qX", "qY", "qZ", "aX", "aY", "aZ"];
opts.SelectedVariableNames = ["Time", "Device", "X", "Y", "Z", "qW", "qX", "qY", "qZ", "aX", "aY", "aZ"];
opts.VariableTypes = ["string", "datetime", "categorical", "double", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double"];

opts.ExtraColumnsRule = "addvars";
opts.EmptyLineRule = "read";

opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Device"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Time", "InputFormat", "HH:mm:ss.SSSS");

fileName = "CROM RawData\Paciente_"+Patient_ID+".txt";

Patient = readtable(fileName, opts);

clear opts prompt

%% Concatenate the values of additional columns

[~,n]=size(Patient);

if n > 12

    Patient.aX = strcat(Patient.ExtraVar5,'.',Patient.ExtraVar6);
    Patient.aY = strcat(Patient.ExtraVar7,'.',Patient.ExtraVar8);
    Patient.aZ = strcat(Patient.ExtraVar9,'.',Patient.ExtraVar10);
end

%% Filtrar datos de la CABEZA

PteH = Patient(Patient.Device=="H",10:12);% Guarda las 3 ultimas col de Paciente de las medidas de la CABEZA

if n > 12
    PteH.aX = str2double(PteH.aX);
    PteH.aY = str2double(PteH.aY);
    PteH.aZ = str2double(PteH.aZ);
end

clear fileName i m n

% Angle between -180 and 180
PteH.aX = wrapTo180(PteH.aX);
PteH.aY = wrapTo180(PteH.aY);
PteH.aZ = wrapTo180(PteH.aZ);

%% FLEXION-EXTENSION

% Delete the OFFSET
mean_aX = mean(PteH.aX(1:9));
PteH_aX_NoOffset = PteH.aX - mean_aX;

% Filter the signal
Freq = 90; % Cuantas muestras se cogen por unidad de tiempo (es una funcion discretizada)
fc = 10; % 
[b,a] = butter(4,fc/(Freq/2));
PteH_aX_Filtered = filtfilt(b,a,PteH_aX_NoOffset);

figure(1);
subplot(1,2,1);
plot(PteH_aX_Filtered);
yline(0, 'm--', 'LineWidth', 1);
title("PteH aX Filtered");

% Obtain the movement complete (warm up included)
[~, locs_MAX_F] = findpeaks(PteH_aX_Filtered,'MinPeakHeight',30,'MinPeakDistance',100,'SortStr','ascend');
aX_inv = -PteH_aX_Filtered;
[~,locs_MIN_F] = findpeaks(aX_inv,'MinPeakHeight',25,'MinPeakDistance',100,'SortStr','ascend');
peak_locs = sort([locs_MAX_F(:); locs_MIN_F(:)]); % Las pos en el eje X de todos los picos
% Find out if the mov. is correctly done
mean_durations = zeros(1, 6);
for i = 1:6
    mean_durations(i) = round(peak_locs(2*i) - peak_locs(2*i - 1));
end
if length(peak_locs) == 12
    start_index = peak_locs(1) - mean_durations(1);
    end_index = peak_locs(12) + mean_durations(6);
    PteH_aX_Mov = PteH_aX_Filtered(start_index:end_index);
else 
    disp("The FLEXO-EXTENSION has been done incorrectly by the participant.");
    return;
end

subplot(1,2,2);
plot(PteH_aX_Mov);
yline(0, 'm--', 'LineWidth', 1);
title("Flexo-Extension WITH WARM UP");

clear i start_index end_index locs_MIN_F locs_MAX_F
%% Separate Flexion-Extension (without warm up)

first_peak = peak_locs(7); % The first peak of useful mov.
PteH_aX_Flexion = PteH_aX_Filtered(first_peak-mean_durations(4):peak_locs(12)+mean_durations(6));

figure(2);
plot(PteH_aX_Flexion);
yline(0, 'm--', 'LineWidth', 1);
title("FLEXO-EXTENSION without warm up");

% The signal without the rest time position
min_segment_length = 150;
contiguous_segments = splitConsecutives(find(abs(diff(PteH_aX_Flexion)) < 0.1));% Buscar zonas sin mucha variacion
rest_position = contiguous_segments(cellfun(@length, contiguous_segments) > min_segment_length);

FERep1 = PteH_aX_Flexion(1:rest_position{1}(1));
FERep2 = PteH_aX_Flexion(rest_position{1}(end):rest_position{end}(1));
FERep3 = PteH_aX_Flexion(rest_position{end}(end):end);

PteH_aX = zeros(length(FERep1)+length(FERep2)+length(FERep3),1);
PteH_aX(1:length(FERep1),1) = FERep1;
PteH_aX(length(FERep1)+1:length(FERep1)+length(FERep2),1) = FERep2;
PteH_aX(length(FERep1)+length(FERep2)+1:length(FERep1)+length(FERep2)+length(FERep3),1) = FERep3;

figure(3);
plot(PteH_aX);
yline(0, 'm--', 'LineWidth', 1);
title("Repetions of Flex-Ext");

clear peak_locs first_peak last_peak mean_durations contiguous_segments rest_position

%% Calculate the CROM FE
[~,locs_MAX_F] = findpeaks(PteH_aX,'MinPeakHeight',30,'MinPeakDistance',20);
aX_inv = -PteH_aX;
[~,locs_MIN_F] = findpeaks(aX_inv,'MinPeakHeight',30,'MinPeakDistance',20);
 Flexion = (PteH_aX(locs_MAX_F)-PteH_aX(locs_MIN_F));
 Flexion = transpose(Flexion);
 mF = mean(Flexion);

 clear locs_MIN_F locs_MAX_F

%% Calculate the Velocity for each repetition FE

% Calculate the time elapsed for each repetition
tF1 = (linspace(1, length(FERep1), length(FERep1))/Freq)';
tF2 = (linspace(1, length(FERep2), length(FERep2))/Freq)';
tF3 = (linspace(1, length(FERep3), length(FERep3))/Freq)';

% Calculate the time increment, how many time has passsed between samples
dtF1 = gradient(tF1);
dtF2 = gradient(tF2);
dtF3 = gradient(tF3);

% Derivates --> Instant velocity
dFERep1 = gradient(FERep1(:));
dFERep2 = gradient(FERep2(:));
dFERep3 = gradient(FERep3(:));

% Velocity
vF1 = dFERep1 ./ dtF1;
vF2 = dFERep2 ./ dtF2;
vF3 = dFERep3 ./ dtF3;

% Filtering
fcV = 4;
[bV,aV] = butter(4,fcV/(Freq/2));
vF1 = filtfilt(bV,aV,vF1);
vF2 = filtfilt(bV,aV,vF2);
vF3 = filtfilt(bV,aV,vF3);

% Concatenate in a single vector
vF = zeros(length(vF1)+length(vF2)+length(vF3),1);
vF(1:length(vF1),1) = vF1;
vF(length(vF1)+1:length(vF1)+length(vF2),1) = vF2;
vF(length(vF1)+length(vF2)+1:length(vF1)+length(vF2)+length(vF3),1) = vF3;

figure(4);
plot(vF);
yline(0, 'm--', 'LineWidth', 1);
title("Velocity of Flex-Ext Movement");

% Calculate the max and min MEAN velocity
vF1Max = max(vF1);
vF2Max = max(vF2);
vF3Max = max(vF3);
vFMaxMean = mean([vF1Max;vF2Max;vF3Max]);
vFMax = [vF1Max,vF2Max,vF3Max,vFMaxMean];

vF1Min = min(vF1);
vF2Min = min(vF2);
vF3Min = min(vF3);
vFMinMean = mean([vF1Min,vF2Min,vF3Min]);
vFMin = [vF1Min,vF2Min,vF3Min,vFMinMean];

% Calculate the middle value of the peak velocity
vFabs = abs(vF);
MPDv = int16(length(vF)/9*0.6);
[~,locs_MAX_Fv] = findpeaks(vFabs,'MinPeakHeight',30,'MinPeakDistance',MPDv);
vFm = mean(vFabs(locs_MAX_Fv));

%% Acceleration calculation FE

% Derivates of velocity to obtain acceleration
aFERep1 = gradient(vF1(:));
aFERep2 = gradient(vF2(:));
aFERep3 = gradient(vF3(:));

% Acceleration
aF1 = aFERep1 ./ dtF1;
aF2 = aFERep2 ./ dtF2;
aF3 = aFERep3 ./ dtF3;

% Filtering 
fcA = 4; 
[bA, aA] = butter(4, fcA / (Freq/2));
aF1 = filtfilt(bA, aA, aF1);
aF2 = filtfilt(bA, aA, aF2);
aF3 = filtfilt(bA, aA, aF3);

% Concatenate in a single vector
aF = zeros(length(aF1)+length(aF2)+length(aF3),1);
aF(1:length(aF1),1) = aF1;
aF(length(aF1)+1:length(aF1)+length(aF2),1) = aF2;
aF(length(aF1)+length(aF2)+1:length(aF1)+length(aF2)+length(aF3),1) = aF3;

figure(5);
plot(aF);
yline(0, 'm--', 'LineWidth', 1);
title('Acceleration of Flex-Ext Movement');

%% Smoothness calculation FE

% Derivation of acceleration to obtain smoothness
sFERep1 = gradient(aF1(:));
sFERep2 = gradient(aF2(:));
sFERep3 = gradient(aF3(:));

% Smoothness
sF1 = sFERep1 ./ dtF1;
sF2 = sFERep2 ./ dtF2;
sF3 = sFERep3 ./ dtF3;

% Filtering 
fcS = 4; 
[bS, aS] = butter(4, fcS / (Freq/2));
sF1 = filtfilt(bS, aS, sF1);
sF2 = filtfilt(bS, aS, sF2);
sF3 = filtfilt(bS, aS, sF3);

% Concatenate in a single vector
sF = zeros(length(sF1)+length(sF2)+length(sF3),1);
sF(1:length(sF1),1) = sF1;
sF(length(sF1)+1:length(sF1)+length(sF2),1) = sF2;
sF(length(sF1)+length(sF2)+1:length(sF1)+length(sF2)+length(sF3),1) = sF3;

figure(6);
plot(sF);
yline(0, 'm--', 'LineWidth', 1);
title('Smoothness of Flex-Ext Movement');

%% ROTATION
% Delete the OFFSET
mean_aY = mean(PteH.aY(1:9));
PteH_aY_NoOffset = PteH.aY - mean_aY;

% Filter the signal
Freq = 90;
fc = 10;
[b,a] = butter(4,fc/(Freq/2));
PteH_aY_Filtered = filtfilt(b,a,PteH_aY_NoOffset);

figure(7);
subplot(1,2,1);
plot(PteH_aY_Filtered);
yline(0, 'm--', 'LineWidth', 1);
title("PteH aY Filtered");

% Obtain the movement complete (warm up included)
[~, locs_MAX_R] = findpeaks(PteH_aY_Filtered,'MinPeakHeight',40,'MinPeakDistance',100,'SortStr','ascend');
aY_inv = -PteH_aY_Filtered;
[~,locs_MIN_R] = findpeaks(aY_inv,'MinPeakHeight',30,'MinPeakDistance',100,'SortStr','ascend');
peak_locs = sort([locs_MAX_R(:); locs_MIN_R(:)]); % Las pos en el eje X de todos los picos
% Find out if the mov. is correctly done
mean_durations = zeros(1, 6);
for i = 1:6
    mean_durations(i) = round(peak_locs(2*i) - peak_locs(2*i - 1));
end
if length(peak_locs) == 12
    start_index = peak_locs(1) - mean_durations(1);
    end_index = peak_locs(12) + mean_durations(6);
    PteH_aY_Mov = PteH_aY_Filtered(start_index:end_index);
else 
    disp("The ROTATION has been done incorrectly by the participant.");
    return;
end
subplot(1,2,2);
plot(PteH_aY_Mov);
yline(0, 'm--', 'LineWidth', 1);
title("Rotation WITH WARM UP");

clear i start_index end_index locs_MIN_R locs_MAX_R

%% Separate Rotation (without warm up)

first_peak = peak_locs(7); % The first peak of useful mov.
PteH_aY_Rotation = PteH_aY_Filtered(first_peak-mean_durations(4):peak_locs(12)+mean_durations(6));

figure(8);
plot(PteH_aY_Rotation);
yline(0, 'm--', 'LineWidth', 1);
title("ROTATION without warm up");

% The signal without the rest time position
min_segment_length = 150;
contiguous_segments = splitConsecutives(find(abs(diff(PteH_aY_Rotation)) < 0.1));% Buscar zonas sin mucha variacion
rest_position = contiguous_segments(cellfun(@length, contiguous_segments) > min_segment_length);

RRep1 = PteH_aY_Rotation(1:rest_position{1}(1));
RRep2 = PteH_aY_Rotation(rest_position{1}(end):rest_position{end}(1));
RRep3 = PteH_aY_Rotation(rest_position{end}(end):end);

PteH_aY = zeros(length(RRep1)+length(RRep2)+length(RRep3),1);
PteH_aY(1:length(RRep1),1) = RRep1;
PteH_aY(length(RRep1)+1:length(RRep1)+length(RRep2),1) = RRep2;
PteH_aY(length(RRep1)+length(RRep2)+1:length(RRep1)+length(RRep2)+length(RRep3),1) = RRep3;

figure(9);
plot(PteH_aY);
yline(0, 'm--', 'LineWidth', 1);
title("Repetions of Rotation");

clear peak_locs first_peak last_peak mean_durations contiguous_segments rest_position

%% Calculate the CROM R

[~,locs_MAX_R] = findpeaks(PteH_aY,'MinPeakHeight',40,'MinPeakDistance',20);
aY_inv = -PteH_aY;
[~,locs_MIN_R] = findpeaks(aY_inv,'MinPeakHeight',30,'MinPeakDistance',20);
 Rotation = (PteH_aY(locs_MAX_R)-PteH_aY(locs_MIN_R));
 Rotation = transpose(Rotation);
 mR = mean(Rotation);

 clear locs_MIN_R locs_MAX_R

%% Calculate the Velocity for each repetition R

tR1 = (linspace(1, length(RRep1), length(RRep1))/Freq)';
tR2 = (linspace(1, length(RRep2), length(RRep2))/Freq)';
tR3 = (linspace(1, length(RRep3), length(RRep3))/Freq)';

dtR1 = gradient(tR1);
dtR2 = gradient(tR2);
dtR3 = gradient(tR3);

% Derivates --> Instant velocity
dRRep1 = gradient(RRep1(:));
dRRep2 = gradient(RRep2(:));
dRRep3 = gradient(RRep3(:));

% Velocity
vR1 = dRRep1 ./ dtR1;
vR2 = dRRep2 ./ dtR2;
vR3 = dRRep3 ./ dtR3;

% Filtering
fcV = 4;
[bV,aV] = butter(4,fcV/(Freq/2));
vR1 = filtfilt(bV,aV,vR1);
vR2 = filtfilt(bV,aV,vR2);
vR3 = filtfilt(bV,aV,vR3);

% Concatenate in a single vector
vR = zeros(length(vR1)+length(vR2)+length(vR3),1);
vR(1:length(vR1),1) = vR1;
vR(length(vR1)+1:length(vR1)+length(vR2),1) = vR2;
vR(length(vR1)+length(vR2)+1:length(vR1)+length(vR2)+length(vR3),1) = vR3;

figure(10);
plot(vR);
yline(0, 'm--', 'LineWidth', 1);
title("Velocity of Rotation Movement");

% Calculate the max and min MEAN velocity
vR1Max = max(vR1);
vR2Max = max(vR2);
vR3Max = max(vR3);
vRMaxMean = mean([vR1Max;vR2Max;vR3Max]);
vRMax = [vR1Max,vR2Max,vR3Max,vRMaxMean];

vR1Min = min(vR1);
vR2Min = min(vR2);
vR3Min = min(vR3);
vRMinMean = mean([vR1Min,vR2Min,vR3Min]);
vRMin = [vR1Min,vR2Min,vR3Min,vRMinMean];

% Calculate the middle value of the peak velocity
vRabs = abs(vR);
MPDv = int16(length(vR)/9*0.6);
[~,locs_MAX_Rv] = findpeaks(vRabs,'MinPeakHeight',30,'MinPeakDistance',MPDv);
vRm = mean(vRabs(locs_MAX_Rv));

%% Acceleration calculation R

aRRep1 = gradient(vR1(:));
aRRep2 = gradient(vR2(:));
aRRep3 = gradient(vR3(:));

% Acceleration
aR1 = aRRep1 ./ dtR1;
aR2 = aRRep2 ./ dtR2;
aR3 = aRRep3 ./ dtR3;

% Filtering 
fcA = 4; 
[bA, aA] = butter(4, fcA / (Freq/2));
aR1 = filtfilt(bA, aA, aR1);
aR2 = filtfilt(bA, aA, aR2);
aR3 = filtfilt(bA, aA, aR3);

% Concatenate in a single vector
aR = zeros(length(aR1)+length(aR2)+length(aR3),1);
aR(1:length(aR1),1) = aR1;
aR(length(aR1)+1:length(aR1)+length(aR2),1) = aR2;
aR(length(aR1)+length(aR2)+1:length(aR1)+length(aR2)+length(aR3),1) = aR3;

figure(11);
plot(aR);
yline(0, 'm--', 'LineWidth', 1);
title('Acceleration of Rotation Movement');

%% Smoothness calculation R

sRRep1 = gradient(aR1(:));
sRRep2 = gradient(aR2(:));
sRRep3 = gradient(aR3(:));

% Smoothness
sR1 = sRRep1 ./ dtR1;
sR2 = sRRep2 ./ dtR2;
sR3 = sRRep3 ./ dtR3;

% Filtering 
fcS = 4; 
[bS, aS] = butter(4, fcS / (Freq/2));
sR1 = filtfilt(bS, aS, sR1);
sR2 = filtfilt(bS, aS, sR2);
sR3 = filtfilt(bS, aS, sR3);

% Concatenate in a single vector
sR = zeros(length(sR1)+length(sR2)+length(sR3),1);
sR(1:length(sR1),1) = sR1;
sR(length(sR1)+1:length(sR1)+length(sR2),1) = sR2;
sR(length(sR1)+length(sR2)+1:length(sR1)+length(sR2)+length(sR3),1) = sR3;

figure(12);
plot(sR);
yline(0, 'm--', 'LineWidth', 1);
title('Smoothness of Rotation Movement');

%% LATERAL INCLINATION
% Delete the OFFSET
mean_aZ = mean(PteH.aZ(1:9));
PteH_aZ_NoOffset = PteH.aZ - mean_aZ;

% Filter the signal
Freq = 90;
fc = 10;
[b,a] = butter(4,fc/(Freq/2));
PteH_aZ_Filtered = filtfilt(b,a,PteH_aZ_NoOffset);

figure(13);
subplot(1,2,1);
plot(PteH_aZ_Filtered);
yline(0, 'm--', 'LineWidth', 1);
title("PteH aZ Filtered");

% Obtain the movement complete (warm up included)
[~, locs_MAX_LI] = findpeaks(PteH_aZ_Filtered,'MinPeakHeight',30,'MinPeakDistance',100,'SortStr','ascend');
aZ_inv = -PteH_aZ_Filtered;
[~,locs_MIN_LI] = findpeaks(aZ_inv,'MinPeakHeight',25,'MinPeakDistance',100,'SortStr','ascend');
peak_locs = sort([locs_MAX_LI(:); locs_MIN_LI(:)]); % Las pos en el eje X de todos los picos
% Find out if the mov. is correctly done
mean_durations = zeros(1, 6);
for i = 1:6
    mean_durations(i) = round(peak_locs(2*i) - peak_locs(2*i - 1));
end
if length(peak_locs) == 12
    start_index = peak_locs(1) - mean_durations(1);
    end_index = peak_locs(12) + mean_durations(6);
    PteH_aZ_Mov = PteH_aZ_Filtered(start_index:end_index);
else 
    disp("The LATERAL INCLINATION has been done incorrectly by the participant.");
    return;
end
subplot(1,2,2);
plot(PteH_aZ_Mov);
yline(0, 'm--', 'LineWidth', 1);
title("Lateral inclination WITH WARM UP");

 clear i start_index end_index locs_MIN_LI locs_MAX_LI

%% Separate Lateral Inclination (without warm up)

first_peak = peak_locs(7); % The first peak of useful mov.
PteH_aZ_Rotation = PteH_aZ_Filtered(first_peak-mean_durations(4):peak_locs(12)+mean_durations(6));

figure(14);
plot(PteH_aZ_Rotation);
yline(0, 'm--', 'LineWidth', 1);
title("LATERAL INCLINATION without warm up");

% The signal without the rest time position
min_segment_length = 150;
contiguous_segments = splitConsecutives(find(abs(diff(PteH_aZ_Rotation)) < 0.12));% Buscar zonas sin mucha variacion
rest_position = contiguous_segments(cellfun(@length, contiguous_segments) > min_segment_length);

LIRep1 = PteH_aZ_Rotation(1:rest_position{1}(1));
LIRep2 = PteH_aZ_Rotation(rest_position{1}(end):rest_position{end}(1));
LIRep3 = PteH_aZ_Rotation(rest_position{end}(end):end);

PteH_aZ = zeros(length(LIRep1)+length(LIRep2)+length(LIRep3),1);
PteH_aZ(1:length(LIRep1),1) = LIRep1;
PteH_aZ(length(LIRep1)+1:length(LIRep1)+length(LIRep2),1) = LIRep2;
PteH_aZ(length(LIRep1)+length(LIRep2)+1:length(LIRep1)+length(LIRep2)+length(LIRep3),1) = LIRep3;

figure(15);
plot(PteH_aZ);
yline(0, 'm--', 'LineWidth', 1);
title("Repetions of Lateral Inclination");

 clear peak_locs first_peak last_peak mean_durations contiguous_segments rest_position

%% Calculate the CROM R

[~,locs_MAX_LI] = findpeaks(PteH_aZ,'MinPeakHeight',40,'MinPeakDistance',20);
aZ_inv = -PteH_aZ;
[~,locs_MIN_LI] = findpeaks(aZ_inv,'MinPeakHeight',30,'MinPeakDistance',20);
 LateralInclination = (PteH_aZ(locs_MAX_LI)-PteH_aZ(locs_MIN_LI));
 LateralInclination = transpose(LateralInclination);
 mLI = mean(LateralInclination);

 clear locs_MIN_LI locs_MAX_LI

%% Calculate the Velocity for each repetition LI

tLI1 = (linspace(1, length(LIRep1), length(LIRep1))/Freq)';
tLI2 = (linspace(1, length(LIRep2), length(LIRep2))/Freq)';
tLI3 = (linspace(1, length(LIRep3), length(LIRep3))/Freq)';

dtLI1 = gradient(tLI1);
dtLI2 = gradient(tLI2);
dtLI3 = gradient(tLI3);

% Derivates --> Instant velocity
dLIRep1 = gradient(LIRep1(:));
dLIRep2 = gradient(LIRep2(:));
dLIRep3 = gradient(LIRep3(:));

% Velocity
vLI1 = dLIRep1 ./ dtLI1;
vLI2 = dLIRep2 ./ dtLI2;
vLI3 = dLIRep3 ./ dtLI3;

% Filtering
fcV = 4;
[bV,aV] = butter(4,fcV/(Freq/2));
vLI1 = filtfilt(bV,aV,vLI1);
vLI2 = filtfilt(bV,aV,vLI2);
vLI3 = filtfilt(bV,aV,vLI3);

% Concatenate in a single vector
vLI = zeros(length(vLI1)+length(vLI2)+length(vLI3),1);
vLI(1:length(vLI1),1) = vLI1;
vLI(length(vLI1)+1:length(vLI1)+length(vLI2),1) = vLI2;
vLI(length(vLI1)+length(vLI2)+1:length(vLI1)+length(vLI2)+length(vLI3),1) = vLI3;

figure(16);
plot(vLI);
yline(0, 'm--', 'LineWidth', 1);
title("Velocity of Lateral Inclination Movement");

% Calculate the max and min MEAN velocity
vLI1Max = max(vLI1);
vLI2Max = max(vLI2);
vLI3Max = max(vLI3);
vLIMaxMean = mean([vLI1Max;vLI2Max;vLI3Max]);
vLIMax = [vLI1Max,vLI2Max,vLI3Max,vLIMaxMean];

vLI1Min = min(vLI1);
vLI2Min = min(vLI2);
vLI3Min = min(vLI3);
vLIMinMean = mean([vLI1Min,vLI2Min,vLI3Min]);
vLIMin = [vLI1Min,vLI2Min,vLI3Min,vLIMinMean];

% Calculate the middle value of the peak velocity
vLIabs = abs(vLI);
MPDv = int16(length(vLI)/9*0.6);
[~,locs_MAX_LIv] = findpeaks(vLIabs,'MinPeakHeight',30,'MinPeakDistance',MPDv);
vLIm = mean(vLIabs(locs_MAX_LIv));

%% Acceleration calculation LI

aLIRep1 = gradient(vLI1(:));
aLIRep2 = gradient(vLI2(:));
aLIRep3 = gradient(vLI3(:));

% Acceleration
aLI1 = aLIRep1 ./ dtLI1;
aLI2 = aLIRep2 ./ dtLI2;
aLI3 = aLIRep3 ./ dtLI3;

% Filtering 
fcA = 4; 
[bA, aA] = butter(4, fcA / (Freq/2));
aLI1 = filtfilt(bA, aA, aLI1);
aLI2 = filtfilt(bA, aA, aLI2);
aLI3 = filtfilt(bA, aA, aLI3);

% Concatenate in a single vector
aLI = zeros(length(aLI1)+length(aLI2)+length(aLI3),1);
aLI(1:length(aLI1),1) = aLI1;
aLI(length(aLI1)+1:length(aLI1)+length(aLI2),1) = aLI2;
aLI(length(aLI1)+length(aLI2)+1:length(aLI1)+length(aLI2)+length(aLI3),1) = aLI3;

figure(17);
plot(aLI);
yline(0, 'm--', 'LineWidth', 1);
title('Acceleration of Lateral Inclination Movement');

%% Smoothness calculation LI

sLIRep1 = gradient(aLI1(:));
sLIRep2 = gradient(aLI2(:));
sLIRep3 = gradient(aLI3(:));

% Smoothness
sLI1 = sLIRep1 ./ dtLI1;
sLI2 = sLIRep2 ./ dtLI2;
sLI3 = sLIRep3 ./ dtLI3;

% Filtering 
fcS = 4; 
[bS, aS] = butter(4, fcS / (Freq/2));
sLI1 = filtfilt(bS, aS, sLI1);
sLI2 = filtfilt(bS, aS, sLI2);
sLI3 = filtfilt(bS, aS, sLI3);

% Concatenate in a single vector
sLI = zeros(length(sLI1)+length(sLI2)+length(sLI3),1);
sLI(1:length(sLI1),1) = sLI1;
sLI(length(sLI1)+1:length(sLI1)+length(sLI2),1) = sLI2;
sLI(length(sLI1)+length(sLI2)+1:length(sLI1)+length(sLI2)+length(sLI3),1) = sLI3;

figure(18);
plot(sLI);
yline(0, 'm--', 'LineWidth', 1);
title('Smoothness of Lateral Inclination Movement');

%% Comparison of the 3 movements 

figure(19);
min_length = min([length(PteH_aX), length(PteH_aY), length(PteH_aZ)]);
PteH_aX = interp1(1:length(PteH_aX), PteH_aX, linspace(1, length(PteH_aX), min_length));
PteH_aY = interp1(1:length(PteH_aY), PteH_aY, linspace(1, length(PteH_aY), min_length));
PteH_aZ = interp1(1:length(PteH_aZ), PteH_aZ, linspace(1, length(PteH_aZ), min_length));

hold on;  
plot(PteH_aX, 'LineWidth', 1, 'Color', 'blue', 'DisplayName', 'Flexion-Extension (X)');
plot(PteH_aY, 'LineWidth', 1, 'Color', 'green', 'DisplayName', 'Rotation (Y)');
plot(PteH_aZ, 'LineWidth', 1, 'Color', 'red', 'DisplayName', 'Lateral Inclination (Z)');
yline(0, 'm--', 'LineWidth', 1, 'HandleVisibility', 'off');
title("Comparison of the 3 different movs.");
xlabel('Seconds');
ylabel('Degrees');
legend('show');
grid on;

%% FUNCTIONS

function segments = splitConsecutives(indices)
    segments = {};
    current_segment = [];
    for i = 1:length(indices)-1
        if indices(i+1) - indices(i) == 1
            current_segment = [current_segment, indices(i)];
        else
            if ~isempty(current_segment)
                current_segment = [current_segment, indices(i)]; % Añadir el último índice al segmento
                segments{end+1} = current_segment;
                current_segment = [];
            end
        end
    end
    if ~isempty(current_segment)
        current_segment = [current_segment, indices(end)]; % Añadir el último índice al segmento
        segments{end+1} = current_segment;
    end
end