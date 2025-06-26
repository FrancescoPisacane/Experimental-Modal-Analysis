clc
close all
clear all 


%% Loading of the files

% selecting main folder where are all the data subfolder sets
mainFolder = 'C:\Users\franc\Desktop\Data Analysis\Progetto\test modal analysis';
% mainFolder = 'C:\Users\Utente\Desktop\LAUREA_MAGISTRALE\SEMESTER_2\DATA ANALYSIS FOR MECHANICAL SYSTEM IDENTIFICATION\PROJECT\experimental_data';
% mainFolder = 'C:\Users\Carlo Tundo\Desktop\PoliMi\Master 1st year\DATA ANALYSIS\assignement\modal analysis\test modal analysis';
folderList = genpath(mainFolder);           % Get all paths
folderListCell = strsplit(folderList, pathsep);  % Split into individual folders
folderListCell = folderListCell(~cellfun('isempty', folderListCell));  % Remove empty entries

folderListCell(1) = []; % deleting main folder path 

% definition of the dimension
ntest = 10;
nset = 5;
test = cell(nset, 1);

for ii = 1:nset % for each set 
    % preallocation 
    test{ii}.data = cell(ntest, 1);
    % selceting folder 
    currentFolder = folderListCell{ii};
    files = dir(fullfile(currentFolder, '*.mat'));
    
    % acquisition data
    for jj = 1:length(files)
        filePath = fullfile(currentFolder, files(jj).name);
        test{ii}.data{jj} = load(filePath);
    end
end

n_acc = size(test{1}.data{1}.Dati, 2) - 1;

%% Sensitivity of the dataset and creation of time vector
hammer_sens_steel = 2.488; % [mV/N]
hammer_sens_none = 2.361; % [mV/N]
sens_hammer = hammer_sens_steel; % default choice

sens_acc = 10.2; % [mV/(m/s^2)]
fmax_acc = 3000; % [Hz]

% definition of EU datas
for ii = 1:nset
    for jj = 1:ntest
        test{ii}.data{jj}.Dati(:,1) = test{ii}.data{jj}.Dati(:,1)./sens_hammer;
        test{ii}.data{jj}.Dati(:,2) = test{ii}.data{jj}.Dati(:,2)./sens_acc;
        test{ii}.data{jj}.Dati(:,3) = test{ii}.data{jj}.Dati(:,3)./sens_acc;
        test{ii}.data{jj}.Dati(:,4) = test{ii}.data{jj}.Dati(:,4)./sens_acc;
        test{ii}.data{jj}.Dati(:,5) = test{ii}.data{jj}.Dati(:,5)./sens_acc;
        test{ii}.data{jj}.Dati(:,6) = test{ii}.data{jj}.Dati(:,6)./sens_acc;
    end
end

% f sampling
safety_factor_aliasing = 2.56;
% fsamp = 810 * safety_factor_aliasing; % seen from mode 6 abaqus 
fsamp = 2500;

% creation of time vector
dt = 1/fsamp;
nsamples = length(test{1}.data{1}.Dati(:,1));
time = (0:dt:(nsamples*dt-dt))';

%% windowing of time domain signals - Tukey and Exponential
treshold = 0.02; % N
duration_test = 25; % seconds
pretrigger_duration = 0.1; % seconds

for ii = 1:nset
    for jj = 1:ntest
        idx = find(test{ii}.data{jj}.Dati(:,1) > treshold, 1); % index at which input signal is higher than treshold
        idx_list = (idx-pretrigger_duration*fsamp):1:(idx - pretrigger_duration*fsamp + duration_test*fsamp);
        test{ii}.data{jj}.Dati = [test{ii}.data{jj}.Dati(idx_list,1), test{ii}.data{jj}.Dati(idx_list,2), test{ii}.data{jj}.Dati(idx_list,3), test{ii}.data{jj}.Dati(idx_list,4), test{ii}.data{jj}.Dati(idx_list,5), test{ii}.data{jj}.Dati(idx_list,6)];
        % data{ii}.Dati = [data{ii}.Dati(idx_list,1), data{ii}.Dati(idx_list,2), data{ii}.Dati(idx_list,3)];
    end
end

N = length(test{1}.data{1}.Dati(: ,1));
% windowing with tukey window
w_tukey = [tukeywin(2*pretrigger_duration*fsamp); zeros((N-2*pretrigger_duration*fsamp), 1)];

% exponential window 
T = N/fsamp;          % seconds duration
t = linspace(0,(T-1), N);    % time axis
t0 = pretrigger_duration;         % initial time sample
P = 0.01;                       % percentage of final value
tau = (T - t0) / log(1/P);      % time constant

% defining exponential window 
w_exp = ones(size(t));
w_exp(t >= t0) = exp(-(t(t >= t0) - t0)/tau);
w_exp = w_exp';

% windowing 
for ii = 1:nset
    for jj = 1:ntest
        for kk = 1:n_acc
            if kk == 1
                test{ii}.data{jj}.Dati(:,kk) = w_tukey.*test{ii}.data{jj}.Dati(:,kk); % se non esegue prova a mettere stoppino su questa riga ed eseguire manualmente. Dovrebbe sbloccarsi
            else
                % data{ii}.Dati(:,jj) = w_exp.*data{ii}.Dati(:,jj);
            end
        end
    end
end

%% plotting time - domain
% optional
% set = 5; % default
% 
% plotting_title = [{'Input'}, {'Acc1'}, {'Acc2'}, {'Acc3'}, {'Acc4'}, {'Acc5'}];
% for jj = 1:ntest
%     figure('Name', ['Test ' num2str(jj)], 'NumberTitle', 'off');
%     plotting_cell = {test{set}.data{jj}.Dati(:, 1), test{set}.data{jj}.Dati(:, 2), test{set}.data{jj}.Dati(:, 3), test{set}.data{jj}.Dati(:, 4), test{set}.data{jj}.Dati(:, 5), test{set}.data{jj}.Dati(:, 6)};
%     for kk = 1:(n_acc + 1)
%         subplot(n_acc+1,1,kk)
%         plot(t, plotting_cell{kk}, 'LineWidth', 1.2)
%         title(plotting_title{kk})
%         xlabel ('f [Hz]')
%         ylabel('EU')
%         grid on 
%     end
%     sgtitle(['Time domain plots: set ', num2str(set)]);
% end

%% Frequency domain calculation
% preallocation of set struct 
FFT = cell(nset, 1);
Pspectra = cell(nset, 1);
Cspectra = cell(nset, 1);
CC_Cspectra = cell(nset, 1);
Phase = cell(nset, 1);
H1 = cell(nset, 1);
H2 = cell(nset, 1);
Coherence = cell(nset, 1);
neg = 0; %to cancel out negative frequencies

for ii = 1:nset             % for each set 
    % preallocation of ffts window 
    FFT{ii}.input_ffts = cell(ntest,1);
    FFT{ii}.acc1_ffts = cell(ntest, 1);
    FFT{ii}.acc2_ffts = cell(ntest, 1);
    FFT{ii}.acc3_ffts = cell(ntest, 1);
    FFT{ii}.acc4_ffts = cell(ntest, 1);
    FFT{ii}.acc5_ffts = cell(ntest, 1);
    
    % population of the tests
    for jj = 1:ntest
        FFT{ii}.input_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 1));
        FFT{ii}.acc1_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 2));
        FFT{ii}.acc2_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 3));
        FFT{ii}.acc3_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 4));
        FFT{ii}.acc4_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 5));
        FFT{ii}.acc5_ffts{jj} = fft(test{ii}.data{jj}.Dati(:, 6));
    end
    
    % Power spectrum calculation
    Pspectra{ii}.input_Ps = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.input_ffts, neg);
    Pspectra{ii}.acc_Ps1 = Ps_Cs(FFT{ii}.acc1_ffts, FFT{ii}.acc1_ffts, neg);
    Pspectra{ii}.acc_Ps2 = Ps_Cs(FFT{ii}.acc2_ffts, FFT{ii}.acc2_ffts, neg);
    Pspectra{ii}.acc_Ps3 = Ps_Cs(FFT{ii}.acc3_ffts, FFT{ii}.acc3_ffts, neg);
    Pspectra{ii}.acc_Ps4 = Ps_Cs(FFT{ii}.acc4_ffts, FFT{ii}.acc4_ffts, neg);
    Pspectra{ii}.acc_Ps5 = Ps_Cs(FFT{ii}.acc5_ffts, FFT{ii}.acc5_ffts, neg);
    
    % Cross spectrum calculation
    Cspectra{ii}.acc_Cs1 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc1_ffts, neg);
    Cspectra{ii}.acc_Cs2 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc2_ffts, neg);
    Cspectra{ii}.acc_Cs3 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc3_ffts, neg);
    Cspectra{ii}.acc_Cs4 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc4_ffts, neg);
    Cspectra{ii}.acc_Cs5 = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc5_ffts, neg);
    
    % complex conjugate cross spectrum 
    CC_Cspectra{ii}.acc_Cs1 = Ps_Cs(FFT{ii}.acc1_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs2 = Ps_Cs(FFT{ii}.acc2_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs3 = Ps_Cs(FFT{ii}.acc3_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs4 = Ps_Cs(FFT{ii}.acc4_ffts, FFT{ii}.input_ffts, neg);
    CC_Cspectra{ii}.acc_Cs5 = Ps_Cs(FFT{ii}.acc5_ffts, FFT{ii}.input_ffts, neg);
    
    % Phase 
    Phase{ii}.phase_acc1 = angle(Cspectra{ii}.acc_Cs1);
    Phase{ii}.phase_acc2 = angle(Cspectra{ii}.acc_Cs2);
    Phase{ii}.phase_acc3 = angle(Cspectra{ii}.acc_Cs3);
    Phase{ii}.phase_acc4 = angle(Cspectra{ii}.acc_Cs4);
    Phase{ii}.phase_acc5 = angle(Cspectra{ii}.acc_Cs5);
    
    
    % ------> Extimation of FRFs ------<
    
    % H1
    H1{ii}.abs_out1 = abs(Cspectra{ii}.acc_Cs1)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out2 = abs(Cspectra{ii}.acc_Cs2)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out3 = abs(Cspectra{ii}.acc_Cs3)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out4 = abs(Cspectra{ii}.acc_Cs4)./Pspectra{ii}.input_Ps;
    H1{ii}.abs_out5 = abs(Cspectra{ii}.acc_Cs5)./Pspectra{ii}.input_Ps;

    % H2
    H2{ii}.abs_out1 = abs(Pspectra{ii}.acc_Ps1)./abs(CC_Cspectra{ii}.acc_Cs1);
    H2{ii}.abs_out2 = abs(Pspectra{ii}.acc_Ps2)./abs(CC_Cspectra{ii}.acc_Cs2);
    H2{ii}.abs_out3 = abs(Pspectra{ii}.acc_Ps3)./abs(CC_Cspectra{ii}.acc_Cs3);
    H2{ii}.abs_out4 = abs(Pspectra{ii}.acc_Ps4)./abs(CC_Cspectra{ii}.acc_Cs4);
    H2{ii}.abs_out5 = abs(Pspectra{ii}.acc_Ps5)./abs(CC_Cspectra{ii}.acc_Cs5);
    
    % coherence 
    Coherence{ii}.gamma_acc1 = (Cspectra{ii}.acc_Cs1.*CC_Cspectra{ii}.acc_Cs1)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps1);
    Coherence{ii}.gamma_acc2 = (Cspectra{ii}.acc_Cs2.*CC_Cspectra{ii}.acc_Cs2)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps2);
    Coherence{ii}.gamma_acc3 = (Cspectra{ii}.acc_Cs3.*CC_Cspectra{ii}.acc_Cs3)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps3);
    Coherence{ii}.gamma_acc4 = (Cspectra{ii}.acc_Cs4.*CC_Cspectra{ii}.acc_Cs4)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps4);
    Coherence{ii}.gamma_acc5 = (Cspectra{ii}.acc_Cs5.*CC_Cspectra{ii}.acc_Cs5)./(Pspectra{ii}.input_Ps.*Pspectra{ii}.acc_Ps5);
end

%% plotting H1, H2, Coherence and Phase for a chosen couple input-sensor
% frequency vector 
f = linspace(0, 900, length(Pspectra{1}.input_Ps)*(900/(fsamp/2))); % stop after six peaks
f2 = linspace(0, fsamp/2, length(Pspectra{1}.input_Ps)); % overall
% user input 
idx_set = input('Select a input (1-5): ');
idx_acc = input('Select a sensor (1-5): ');

% preallocating plotting vector
plotting_vector_h1 = {H1{idx_set}.abs_out1, H1{idx_set}.abs_out2, H1{idx_set}.abs_out3, H1{idx_set}.abs_out4, H1{idx_set}.abs_out5};
plotting_vector_h2 = {H2{idx_set}.abs_out1, H2{idx_set}.abs_out2, H2{idx_set}.abs_out3, H2{idx_set}.abs_out4, H2{idx_set}.abs_out5};
plotting_vector_coh = {Coherence{idx_set}.gamma_acc1, Coherence{idx_set}.gamma_acc2, Coherence{idx_set}.gamma_acc3, Coherence{idx_set}.gamma_acc4, Coherence{idx_set}.gamma_acc5};
plotting_vector_phase = {Phase{idx_set}.phase_acc1, Phase{idx_set}.phase_acc2, Phase{idx_set}.phase_acc3, Phase{idx_set}.phase_acc4, Phase{idx_set}.phase_acc5};

figure
subplot(2,1,1)
semilogy(f, plotting_vector_h1{idx_acc}(1:length(f)), 'LineWidth', 1.2,'Color', 'b')
hold on 
semilogy(f, plotting_vector_h2{idx_acc}(1:length(f)), 'LineWidth', 1.2, 'Color', 'r')
xlabel('frequency [Hz]')
ylabel('abs')
title(['H1, H2 and Coherence for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])
yyaxis right
ax = gca;
ax.YAxis(2).Color = 'k';
plot(f,plotting_vector_coh{idx_acc}(1:length(f)), 'LineWidth', 1.2, 'Color', 'k')
ylim([0, 3])
yline(1, '--k', 'LineWidth', 1.2)
grid on 
legend('H1', 'H2', 'Coherence', 'Location','northwest')

ylabel('Coherence')
hold off

subplot(2,1,2)
plot(f,plotting_vector_phase{idx_acc}(1:length(f)), 'LineWidth', 1.2, 'Color', 'b')
grid on
xlabel('frequency [Hz]')
ylabel('angle [rad]')
title(['Phase for Input location ', num2str(idx_set), ' and Output accelerometer ', num2str(idx_acc)])


%% plotting power and cross spectrum
% % power spectra plot -- OPTIONAL 
% plotting_cell = {input_Ps, acc_Ps1, acc_Ps2, acc_Ps3, acc_Ps4, acc_Ps5};
% % plotting_cell = {input_Ps, acc_Ps1, acc_Ps2};
% plotting_title = [{'PS-input'}, {'PS-out1'}, {'PS-out2'}, {'PS-out3'}, {'PS-out4'}, {'PS-out5'}];
% 
% figure('Name', 'Power Spectra', 'NumberTitle', 'off');
% for jj = 1:(n_acc+1)
%     subplot(n_acc+1,1,jj)
%     semilogy(f,abs(plotting_cell{jj}), 'LineWidth', 1.2)
%     grid on
%     hold on 
%     xlabel('f')
%     ylabel('EU')
%     title(plotting_title(jj))
% end
% 
% % cross spectra plot
% % plotting_cell = {input_Ps, acc_Ps1, acc_Ps2, acc_Ps3, acc_Ps4, acc_Ps5};
% plotting_cell = {input_Ps, acc_Cs1, acc_Cs2, acc_Cs3, acc_Cs4, acc_Cs5};
% plotting_title = [{'PS-input'}, {'CS-out1'}, {'CS-out2'}, {'CS-out3'}, {'CS-out4'}, {'CS-out5'}];
% 
% figure('Name', 'Cross Spectra', 'NumberTitle', 'off');
% for jj = 1:(n_acc+1)
%     subplot(n_acc+1,1,jj)
%     semilogy(f,abs(plotting_cell{jj}), 'LineWidth', 1.2)
%     grid on
%     hold on 
%     xlabel('f')
%     ylabel('EU')
%     title(plotting_title(jj))
% end


%% MDOF 

% We choose H1 because we believe we have more noise on the output
n_inp = n_acc;
Pspectra_n = cell(nset, 1);
Cspectra_n = cell(nset, 1);
neg = 1; %now we want to keep negative frequencies for doing the inverse F.T.
H1_C = cell(n_inp, 1);

for ii = 1:nset 
    % Power spectrum calculation
    Pspectra_n{ii}.input_Ps_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.input_ffts, neg);
    Pspectra_n{ii}.acc_Ps1_n = Ps_Cs(FFT{ii}.acc1_ffts, FFT{ii}.acc1_ffts, neg);
    Pspectra_n{ii}.acc_Ps2_n = Ps_Cs(FFT{ii}.acc2_ffts, FFT{ii}.acc2_ffts, neg);
    Pspectra_n{ii}.acc_Ps3_n = Ps_Cs(FFT{ii}.acc3_ffts, FFT{ii}.acc3_ffts, neg);
    Pspectra_n{ii}.acc_Ps4_n = Ps_Cs(FFT{ii}.acc4_ffts, FFT{ii}.acc4_ffts, neg);
    Pspectra_n{ii}.acc_Ps5_n = Ps_Cs(FFT{ii}.acc5_ffts, FFT{ii}.acc5_ffts, neg);

    % Cross spectrum calculation
    Cspectra_n{ii}.acc_Cs1_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc1_ffts, neg);
    Cspectra_n{ii}.acc_Cs2_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc2_ffts, neg);
    Cspectra_n{ii}.acc_Cs3_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc3_ffts, neg);
    Cspectra_n{ii}.acc_Cs4_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc4_ffts, neg);
    Cspectra_n{ii}.acc_Cs5_n = Ps_Cs(FFT{ii}.input_ffts, FFT{ii}.acc5_ffts, neg);

    % H1-COMPLEX
    H1_C{ii, 1}.out1 = Cspectra_n{ii}.acc_Cs1_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out2 = Cspectra_n{ii}.acc_Cs2_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out3 = Cspectra_n{ii}.acc_Cs3_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out4 = Cspectra_n{ii}.acc_Cs4_n./Pspectra_n{ii}.input_Ps_n;
    H1_C{ii, 1}.out5 = Cspectra_n{ii}.acc_Cs5_n./Pspectra_n{ii}.input_Ps_n;

end

%% Matrix h for each time instant

h = zeros(n_inp, n_acc, N);

for ii = 1:n_inp
    h(1, ii, :) = ifft(H1_C{ii, 1}.out1);
    h(2, ii, :) = ifft(H1_C{ii, 1}.out2);
    h(3, ii, :) = ifft(H1_C{ii, 1}.out3);
    h(4, ii, :) = ifft(H1_C{ii, 1}.out4);
    h(5, ii, :) = ifft(H1_C{ii, 1}.out5);
end

% if you want to see the unit impulse response 
% figure
% B = squeeze(h(1,5,:));
% plot(t, B)

%% Now we cut the h at t = 20sec

t_cut = 0:dt:((20*length(t)/(25))*dt-dt);
h_cut = h(:, :, 1:(20*length(t)/(25)));
figure
B_cut = squeeze(h_cut(2,5,:));
plot(t_cut, B_cut)

%% Now we define Hmn

% NB!!!!!!!!!!!
%THIS PART HAS NOT BEEN CANCELLED BECAUSE I'M NOT SURE, SO I KEEP IT BUT
%THE PART THAT WE HAVE TO RUN IS Whole stabilization diagram


%tot = length(t_cut);
% m = 1000;
% n = tot - m;
% 
% Hmn = zeros(m*n_acc, n*n_inp);
% 
% for ii = 1:m
%     for jj = 1:n
%         Hmn((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2);
%     end
% end
% 
% 
% 
% 
% Hmn_next= zeros(m*n_acc, n*n_inp);
% 
% for ii = 1:m
%     for jj = 1:n
%         Hmn_next((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2 + 1);
%     end
% end

%% Now we do the diagonalization


% NB!!!!!!!!!!!
%THIS PART HAS NOT BEEN CANCELLED BECAUSE I'M NOT SURE, SO I KEEP IT BUT
%THE PART THAT WE HAVE TO RUN IS Whole stabilization diagram


% Hmn_pseudo = pinv(Hmn);
% Hmm_final = Hmn_next*Hmn_pseudo;
% 
% [phi, eigenvalues] = eig(Hmm_final);
% 
% poles = (log(diag(eigenvalues)))./dt;
% freq = abs(poles)./(2*pi);
% 
% 
% x = 1;                 % larghezza desiderata per ogni rettangolo
% 
% % Calcolo i bin edges con passo x
% edges = min(freq):x:810;  % +x per includere l'ultimo dato
% 
% figure
% % Crea l'istogramma
% histogram(freq, 'BinEdges', edges);


%% Frequency stabilization:

tot = length(t_cut);
m_values = [12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132];
freq_stab = [];
poles_tot = [];
phi = [];

for idx = 1:length(m_values)
    m = m_values(idx);
    n = tot - m;

    Hmn = zeros(m*n_acc, n*n_inp);
    Hmn_next = zeros(m*n_acc, n*n_inp);

    for ii = 1:m % filling the matrix inserting a 5x5 each time
        for jj = 1:n
            Hmn((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2);
            Hmn_next((ii*5-4):ii*5, (jj*5-4):jj*5) = h_cut(:, :, 1 + ii + jj-2 + 1);
        end
    end

    H_pinv = pinv(Hmn);
    Hmm = Hmn_next * H_pinv;
    [eigenvectors, eigenvalues] = eig(Hmm);
    poles = (log(diag(eigenvalues)))./dt;
    freq = abs(poles)./(2*pi);
    poles_tot = [poles_tot; poles];
    freq_stab = [freq_stab; freq];
    phi_i = eigenvectors(1:5, :);
    phi = [phi, phi_i];
end

freq_res = 0.6; %width of the frequency range for the stabilization diagram 0.1%
edges = min(freq):freq_res:820;

figure
histogram(freq_stab, 'BinEdges', edges)



%% Damping stabilization:
freq_stable_poles = cell(6, 1); %we put here all the poles that are stable in frequency
freq_stable_modes = cell(6, 1);
damping_cell = cell(6, 1);
freq_range = [247.2, 327.6, 391.8, 596.4, 620.4, 806.4];


for ii = 1:length(freq_range)
    kk = 0;
    for jj = 1:length(poles_tot)
        
        current_mode = phi(:, jj);
        current_pole = poles_tot(jj);
        current_freq = abs(current_pole)./(2*pi);

        if(current_freq >= freq_range(ii) && current_freq <= (freq_range(ii)+freq_res))
            kk = kk + 1;
            freq_stable_poles{ii, 1}(kk) = current_pole;
            freq_stable_modes{ii, 1}(:, kk) = current_mode;
            %we compile the damping cell 
            current_damping = (-1)*real(current_pole)./(current_freq*2*pi);
            damping_cell{ii, 1}(kk) = current_damping;
        end
    end
end


mostFrequentInterval_damping = zeros(6, 2);
%now we plot the stabilization diagrams for each mode for the damping
for ii = 1:6
    max_damping = max(damping_cell{ii, 1}(:));
    damping_res = (0.5/100)*max_damping;
    edges = min(damping_cell{ii, 1}(:)) : damping_res : max_damping;

    figure
    histogram(damping_cell{ii, 1}(:), 'BinEdges', edges)

    % Calcolo dell'istogramma (senza plottarlo)
    [counts, ~] = histcounts(damping_cell{ii, 1}(:), edges);

    % Trovo il massimo conteggio
    [maxCount, idxMax] = max(counts);

    % Estrazione dell'intervallo associato al valore più frequente
    mostFrequentInterval_damping(ii,:) = [edges(idxMax), edges(idxMax+1)];
    
end

freq_damping_stable_poles = cell(6, 1);
freq_damping_stable_modes = cell(6, 1);

for ii = 1:6
    kk = 0;
    for jj = 1:length(freq_stable_poles{ii, 1}(:))

        current_mode = freq_stable_modes{ii, 1}(:, jj);
        current_freq = abs(freq_stable_poles{ii, 1}(jj))./(2*pi);
        current_damping = (-1)*real(freq_stable_poles{ii, 1}(jj))./(current_freq*2*pi);

        if(current_damping >= mostFrequentInterval_damping(ii,1) && current_damping <= mostFrequentInterval_damping(ii,2))
            kk = kk + 1;

            freq_damping_stable_modes{ii, 1}(:, kk) = freq_stable_modes{ii, 1}(:, jj);
            freq_damping_stable_poles{ii, 1}(kk) = freq_stable_poles{ii, 1}(jj);
        end
    end
end


%% Modes

%siamo in caso di damping proporzionale perchè abbiamo un oggetto metallico
%e quindi ci aspettiamo di avere una parte immaginaria bassa nei modi,
%infatti è così!! Quindi quello che facciamo è scartare la parte
%immaginaria (un ordine di grandezza in meno) e tenere solo quella reale.
% Se avessimo avuto modi complessi con parte immaginaria non trascurabile 
%avremmo trovato una massa modale complessa e il rapporto tornava sensato
%di nuovo 

freq_damping_stable_modes_real = cell(6, 1);
for ii = 1:6
    freq_damping_stable_modes_real{ii, 1}(:, :) = real(freq_damping_stable_modes{ii, 1}(:, :));
end


MAC = [];
for ii = 1:6
    %we choose as a reference one mode and we compute the mac with respect
    %to this reference mode for all the other modes of the ii-th component
    reference_mode = freq_damping_stable_modes_real{ii, 1}(:, 1);
    k = 0;
    
    numerator = [];
    denominator = [];
    for jj = 1:length(freq_damping_stable_modes_real{ii, 1}(1, :))
        current_mode = freq_damping_stable_modes_real{ii, 1}(:, jj);

        numerator = abs(reference_mode' * current_mode)^2;
        denominator = (reference_mode' * reference_mode) * (current_mode' * current_mode);
        MAC(ii, jj) = numerator / denominator;


        if(MAC > 0.9)
            k = k + 1;
        end
    end
end
        

%SINCE WE DID A GREAT JOB BEFORE THE freq_damping_stable_poles ARE EXACTLY
%EQUAL TO THE POLES ASSOCIATED WITH THE POLES THAT HAVE STABLE MODES (MAC
%BIG ENOUGH FOR ALL) IL MAC NON SI ABBASSA MAI, RIMANE SEMPRE ALTO PERCHE
%ABBIAMO FILTRATO MOLTO NEI PASSAGGI PRECEDENTI DI FREQUENZA E DAMPING.

final_stable_poles = freq_damping_stable_poles;
final_stable_modes = freq_damping_stable_modes_real;


%% Averages:

f0_average = [];
csi_average = [];
modes_average = [];
poles_average = [];


for ii = 1:6
    
    %we have complex and conjugated poles, in order to do the average we
    %need to split in 2 sets of poles otherwise we get a real number
    %because the imaginary parts cancel out each other
    current_pole1 = final_stable_poles{ii, 1}(1:2:end);
    current_pole2 = final_stable_poles{ii, 1}(2:2:end);
    poles_average(2*ii - 1,1) = mean(current_pole1);
    poles_average(2*ii,1) = mean(current_pole2);


    current_f0 = abs(final_stable_poles{ii, 1}(:))./(2*pi);
    f0_average(2*ii - 1, 1) = mean(current_f0);
    f0_average(2*ii, 1) = mean(current_f0);

    current_csi = (-1)*real(final_stable_poles{ii, 1}(:))./(current_f0*2*pi);
    csi_average(2*ii-1, 1) = mean(current_csi);
    csi_average(2*ii, 1) = mean(current_csi);

    current_mode = final_stable_modes{ii, 1}(:, :);
    modes_average(2*ii - 1, :) = mean(current_mode, 2);
    modes_average(2*ii, :) = mean(current_mode, 2);

end



%% Modal mass:

N_modes = 6;
f_new = f(4:end);
A = zeros(length(f_new),(2*N_modes+2));
X = zeros((2*N_modes + 2), 1);
B = zeros(N_modes,1);

B = H1_C{1, 1}.out5(1:length(f_new)); %this is Hpq where p = sensor and q = input

for ii = 1:(2*N_modes)

    %in this case we use 5 because the sensor is the 5-th
    column_numerator = ones(length(f_new), 1).*modes_average(ii, 5).^2;
    column_denominator = 1i.*(2*pi.*f_new') - (ones(length(f_new), 1).*poles_average(ii));

    A(:, ii) = column_numerator./column_denominator;

end

A(:, end-1) = ones(length(f_new), 1);

A(:, end) = 1./((f_new.*2*pi).^2);


X = pinv(A) * B;  % pseudoinversa più robusta

%% Finding m_modal:

Q_vector = X(1:2:12, 1);
modal_masses = [];

for ii = 1 : N_modes

    modal_masses(ii) = 1/(Q_vector(ii)*2*1i*f0_average(ii*2-1)*2*pi*sqrt(round(1-(csi_average(ii*2-1))^2)));

end


%% FRF estimation:

H_pred_modal = 0;
for ii = 1:6
    H_pred_i = (modes_average(2*ii, 5)*modes_average(2*ii, 5))./(modal_masses(ii)*(-(2*pi.*f).^2 + 2*1i.*(2*pi.*f).*(2*pi*f0_average(2*ii)).*csi_average(2*ii) + (2*pi*f0_average(2*ii))^2));
    H_pred_modal = H_pred_modal + H_pred_i;
end

H_pred = A*X;

figure
semilogy(f_new, abs(H_pred_modal(1:length(f_new))), 'LineWidth', 1.2, 'Color', 'r');
hold on
semilogy(f_new, abs(H1_C{1, 1}.out5(1:length(f_new))), 'LineWidth', 1.2, 'Color', 'b');
legend('Mdof', 'Experimental', 'Location','northwest')
semilogy(f_new, abs(H_pred(1:end))', 'LineWidth', 1.2, 'Color', 'k');

 %% circle calculation
% npeaks = 6;
% ii = 3;
% cell = {H1_C{ii}.out1, H1_C{ii}.out2, H1_C{ii}.out3, H1_C{ii}.out4, H1_C{ii}.out5};
% 
%     frf = cell{3};
%     real_part_frf = real(frf);
%     imag_part_frf = imag(frf);
% 
%     [pks, idx_pks] = findpeaks(abs(frf), "MinPeakHeight", 40);
% 
%     pks = pks(1:7);
%     idx_pks = idx_pks(1:7);
% 
%     pks(2) = [];
%     idx_pks(2) = [];
% 
%     freq0 = f2(idx_pks);
%     min_f = freq0 - 20;
%     max_f = freq0 + 20;
% 
%     for zz = 1:npeaks
%         [~, idx_min] = min(abs(f2 - min_f(zz)));
%         [~, idx_max] = min(abs(f2 - max_f(zz)));
% 
%         frf_ranged = frf(idx_min:idx_max);
%         real_ranged_frf = real_part_frf(idx_min:idx_max);
%         imag_ranged_frf = imag_part_frf(idx_min:idx_max);
% 
%         % secondo lui sono cerchi, cosa potrebbe cambiare per ottenere piu punti li? ragionare
%         % cosa sto plottando li?
%         % non è un problema di sensori e di misurazioni
%         % cosa dovrei cambiare nelle misurazioni per poter cambiare quella cosa li?
% 
%         figure
%         scatter(real_ranged_frf, imag_ranged_frf);
%         hold on
%         grid on
%     end
