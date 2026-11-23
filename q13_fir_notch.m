%% BEC504 - DIGITAL SIGNAL PROCESSING : EXPERIENTIAL LEARNING (V Sem, A Sec, AY 2026-27)
%  Question 13 : Suppression of 50 Hz Power-Line Interference Using an FIR Notch Filter
%  Team 2      : Gagan Manjunath (4NI24EC047), Manoj Girish Bhat (4NI24EC071), B Srivatsa Joshi (4NI24EC026)
%
%  Real input  : PhysioNet ECG-ID Database, Person_25, rec_1 (ECG lead I, raw channel)
%                fs = 500 Hz, 20 s, 12-bit, gain 200 ADU/mV.  The raw channel contains
%                visible 50 Hz mains interference (checked in Section 2).
%  Source link : https://physionet.org/content/ecgiddb/1.0.0/Person_25/
%
%  Run this file from its own folder (q13_fir_notch.m); it needs data/Person25_rec_1.dat
%  and writes figures + text results into ./results .
%  Needs: MATLAB + Signal Processing Toolbox (fir1, freqz, findpeaks, grpdelay).

clear; clc; close all;
set(0,'DefaultFigureVisible','off');                 % figures are saved to PNG files
here   = fileparts(mfilename('fullpath'));
outDir = fullfile(here,'results');
if ~exist(outDir,'dir'), mkdir(outDir); end
logf = fullfile(outDir,'console_output.txt');
if exist(logf,'file'), delete(logf); end
diary(logf); diary on;

%% 1. LOAD THE REAL ECG RECORD ------------------------------------------------
fs   = 500;                 % sampling frequency from the .hea header (Hz)
f0   = 50;                  % mains frequency to be suppressed (Hz)
gain = 200;                 % ADU per mV (from the .hea header)
fid  = fopen(fullfile(here,'data','Person25_rec_1.dat'),'r');
raw  = fread(fid,[2 Inf],'int16=>double');          % row 1 = raw ECG, row 2 = dataset's own filtered ECG
fclose(fid);
x    = raw(1,:).'/gain;     % raw ECG in mV
x    = x - mean(x);         % remove DC offset
L    = numel(x);
t    = (0:L-1).'/fs;
fprintf('Record: %d samples, fs = %d Hz, duration = %.1f s\n', L, fs, L/fs);
fprintf('Digital frequency of the mains: w0 = 2*pi*%d/%d = %.4f rad/sample (= 0.2*pi)\n\n', f0, fs, 2*pi*f0/fs);

%% 2. SPECTRUM OF THE INPUT (is the 50 Hz interference really there?) ---------
[fx, Ax] = amp_spectrum(x, fs);
q50 = find(fx>=47 & fx<=53);                           % search around the nominal mains frequency
[pk50, ip] = max(Ax(q50));  fpk = fx(q50(ip));         % mains tones are rarely exactly 50.000 Hz
nb = fx>=40 & fx<=46;                                  % neighbouring noise floor
fprintf('Raw ECG: interference peak = %.4f mV at %.2f Hz (nominal 50 Hz), median floor (40-46 Hz) = %.4f mV  (ratio %.1f)\n\n', ...
        pk50, fpk, median(Ax(nb)), pk50/median(Ax(nb)));

%% 3. (a) DESIGN THE FIR NOTCH FILTER (windowed-sinc band-stop, Hamming window) -
N  = 200;                   % filter order  (N+1 = 201 taps, Type-I linear phase)
BW = 10;                    % stop-band width of the notch (Hz): 45 - 55 Hz
M  = N/2;                   % group delay in samples
h  = design_notch(N, BW, fs, f0);
hf = fir1(N, [f0-BW/2 f0+BW/2]/(fs/2), 'stop', hamming(N+1));   % MATLAB reference design
fprintf('Design: N = %d (%d taps), stop-band %g-%g Hz, Hamming window, group delay = %d samples = %.0f ms\n', ...
        N, N+1, f0-BW/2, f0+BW/2, M, 1000*M/fs);
fprintf('Max |h_manual - h_fir1| = %.3e  (own design formula agrees with fir1)\n', max(abs(h-hf(:))));
fprintf('Symmetry check max|h[n]-h[N-n]| = %.1e, sum(h) = %.6f (DC gain)\n\n', max(abs(h-flipud(h))), sum(h));

%% 4. (b) MAGNITUDE RESPONSE NEAR 50 Hz ---------------------------------------
[H, fr] = freqz(h, 1, 2^16, fs);
magdB   = 20*log10(abs(H)+eps);
att50   = magdB(find(fr>=f0,1));
pb      = magdB(fr>=0.5 & fr<=40);
fprintf('|H(50 Hz)| = %.2f dB   (attenuation %.1f dB, linear gain %.4f)\n', att50, -att50, 10^(att50/20));
fprintf('|H| at the measured interference peak (%.2f Hz) = %.1f dB\n', fpk, magdB(find(fr>=fpk,1)));
fprintf('Gain over the ECG band 0.5-40 Hz : min %.4f dB, max %.4f dB\n', min(pb), max(pb));
fprintf('-3 dB notch width = %.2f Hz, -20 dB notch width = %.2f Hz, -40 dB width = %.2f Hz\n\n', ...
        stopwidth(h,fs,-3), stopwidth(h,fs,-20), stopwidth(h,fs,-40));

fig = figure('Position',[50 50 1200 800]);
subplot(2,2,1); stem(0:N, h, 'Marker','none'); grid on;
title('Impulse response h[n] (201 taps, symmetric)'); xlabel('n'); ylabel('h[n]');
subplot(2,2,2); plot(fr, magdB, 'LineWidth',1.2); grid on; ylim([-80 5]);
title('Magnitude response (full band)'); xlabel('Frequency (Hz)'); ylabel('|H| (dB)');
subplot(2,2,3); plot(fr, magdB, 'LineWidth',1.5); grid on; xlim([35 65]); ylim([-80 5]); hold on;
xline(f0,'r--','50 Hz'); yline(-3,'k:'); yline(-20,'k:');
title(sprintf('Zoom near 50 Hz: attenuation at 50 Hz = %.1f dB', -att50)); xlabel('Frequency (Hz)'); ylabel('|H| (dB)');
subplot(2,2,4); plot(fr, abs(H), 'LineWidth',1.5); grid on; xlim([0 100]); ylim([0 1.1]); hold on; xline(f0,'r--');
title('Linear magnitude, 0-100 Hz'); xlabel('Frequency (Hz)'); ylabel('|H|');
sgtitle('Fig. 1  FIR notch filter design (N = 200, stop-band 45-55 Hz, Hamming)');
savefig_png(fig, fullfile(outDir,'fig1_filter_design.png'));

%% 5. (c) FILTER THE ECG : manual FIR difference equation vs filter() ---------
y_man = fir_manual(h, x);                     % y[n] = sum_k h[k] x[n-k]  (nested loops)
y     = filter(h, 1, x);                      % MATLAB built-in
fprintf('Max |y_manual - filter()| = %.3e  (manual difference-equation implementation verified)\n', max(abs(y_man-y)));
ya    = [y(M+1:end); zeros(M,1)];             % delay-compensated output (linear phase => pure delay of M samples)

% Ideal reference: remove 50 Hz +/- 1 Hz directly in the DFT domain (not a causal filter, used only as a yardstick)
X = fft(x); kf = (0:L-1).'*fs/L; kf(kf>fs/2) = kf(kf>fs/2)-fs;
X(abs(abs(kf)-f0)<=1) = 0;  xref = real(ifft(X));

idx = 401:L-400;                              % common segment (free of start/end transients for every N used below)
[f1, A1] = amp_spectrum(x(idx),  fs);         % before
[f2, A2] = amp_spectrum(ya(idx), fs);         % after
a50b = max(A1(abs(f1-f0)<=1));  a50a = max(A2(abs(f2-f0)<=1));
fprintf('\n50 Hz (mains peak, +/-1 Hz) spectral component  : before = %.5f mV, after = %.6f mV  ->  reduction = %.1f dB\n', a50b, a50a, 20*log10(a50b/a50a));
rms_b = rms(bandpass_ms(x(idx),fs,f0)); rms_a = rms(bandpass_ms(ya(idx),fs,f0));
fprintf('Power in the 49-51 Hz band : before = %.3e mV^2, after = %.3e mV^2 -> %.1f dB lower\n\n', rms_b^2, rms_a^2, 20*log10(rms_b/rms_a));

% R-peak detection on the clean reference
[~,locR] = findpeaks(xref, 'MinPeakHeight', 0.5*max(xref(idx)), 'MinPeakDistance', round(0.4*fs));
locR = locR(locR>=idx(1) & locR<=idx(end));
fprintf('Detected %d R peaks in the common segment, mean heart rate = %.1f bpm\n', numel(locR), 60*fs/mean(diff(locR)));
fprintf('R-peak amplitude: reference %.4f mV, filtered %.4f mV (change %.2f %%)\n', mean(xref(locR)), mean(ya(locR)), 100*(mean(ya(locR))/mean(xref(locR))-1));
nr = norm(ya(idx)-xref(idx))/norm(xref(idx));
fprintf('Normalised RMS difference (filtered vs ideal-cleaned reference) = %.2f %%\n', 100*nr);
fprintf('Correlation coefficient with reference = %.5f\n\n', corr(ya(idx), xref(idx)));

% ---- Figure 2 : time domain
k0 = find(t>=4,1); k1 = find(t>=9,1);
fig = figure('Position',[50 50 1200 800]);
subplot(3,1,1); plot(t(k0:k1), x(k0:k1),'Color',[0.75 0.1 0.1]); grid on; xlim([4 9]);
title('Raw ECG with 50 Hz interference'); ylabel('mV');
subplot(3,1,2); plot(t(k0:k1), ya(k0:k1),'Color',[0.1 0.3 0.7]); grid on; xlim([4 9]);
title('After FIR notch filter (delay-compensated)'); ylabel('mV');
kk = find(t>=locR(5)/fs-0.25,1):find(t>=locR(5)/fs+0.35,1);
subplot(3,1,3); plot(t(kk),x(kk),'Color',[0.75 0.1 0.1],'LineWidth',0.8); hold on;
plot(t(kk),ya(kk),'Color',[0.1 0.3 0.7],'LineWidth',1.5); plot(t(kk),xref(kk),'k--','LineWidth',1);
grid on; legend('raw','filtered','ideal reference','Location','northeast'); xlabel('Time (s)'); ylabel('mV');
title('One heart beat (zoom): QRS shape is preserved');
sgtitle('Fig. 2  Time-domain result'); savefig_png(fig, fullfile(outDir,'fig2_time_domain.png'));

% ---- Figure 3 : spectra
fig = figure('Position',[50 50 1200 650]);
subplot(1,2,1); plot(f1, 20*log10(A1+eps),'Color',[0.75 0.1 0.1]); hold on; plot(f2, 20*log10(A2+eps),'Color',[0.1 0.3 0.7]);
grid on; xlim([0 150]); xline(f0,'k--'); legend('before','after','50 Hz','Location','northeast');
xlabel('Frequency (Hz)'); ylabel('Amplitude (dB re 1 mV)'); title('ECG spectrum, 0-150 Hz');
subplot(1,2,2); plot(f1, 20*log10(A1+eps),'Color',[0.75 0.1 0.1],'LineWidth',1.3); hold on; plot(f2, 20*log10(A2+eps),'Color',[0.1 0.3 0.7],'LineWidth',1.3);
grid on; xlim([40 60]); xline(f0,'k--'); legend('before','after','Location','northeast');
xlabel('Frequency (Hz)'); ylabel('Amplitude (dB re 1 mV)'); title(sprintf('Zoom at 50 Hz: %.1f dB reduction', 20*log10(a50b/a50a)));
sgtitle('Fig. 3  Spectrum before and after filtering (Hann window, 9200 samples)'); savefig_png(fig, fullfile(outDir,'fig3_spectra.png'));

%% 6. MANUAL CALCULATION DATA (for the report) --------------------------------
fprintf('================ DATA FOR MANUAL CALCULATIONS ================\n');
wc1 = 2*pi*(f0-BW/2)/fs;  wc2 = 2*pi*(f0+BW/2)/fs;
fprintf('wc1 = 2*pi*45/500 = %.5f rad = %.2f*pi ; wc2 = 2*pi*55/500 = %.5f rad = %.2f*pi\n', wc1, wc1/pi, wc2, wc2/pi);
fprintf('Centre tap (before scaling): hd[M] = 1 - (wc2-wc1)/pi = %.4f\n', 1-(wc2-wc1)/pi);
for m = 1:3
    hd = -(sin(wc2*m)-sin(wc1*m))/(pi*m);  w = 0.54-0.46*cos(2*pi*(M+m)/N);
    fprintf('m=%d: hd = -(sin(%.4f)-sin(%.4f))/(pi*%d) = %.6f ; Hamming w[%d] = %.6f ; hd*w = %.6f ; final h[%d] = %.6f (MATLAB)\n', ...
        m, wc2*m, wc1*m, m, hd, M+m, w, hd*w, M+m, h(M+m+1));
end
fprintf('h[0..4]      = %s\n', mat2str(h(1:5).',6));
fprintf('h[M]=h[100]  = %.6f\n', h(M+1));
fprintf('x[0..4] (mV) = %s\n', mat2str(x(1:5).',6));
fprintf('y[0] = h[0]x[0]                          = %.8f ; MATLAB filter() = %.8f\n', h(1)*x(1), y(1));
fprintf('y[1] = h[0]x[1] + h[1]x[0]               = %.8f ; MATLAB filter() = %.8f\n', h(1)*x(2)+h(2)*x(1), y(2));
fprintf('y[2] = h[0]x[2] + h[1]x[1] + h[2]x[0]    = %.8f ; MATLAB filter() = %.8f\n', h(1)*x(3)+h(2)*x(2)+h(3)*x(1), y(3));
% a steady-state sample: show the large terms of the 201-term sum
n0 = 1000; terms = h(1:N+1).*x(n0+1:-1:n0+1-N);
fprintf('Steady-state sample n = %d: y[n] = sum of 201 products = %.8f ; filter() = %.8f\n', n0, sum(terms), y(n0+1));
[~,ord] = sort(abs(terms),'descend');
fprintf('   largest products (k, h[k], x[n-k], h[k]x[n-k]):\n');
for q = 1:5, k = ord(q)-1; fprintf('   k=%3d  h=%+.6f  x=%+.5f  prod=%+.6f\n', k, h(k+1), x(n0+1-k), terms(ord(q))); end

% Simple 3-tap FIR notch (zeros on the unit circle at +/- w0) for hand calculation and comparison
w0  = 2*pi*f0/fs;  h3 = [1 -2*cos(w0) 1]/(2-2*cos(w0));
y3  = filter(h3,1,x);
fprintf('\n3-tap notch: h3 = [1, -2cos(w0), 1]/(2-2cos(w0)) = %s ; |H3(50 Hz)| = %.3e\n', mat2str(h3,6), abs(polyval(h3,exp(1j*w0))));
n1 = 1000; fprintf('y3[%d] = h3[0]x[%d]+h3[1]x[%d]+h3[2]x[%d] = %.8f (MATLAB %.8f)\n', n1,n1,n1-1,n1-2, ...
    h3(1)*x(n1+1)+h3(2)*x(n1)+h3(3)*x(n1-1), y3(n1+1));
fprintf('x[%d..%d] = %s\n\n', n1-2, n1, mat2str(x(n1-1:n1+1).',6));

%% 7. (d) ATTENUATION AT 50 Hz AND DISTORTION OF THE USEFUL SIGNAL --------------
% quantities per design: attenuation at 50 Hz, passband ripple 0.5-40 Hz, delay, 50 Hz reduction, NRMSE, R-peak change
%% 8. (e) COMPARE FILTER ORDERS (BW = 10 Hz) AND NOTCH BANDWIDTHS (N = 200) ------
orders = [50 100 200 300 400];
bws    = [2 4 6 10 16];
Tn = sweep_table(orders, repmat(10,1,numel(orders)), x, xref, locR, idx, fs, f0);
Tb = sweep_table(repmat(200,1,numel(bws)), bws, x, xref, locR, idx, fs, f0);
fprintf('================ EFFECT OF FILTER ORDER (stop-band width 10 Hz) ================\n'); disp(Tn);
fprintf('================ EFFECT OF NOTCH BANDWIDTH (N = 200) ================\n'); disp(Tb);
% extra: 3-tap notch in the same metrics
T3 = metrics_row(h3(:), 2, 0, x, xref, locR, idx, fs, f0);
fprintf('3-tap notch (N = 2): att50 = %.1f dB, ripple 0.5-40 Hz = [%.2f, %.2f] dB, 50 Hz reduction = %.1f dB, NRMSE = %.1f %%, R-peak change = %.1f %%\n\n', ...
    T3.Att50_dB, T3.PB_min_dB, T3.PB_max_dB, T3.Red50_dB, T3.NRMSE_pct, T3.Rchange_pct);
writetable(Tn, fullfile(outDir,'table_order_comparison.csv'));
writetable(Tb, fullfile(outDir,'table_bandwidth_comparison.csv'));

col = lines(5);
fig = figure('Position',[50 50 1250 820]);
subplot(2,2,1); hold on;
for q=1:numel(orders), hq = design_notch(orders(q),10,fs,f0); [Hq,fq]=freqz(hq,1,2^15,fs); plot(fq,20*log10(abs(Hq)+eps),'Color',col(q,:),'LineWidth',1.3); end
grid on; xlim([35 65]); ylim([-80 5]); xline(f0,'k--'); legend(arrayfun(@(o) sprintf('N=%d',o),orders,'Un',0),'Location','southwest');
xlabel('Frequency (Hz)'); ylabel('|H| (dB)'); title('Effect of order (stop-band 45-55 Hz)');
subplot(2,2,2); hold on;
for q=1:numel(bws), hq = design_notch(200,bws(q),fs,f0); [Hq,fq]=freqz(hq,1,2^15,fs); plot(fq,20*log10(abs(Hq)+eps),'Color',col(q,:),'LineWidth',1.3); end
grid on; xlim([35 65]); ylim([-80 5]); xline(f0,'k--'); legend(arrayfun(@(b) sprintf('BW=%d Hz',b),bws,'Un',0),'Location','southwest');
xlabel('Frequency (Hz)'); ylabel('|H| (dB)'); title('Effect of notch bandwidth (N = 200)');
subplot(2,2,3); yyaxis left; plot(Tn.N, Tn.Att50_dB,'o-','LineWidth',1.5); ylabel('Attenuation at 50 Hz (dB)');
yyaxis right; plot(Tn.N, Tn.Delay_ms,'s--','LineWidth',1.5); ylabel('Delay (ms)'); grid on; xlabel('Order N'); title('Order: selectivity vs delay');
subplot(2,2,4); yyaxis left; plot(Tb.BW_Hz, Tb.Att50_dB,'o-','LineWidth',1.5); ylabel('Attenuation at 50 Hz (dB)');
yyaxis right; plot(Tb.BW_Hz, Tb.NRMSE_pct,'s--','LineWidth',1.5); ylabel('NRMSE vs reference (%)'); grid on; xlabel('Stop-band width (Hz)'); title('Bandwidth: attenuation vs distortion');
sgtitle('Fig. 4  Comparison of filter orders and notch bandwidths'); savefig_png(fig, fullfile(outDir,'fig4_order_bandwidth.png'));

% Figure 5: time-domain effect of the three extremes on one beat + 3-tap
kk = find(t>=locR(5)/fs-0.2,1):find(t>=locR(5)/fs+0.3,1);
fig = figure('Position',[50 50 1200 600]); hold on;
plot(t(kk), xref(kk),'k','LineWidth',2);
cfg = {[200 10],[200 16],[400 16],[100 16]}; cc = lines(4);
for q = 1:numel(cfg)
    hq = design_notch(cfg{q}(1),cfg{q}(2),fs,f0); yq = filter(hq,1,x); Mq = cfg{q}(1)/2; yq = [yq(Mq+1:end); zeros(Mq,1)];
    plot(t(kk), yq(kk),'Color',cc(q,:),'LineWidth',1.1);
end
y3a = [y3(2:end);0]; plot(t(kk), y3a(kk),'m:','LineWidth',1.3);
grid on; xlabel('Time (s)'); ylabel('mV'); legend('ideal reference','N=200, BW=10','N=200, BW=16','N=400, BW=16','N=100, BW=16','3-tap notch','Location','northeast');
title('Fig. 5  QRS complex after different notch filters (delay-compensated)'); savefig_png(fig, fullfile(outDir,'fig5_qrs_distortion.png'));

% Figure 6: group delay + phase of main filter
[gd, fg] = grpdelay(h,1,4096,fs);
fig = figure('Position',[50 50 1000 400]);
subplot(1,2,1); plot(fg, gd,'LineWidth',1.5); grid on; ylim([0 M+20]); xlabel('Frequency (Hz)'); ylabel('Samples'); title(sprintf('Group delay = %d samples (%.0f ms)',M,1000*M/fs));
subplot(1,2,2); zplane(h,1); title('Pole-zero plot (zeros only, FIR)');
sgtitle('Fig. 6  Linear-phase property'); savefig_png(fig, fullfile(outDir,'fig6_delay_zplane.png'));

save(fullfile(outDir,'q13_workspace.mat'),'h','h3','x','y','ya','xref','fs','Tn','Tb');
fprintf('Done. Figures, tables and this log are in: %s\n', outDir);
diary off;

%% ============================ LOCAL FUNCTIONS =================================
function h = design_notch(N, BW, fs, f0)
% Windowed-sinc band-stop FIR:  h[n] = w[n]*( delta[m] - (sin(wc2 m)-sin(wc1 m))/(pi m) ),  m = n - N/2
    M = N/2;  n = (0:N).';  m = n - M;
    wc1 = 2*pi*(f0-BW/2)/fs;  wc2 = 2*pi*(f0+BW/2)/fs;
    hd = -(sin(wc2*m)-sin(wc1*m))./(pi*m);
    hd(M+1) = 1 - (wc2-wc1)/pi;
    w  = 0.54 - 0.46*cos(2*pi*n/N);          % Hamming window
    h  = hd.*w;
    h  = h/sum(h);                           % unity gain at DC (same normalisation as fir1)
end

function y = fir_manual(h, x)
% Direct-form FIR difference equation  y[n] = sum_{k=0}^{N} h[k] x[n-k]  with nested loops (zero initial conditions)
    L = numel(x);  K = numel(h);  y = zeros(L,1);
    for n = 1:L
        acc = 0;
        for k = 1:min(n,K)
            acc = acc + h(k)*x(n-k+1);
        end
        y(n) = acc;
    end
end

function [f, A] = amp_spectrum(s, fs)
% single-sided amplitude spectrum with a Hann window (amplitude-corrected)
    s = s(:);  L = numel(s);  w = hann(L);
    S = fft((s-mean(s)).*w);
    A = 2*abs(S(1:floor(L/2)+1))/sum(w);
    f = (0:floor(L/2)).'*fs/L;
end

function r = bandpass_ms(s, fs, f0)
% 49-51 Hz band component via DFT masking (used only to measure residual interference power)
    L = numel(s);  S = fft(s(:)-mean(s));  k = (0:L-1).'*fs/L;  k(k>fs/2) = k(k>fs/2)-fs;
    S(abs(abs(k)-f0)>1) = 0;  r = real(ifft(S));
end

function bw = stopwidth(h, fs, thr)
    [H,f] = freqz(h,1,2^16,fs);  m = 20*log10(abs(H)+eps);  q = find(m<thr & f>20 & f<80);
    if isempty(q), bw = NaN; else, bw = f(q(end))-f(q(1)); end
end

function T = sweep_table(Ns, BWs, x, xref, locR, idx, fs, f0)
    T = table();
    for q = 1:numel(Ns)
        T = [T; metrics_row(design_notch(Ns(q),BWs(q),fs,f0), Ns(q), BWs(q), x, xref, locR, idx, fs, f0)]; %#ok<AGROW>
    end
end

function R = metrics_row(h, N, BW, x, xref, locR, idx, fs, f0)
    M  = N/2;  [H,f] = freqz(h,1,2^16,fs);  m = 20*log10(abs(H)+eps);
    y  = filter(h,1,x);  if N>=2, ya = [y(M+1:end); zeros(M,1)]; end
    [f1,A1] = amp_spectrum(x(idx),fs);  [f2,A2] = amp_spectrum(ya(idx),fs);
    a_b = max(A1(abs(f1-f0)<=1));  a_a = max(A2(abs(f2-f0)<=1));
    pb  = m(f>=0.5 & f<=40);
    R = table(N, BW, -m(find(f>=f0,1)), min(pb), max(pb), 1000*M/fs, 20*log10(a_b/a_a), ...
              100*norm(ya(idx)-xref(idx))/norm(xref(idx)), 100*(mean(ya(locR))/mean(xref(locR))-1), ...
              'VariableNames',{'N','BW_Hz','Att50_dB','PB_min_dB','PB_max_dB','Delay_ms','Red50_dB','NRMSE_pct','Rchange_pct'});
end

function savefig_png(fig, file)
    exportgraphics(fig, file, 'Resolution', 170);  close(fig);
end
