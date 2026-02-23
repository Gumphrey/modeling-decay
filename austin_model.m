close all; clear all;

%% =======================
%  Global Parameters
%  =======================
fs        = 1000;        % Sampling rate (samples/s)
T         = 30;          % Total duration (s)
T_change  = 15;          % Change time (s) -- switch from initial to post condition
tau_mean  = 3.4;         % Time constant (s) for mean adaptation (offset)
tau_recovery = 5.0;      % Fixed recovery time constant (s) for gains when no stimulus
D         = 0.5;         % Duty cycle for both segments (unless customized)
mean0     = 0.5;         % Mean for initial segment
mean_post = 0.5;         % Mean for post segment (can differ)
vis_kPeriods = 3;        % Visibility smoothing window = k periods of the slowest frequency

% Pathway time constants used in transfer magnitudes
tauP_hp = 0.2;
tauP_lp = 0.025;         % Parvo LowPass
dBP_hp  = .5;             % Parvo steepness of drop off
dBP_lp = 3;
tau_hp = 0.038;         % Magno HP ~ 10 ms
tau_lp = 0.01;          % Magno LP ~ 40 ms
dBM_hp = 2;             % Magno HP steepness of drop off
dBM_lp = 2.9;           % Magno LP steepness of drop off

% Amplitudes used in transfer magnitudes
Ap = 1;          % Parvo transfer function amplitude
Am = 2;          % Magno transfer function amplitude
normalizeBool = false;     % Show transfer functions as normalized or not in Fig3

% Desensitization floors
AminP  = 0.0;
AminM  = 0.0;

% Frequency lists to test
%f_inits = [60];      % Hz
%f_posts = [0];      % Hz 
f_inits = [0, 2, 5, 15, 30, 60];      % Hz (rows)
f_posts = [0, 2, 5, 15, 30, 60];      % Hz (columns)

% Duty for post segment (can be customized per condition if needed)
D_post = 0.5;

%% =======================
%  Preallocate Results
%  =======================
nI = numel(f_inits);
nJ = numel(f_posts);
maxVIS_afterChange = nan(nI, nJ);

%% =======================
%  Run All Conditions
%  =======================
for ii = 1:nI
    f_init = f_inits(ii);

    for jj = 1:nJ
        f_post = f_posts(jj);

        % ---- Run one condition, update Figure 1, and return maximum post-change visibility
        maxVIS_afterChange(ii, jj) = simulate_and_plot_one_condition( ...
            fs, T, T_change, ...
            f_init, f_post, ...
            D, D_post, ...
            mean0, mean_post, ...
            tau_mean, tau_recovery, ...
            tauP_hp, tauP_lp, dBP_hp, dBP_lp, tau_hp, tau_lp, dBM_hp, dBM_lp, ...
            Ap, Am, AminP, AminM, ...
            vis_kPeriods);

        drawnow;
        % pause(0.05); % <- optional: slow down to watch each condition update
    end
end

%% =======================
%  Figure 2: Heatmap of max visibility after change (rows flipped)
%  =======================
figure(2); clf;

% Flip the matrix vertically so the top row corresponds to the highest initial frequency
Z = flipud(maxVIS_afterChange);

imagesc(Z);
axis image; colormap(bone); colorbar;

% Axis ticks/labels: columns are post freqs (unchanged), rows are initial freqs (reversed)
set(gca, 'XTick', 1:nJ, 'XTickLabel', string(f_posts), ...
         'YTick', 1:nI, 'YTickLabel', string(flip(f_inits)), ...
         'YDir','normal');

xlabel('Post frequency (Hz)');
ylabel('Initial frequency (Hz)');
title('Max visibility after change');

% (Optional) overlay numeric values in each cell
for ii = 1:nI
    for jj = 1:nJ
        text(jj, ii, sprintf('%.2f', Z(ii, jj)), ...
            'HorizontalAlignment','center', 'Color','w', 'FontSize',8, 'FontWeight','bold');
    end
end

%% =======================
%  Figure 3: Temporal tuning curves (Kelly CSF, Parvo, Magno) on log-log axes
%  =======================
normalize_tuning = normalizeBool;     % set false to show raw magnitudes
f_min_plot = 0.1;            % Hz (must be > 0 for log scale)
f_max_plot = 60;             % Hz
Nf = 2000;                   % dense sampling for smooth curves

f_vec = logspace(log10(f_min_plot), log10(f_max_plot), Nf);

% Compute magnitudes across frequency
K      = kelly_tCSF(f_vec);                      % Kelly temporal CSF
H_par  = abs(parvo_transfer(f_vec, tauP_hp, tauP_lp, Ap, dBP_hp, dBP_lp));    % Parvo |H_P(f)|
H_mag  = abs(magno_transfer(f_vec, tau_hp, tau_lp, Am, dBM_hp, dBM_lp));  % Magno |H_M(f)|

% Optional normalization for comparability
if normalize_tuning
    Kn     = K    ./ max(K    + eps);
    H_parn = H_par./ max(H_par+ eps);
    H_magn = H_mag./ max(H_mag+ eps);
end

figure(3); clf; set(gcf,'Color','w','Position',[1150 100 700 480]); hold on;

if normalize_tuning
    plot(f_vec, Kn,     'k-', 'LineWidth', 2);
    plot(f_vec, H_parn, 'g-', 'LineWidth', 2);
    plot(f_vec, H_magn, 'b-', 'LineWidth', 2);
    ylabel('Normalized amplitude (unitless)');
else
    plot(f_vec, K,     'k-', 'LineWidth', 2);
    plot(f_vec, H_par, 'g-', 'LineWidth', 2);
    plot(f_vec, H_mag, 'b-', 'LineWidth', 2);
    ylabel('Amplitude (unitless)');
end

set(gca, 'XScale','log', 'YScale','log');
xlim([f_min_plot, f_max_plot]);

% Nice log ticks (adjust as you like)
xticks([0.1 0.2 0.5 1 2 5 10 20 50]);
xticklabels(string([0.1 0.2 0.5 1 2 5 10 20 50]));
grid on; box on;

xlabel('Temporal frequency (Hz, log scale)');
title('Temporal Tuning (log–log): Kelly CSF, Parvo |H_P(f)|, Magno |H_M(f)|');
legend({'Kelly CSF','Parvo','Magno'}, 'Location','southwest');

%% =======================
%  Local Function: Simulate + Plot Condition
%  =======================
function maxVIS_post = simulate_and_plot_one_condition( ...
    fs, T, T_change, ...
    f_init, f_post, ...
    D, D_post, ...
    mean0, mean_post, ...
    tau_mean, tau_recovery, ...
    tauP_hp, tauP_lp, dBP_hp, dBP_lp, tau_hp, tau_lp, dBM_hp, dBM_lp,...
    Ap, Am, AminP, AminM, ...
    vis_kPeriods)
% SIMULATE_AND_PLOT_ONE_CONDITION
%   Simulates one experimental condition of a temporal-adaptation model.
%
%   The model has three stages:
%     1. Mean adaptation   – a slow LowPass filter tracks the running mean of the
%                            stimulus; subtracting it yields an AC-coupled
%                            "adaptation-only" signal x(t).
%     2. Desensitization   – two gain channels (Parvo & Magno) that decay
%                            toward a floor while a stimulus is present and
%                            recover toward 1 when the stimulus is absent.
%     3. Visibility output – the product of the two gains multiplied by
%                            x(t), then smoothed with a sliding window.
%
%   SPECIAL CASE: f_init = 0 or f_post = 0 is interpreted as "stimulus OFF"
%   (all zeros) for that segment. This replaces the need for post_mode='off'.
%
%   INPUTS (grouped by role):
%     Timing & sampling
%       fs          – sampling rate (Hz / samples per second)
%       T           – total simulation duration (s)
%       T_change    – time of the stimulus change (s); separates the
%                     "initial" (adapting) segment from the "post" segment
%
%     Stimulus frequencies
%       f_init      – square-wave frequency during the initial segment (Hz)
%                     0 Hz = stimulus OFF (all zeros) for this segment
%       f_post      – square-wave frequency during the post segment (Hz)
%                     0 Hz = stimulus OFF (all zeros) for this segment
%
%     Square-wave shape
%       D           – duty cycle for the initial segment (0–1; fraction of
%                     each period spent in the HIGH state)
%       D_post      – duty cycle for the post segment
%
%     Luminance / mean levels
%       mean0       – mean luminance around which the initial square wave is
%                     centred (only used when f_init > 0)
%       mean_post   – mean luminance for the post square wave
%                     (only used when f_post > 0)
%
%     Adaptation time constants
%       tau_mean    – time constant (s) of the 1st-order LP that tracks
%                     the running mean (controls how fast the offset adapts)
%       tau_recovery– time constant (s) at which Parvo/Magno gains recover
%                     back toward 1 when no stimulus is driving them
%
%     Pathway filter parameters (used to compute desensitization rates)
%       tauP        – Parvo pathway low-pass time constant (s); 
%       dbP         - Parvo drop off steepness (dB);
%       tau_hp      – Magno pathway high-pass time constant (s); ~10 ms
%       tau_lp      – Magno pathway low-pass time constant (s); ~40 ms
%       Ap          - Parvo transfer function amplitude
%       Am          - Magno transfer function amplitude
%
%     Desensitization floors
%       AminP       – minimum gain the Parvo channel can reach (0–1)
%       AminM       – minimum gain the Magno channel can reach (0–1)
%
%     Visibility smoothing
%       vis_kPeriods – number of periods of the slowest nonzero frequency in
%                      the stimulus pair used as the smoothing window length
%
%   OUTPUT:
%       maxVIS_post – the peak smoothed visibility value occurring AFTER the
%                     stimulus change at t_change (scalar)

    % =====================================================================
    %  1. TIME BASE AND SEGMENT MASKS
    % =====================================================================
    t  = (0:1/fs:T)';          % Column vector of sample times (s)
    dt = 1/fs;                  % Sample period (s)
    pre_idx  = (t < T_change);   % Logical mask: TRUE for samples in the initial segment
    post_idx = ~pre_idx;        % Logical mask: TRUE for samples in the post-change segment

    % =====================================================================
    %  2. BUILD THE INITIAL-SEGMENT SQUARE WAVE  (s1)
    % =====================================================================
    %  Special case: 0 Hz means stimulus OFF → all zeros.
    %  Otherwise, the Kelly temporal CSF gives a frequency-dependent
    %  peak-to-peak amplitude (App0). HIGH and LOW are placed around mean0
    %  weighted by the duty cycle so the time-average equals mean0.
    if f_init == 0
        % --- 0 Hz: stimulus is OFF for the initial segment ---
        s1 = zeros(size(t));
    else
        Amp0   = kelly_tCSF(f_init);              % Peak-to-peak amplitude from Kelly CSF
        high0  = mean0 + (1 - D)*Amp0;            % HIGH level of the initial square wave
        low0   = mean0 - D*Amp0;                  % LOW  level
        phase1 = mod(2*pi*f_init*t + 0, 2*pi);    % Instantaneous phase, set to 0 rad (wraps each cycle)
        sq1    = double(phase1 < 2*pi*D);         % Binary square wave: 1 during HIGH
        s1     = low0 + Amp0 * sq1;               % Scaled square wave
    end

    % =====================================================================
    %  3. BUILD THE POST-SEGMENT STIMULUS  (s2)
    % =====================================================================
    %  Same logic: 0 Hz → all zeros; otherwise build a square wave at f_post.
    if f_post == 0
        % --- 0 Hz: stimulus is OFF for the post segment ---
        s2 = zeros(size(t));
    else
        Amp_post = kelly_tCSF(f_post);                    % Peak-to-peak amplitude from Kelly CSF
        high2    = mean_post + (1 - D_post)*Amp_post;     % HIGH level of post square wave
        low2     = mean_post - D_post*Amp_post;           % LOW  level
        phase2   = mod(2*pi*f_post*t + 0, 2*pi);          % Instantaneous phase, set to 0 rad (wraps each cycle)
        sq2      = double(phase2 < 2*pi*D_post);
        s2       = low2 + Amp_post * sq2;
    end

    % =====================================================================
    %  4. COMBINE INTO THE FULL STIMULUS WAVEFORM  (s0)
    % =====================================================================
    s0 = s1;
    s0(post_idx) = s2(post_idx);

    % =====================================================================
    %  5. MEAN ADAPTATION (OFFSET TRACKING)
    % =====================================================================
    %  A causal 1st-order IIR low-pass filter tracks the running mean of s0.
    %  offset(n) = a_mean · offset(n-1) + (1-a_mean) · s0(n)
    a_mean = exp(-dt / tau_mean);
    offset = zeros(size(s0));
    offset(1) = 0;                      % Start unadapted
    for n = 2:numel(s0)
        offset(n) = a_mean*offset(n-1) + (1-a_mean)*s0(n);
    end

    % =====================================================================
    %  6. ADAPTATION-ONLY SIGNAL  (x)
    % =====================================================================
    x = s0 - offset;

    % =====================================================================
    %  7. INSTANTANEOUS EFFECTIVE FREQUENCY  (f_eff)
    % =====================================================================
    %  0 Hz in either segment signals "no stimulus" → gains recover.
    f_eff = zeros(size(t));
    f_eff(pre_idx)  = f_init;       % Could be 0 → recovery during initial segment
    f_eff(post_idx) = f_post;       % Could be 0 → recovery during post segment

    % =====================================================================
    %  8. DYNAMIC DESENSITIZATION GAINS  (A1 = Parvo, A2 = Magno)
    % =====================================================================
    %  f_eff > 0 → desensitize at a frequency-dependent rate
    %  f_eff == 0 → recover toward 1 at tau_recovery
    A1    = ones(size(t));
    A2    = ones(size(t));
    a_rec = exp(-dt / tau_recovery);

    for n = 2:numel(t)
        if f_eff(n) > 0
            tCSF_n = kelly_tCSF(f_eff(n));
            magP_n = abs(parvo_transfer(f_eff(n), tauP_hp, tauP_lp, Ap, dBP_hp, dBP_lp));
            magM_n = abs(magno_transfer(f_eff(n), tau_hp, tau_lp, Am, dBM_hp, dBM_lp));

            tau_desP_n = abs(5 / max(eps, (tCSF_n * magP_n)));
            tau_desM_n = abs(5 / max(eps, (tCSF_n * magM_n)));

            a_desP_n = exp(-dt / tau_desP_n);
            a_desM_n = exp(-dt / tau_desM_n);

            A1(n) = a_desP_n * A1(n-1) + (1 - a_desP_n) * AminP;
            A2(n) = a_desM_n * A2(n-1) + (1 - a_desM_n) * AminM;
        else
            % No stimulus → recover toward 1
            A1(n) = a_rec * A1(n-1) + (1 - a_rec) * 1;
            A2(n) = a_rec * A2(n-1) + (1 - a_rec) * 1;
        end
    end

    % Clamp gains to valid range [Amin, 1]
    A1 = max(min(A1, 1), AminP);
    A2 = max(min(A2, 1), AminM);

    % =====================================================================
    %  9. PERCEIVED SIGNAL AND VISIBILITY
    % =====================================================================
    G = A1 .* A2;           % Combined desensitization gain (multiplicative)
    y = G .* x;             % Perceived signal

    % --- Smoothing window for visibility ---
    %  Use the slowest NONZERO frequency in the pair for the window.
    %  If both are 0, fall back to a 1-sample window (no smoothing).
    nonzero_freqs = [f_init, f_post];
    nonzero_freqs = nonzero_freqs(nonzero_freqs > 0);
    if isempty(nonzero_freqs)
        f_min = Inf;    % Both segments OFF → 1-sample window (no smoothing)
    else
        f_min = min(nonzero_freqs);
    end
    winSamps = max(1, round(vis_kPeriods * fs / max(f_min, eps)));

    % Visibility as smoothed envelope via moving max of |y|.
    %  |y| captures both oscillatory flicker (crests and troughs as
    %  deviations from the adapted mean) AND slow one-sided deflections
    %  like afterimages. Using movmax rather than movmean avoids the
    %  underestimation that occurs when averaging across zero-crossings
    %  during flicker, while still faithfully tracking sustained
    %  afterimage deflections. The window provides temporal smoothing
    %  consistent with subjects reporting visibility ~once per second.
    VIS = movmax(abs(y), winSamps, 'Endpoints','shrink');

    % =====================================================================
    %  10. EXTRACT THE OUTPUT METRIC
    % =====================================================================
    maxVIS_post = max(VIS(post_idx));

    % =====================================================================
    %  11. PLOTTING (Figure 1) — six-panel time-series overview
    % =====================================================================
    figure(1); clf; set(gcf,'Color','w','Position',[100 100 1000 1200]);

    % --- Panel 1: Raw stimulus waveform ---
    subplot(6,1,1);
    plot(t, s0, 'g', 'LineWidth', 1); hold on
    xline(T_change,'--k','Change');
    ylabel('Stimulus s_0(t)');
    % Build a title that labels 0 Hz segments as "OFF"
    if f_init == 0
        init_str = 'OFF (0 Hz)';
    else
        init_str = sprintf('f=%.2f Hz (mean=%.2f, Amp=%.2f)', ...
                           f_init, mean0, kelly_tCSF(f_init));
    end
    if f_post == 0
        post_str = 'OFF (0 Hz)';
    else
        post_str = sprintf('f=%.2f Hz (mean=%.2f, Amp=%.2f)', ...
                           f_post, mean_post, kelly_tCSF(f_post));
    end
    title(sprintf('INIT: %s\nPOST: %s', init_str, post_str));
    grid on;

    % --- Panel 2: Mean-adaptation offset signal ---
    subplot(6,1,2);
    plot(t, offset, 'b', 'LineWidth', 1.2); hold on
    xline(T_change,'--k');
    ylabel('offset(t)');
    title(sprintf('ADAPTATION STATE (mean tracker): \\tau_{mean} = %.2f s', tau_mean));
    grid on;

    % --- Panel 3: Adaptation-only (AC-coupled) signal ---
    subplot(6,1,3);
    plot(t, x, 'm', 'LineWidth', 1); hold on
    xline(T_change,'--k'); yline(0,'k:');
    ylabel('s_0 - offset');
    title('PREDICTED VISUAL SIGNAL FOR ADAPTATION ONLY');
    grid on;

    % --- Panel 4: Parvo & Magno desensitization gains over time ---
    subplot(6,1,4);
    plot(t, A1, 'g', 'LineWidth', 1.2); hold on
    plot(t, A2, 'b', 'LineWidth', 1.2);
    xline(T_change,'--k');
    ylabel('A(t)');
    title(sprintf('DESENSITIZATION GAINS (\\tau_{rec}=%.1f s): init f=%g Hz, post f=%g Hz', ...
          tau_recovery, f_init, f_post));
    ylim([-0.05 1.05]); grid on;

    % --- Panel 5: Final perceived signal ---
    subplot(6,1,5);
    plot(t, y, 'r', 'LineWidth', 1); hold on
    xline(T_change,'--k'); yline(0,'k:');
    xlabel('Time (s)'); ylabel('Resultant y(t)');
    title('PREDICTED VISUAL SIGNAL WITH ADAPTATION AND DESENSITIZATION');
    grid on;

    % --- Panel 6: Smoothed visibility envelope ---
    subplot(6,1,6);
    plot(t, VIS, 'k', 'LineWidth', 1.5); hold on
    xline(T_change,'--k');
    ylabel(sprintf('|y(t)| (avg over %d periods)', vis_kPeriods));
    title(sprintf('VISIBILITY (smoothed). Max after change = %.3f', maxVIS_post));
    grid on;
end

%% =======================
%  Utility: Kelly temporal CSF amplitude scaling
%  =======================
function val = kelly_tCSF(f)
    % Kelly equation used previously for temporal CSF-based amplitude scaling
    val = (1/3.67879) * f .* exp(-f/10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local transfer functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = parvo_transfer(f, tau_hp, tau_lp, A, n_hp, n_lp)
% PARVO_TRANSFER  n-th order: NOW IDENTICAL TO MAGNO FUNCTION: band-pass:
% KEPT THE NAME FOR BACKWARDS COMPATIBILITY
% H(f) = A *
%   [ (i*w*tau_hp)/(1 + i*w*tau_hp) ]^n_hp *
%   [ 1/(1 + i*w*tau_lp) ]^n_lp
%
% n_hp controls low-frequency slope  (+20*n_hp dB/dec)
% n_lp controls high-frequency slope  (-20*n_lp dB/dec)

    if nargin < 5
        n_hp = 1;
    end
    if nargin < 6
        n_lp = 1;
    end

    w = 2*pi*f(:);

    H_hp = ((1i*w*tau_hp) ./ (1 + 1i*w*tau_hp)).^n_hp;
    H_lp = (1 ./ (1 + 1i*w*tau_lp)).^n_lp;

    H = A .* H_hp .* H_lp;
end

function H = magno_transfer(f, tau_hp, tau_lp, A, n_hp, n_lp)
% MAGNO_TRANSFER  Adjustable-order band-pass
%
% H(f) = A *
%   [ (i*w*tau_hp)/(1 + i*w*tau_hp) ]^n_hp *
%   [ 1/(1 + i*w*tau_lp) ]^n_lp
%
% n_hp controls low-frequency slope  (+20*n_hp dB/dec)
% n_lp controls high-frequency slope  (-20*n_lp dB/dec)

    if nargin < 5
        n_hp = 1;
    end
    if nargin < 6
        n_lp = 1;
    end

    w = 2*pi*f(:);

    H_hp = ((1i*w*tau_hp) ./ (1 + 1i*w*tau_hp)).^n_hp;
    H_lp = (1 ./ (1 + 1i*w*tau_lp)).^n_lp;

    H = A .* H_hp .* H_lp;
end