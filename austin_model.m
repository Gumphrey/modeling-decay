close all; clear all;

%% =======================
%  Global Parameters
%  =======================
fs        = 1000;        % Sampling rate (samples/s)
T         = 30;          % Total duration (s)
T_stim    = 15;          % Change time (s) -- switch from initial to post condition
tau_mean  = 3;           % Time constant (s) for mean adaptation (offset)
tau_recovery = 5.0;      % Fixed recovery time constant (s) for gains when no stimulus
D         = 0.5;         % Duty cycle for both segments (unless customized)
phi       = 0;           % Initial phase
phi_post  = 0;           % Post phase
mean0     = 0.5;         % Mean for initial segment
mean_post = 0.5;         % Mean for post segment (can differ)
I_bg      = 0;           % Background when post_mode = 'off'
vis_kPeriods = 3;        % Visibility smoothing window = k periods of the slowest frequency

% Pathway time constants used in transfer magnitudes
tauP   = 0.040;          % Parvo LP ~ 40 ms
tau_hp = 0.010;          % Magno HP ~ 10 ms
tau_lp = 0.040;          % Magno LP ~ 40 ms

% Desensitization floors
AminP  = 0.0;
AminM  = 0.0;

% Choose what happens after T_stim:
%   'off'     -> constant I_bg (heatmap not meaningful for 'off', but simulation still runs)
%   'square'  -> new square wave with its own frequency
post_mode  = 'square';

% Frequency lists to test
f_inits = [15];      % Hz (rows)
f_posts = [2];      % Hz (columns)
%f_inits = [0, 2, 5, 15, 30, 60];      % Hz (rows)
%f_posts = [0, 2, 5, 15, 30, 60];      % Hz (columns)

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
            fs, T, T_stim, ...
            f_init, f_post, ...
            D, D_post, phi, phi_post, ...
            mean0, mean_post, I_bg, post_mode, ...
            tau_mean, tau_recovery, ...
            tauP, tau_hp, tau_lp, ...
            AminP, AminM, ...
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
title(sprintf('Max visibility after change (post\\_mode = %s) — rows flipped', post_mode));

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
normalize_tuning = false;     % set false to show raw magnitudes
f_min_plot = 0.1;            % Hz (must be > 0 for log scale)
f_max_plot = 60;             % Hz
Nf = 2000;                   % dense sampling for smooth curves

f_vec = logspace(log10(f_min_plot), log10(f_max_plot), Nf);

% Compute magnitudes across frequency
K      = kelly_tCSF(f_vec);                      % Kelly temporal CSF
H_par  = abs(parvo_transfer(f_vec, tauP, 1));    % Parvo |H_P(f)|
H_mag  = abs(magno_transfer(f_vec, tau_hp, tau_lp, 6));  % Magno |H_M(f)|

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
    fs, T, T_stim, ...
    f_init, f_post, ...
    D, D_post, phi, phi_post, ...
    mean0, mean_post, I_bg, post_mode, ...
    tau_mean, tau_recovery, ...
    tauP, tau_hp, tau_lp, ...
    AminP, AminM, ...
    vis_kPeriods)

    % ----- Time base and indices
    t  = (0:1/fs:T)'; 
    dt = 1/fs;
    pre_idx  = (t < T_stim);
    post_idx = ~pre_idx;

    % ----- Initial segment: valued square s1 at f_init
    App0     = kelly_tCSF(f_init);    % peak-to-peak amplitude based on Kelly equation
    high0    = mean0 + (1 - D)*App0;
    low0     = mean0 - D*App0;
    phase1   = mod(2*pi*f_init*t + phi, 2*pi);
    sq1      = double(phase1 < 2*pi*D);
    s1       = low0 + App0 * sq1;

    % ----- Post segment: either OFF (I_bg) or a new square at f_post
    if strcmpi(post_mode, 'off')
        s2 = I_bg * ones(size(t));
    elseif strcmpi(post_mode, 'square')
        App_post = kelly_tCSF(f_post);
        high2    = mean_post + (1 - D_post)*App_post;
        low2     = mean_post - D_post*App_post;
        phase2   = mod(2*pi*f_post*t + phi_post, 2*pi);
        sq2      = double(phase2 < 2*pi*D_post);
        s2       = low2 + App_post * sq2;
    else
        error('post_mode must be ''off'' or ''square''.');
    end

    % ----- Combined stimulus
    s0 = s1;
    s0(post_idx) = s2(post_idx);

    % ----- Dynamic mean adaptation (offset) via 1st-order low-pass that tracks s0
    a_mean = exp(-dt / tau_mean);
    offset = zeros(size(s0));
    offset(1) = 0;   % start unadapted so adaptation-only initially matches stimulus
    for n = 2:numel(s0)
        offset(n) = a_mean*offset(n-1) + (1-a_mean)*s0(n);
    end

    % ----- Adaptation-only signal
    x = s0 - offset;

    % ----- Dynamic desensitization gains with time-varying frequency
    % Instantaneous frequency vector
    f_eff = zeros(size(t));
    f_eff(pre_idx) = f_init;
    if strcmpi(post_mode, 'square')
        f_eff(post_idx) = f_post;
    else
        f_eff(post_idx) = 0;      % off -> recover
    end

    A1 = ones(size(t));           % Parvo gain
    A2 = ones(size(t));           % Magno gain
    a_rec = exp(-dt / tau_recovery);

    for n = 2:numel(t)
        if f_eff(n) > 0
            tCSF_n = kelly_tCSF(f_eff(n));
            magP_n = abs(parvo_transfer(f_eff(n), tauP, 1));
            magM_n = abs(magno_transfer(f_eff(n), tau_hp, tau_lp, 4));
            tau_desP_n = abs(5 / max(eps, (tCSF_n * magP_n)));
            tau_desM_n = abs(5 / max(eps, (tCSF_n * magM_n)));
            a_desP_n = exp(-dt / tau_desP_n);
            a_desM_n = exp(-dt / tau_desM_n);
            A1(n) = a_desP_n * A1(n-1) + (1 - a_desP_n) * AminP;
            A2(n) = a_desM_n * A2(n-1) + (1 - a_desM_n) * AminM;
        else
            A1(n) = a_rec * A1(n-1) + (1 - a_rec) * 1;
            A2(n) = a_rec * A2(n-1) + (1 - a_rec) * 1;
        end
    end

    % Clamp gains to [Amin, 1]
    A1 = max(min(A1, 1), AminP);
    A2 = max(min(A2, 1), AminM);

    % ----- Perceived signal and visibility
    G  = A1 .* A2; % TODO: SHOULD THIS BE AVERAGE OR ADDITIVE? RIGHT NOW
                    % IT IS ZERO IF ONE PATHWAY IS ZERO
    y  = G .* x; % y is visibility in second to last curve

    % Instantaneous visibility and smoothed visibility over k periods of the slowest frequency in pair
    if strcmpi(post_mode,'square')
        f_min = min([f_init, f_post]);
    else
        f_min = max(eps, f_init);  % only initial segment has flicker
    end
    winSamps = max(1, round(vis_kPeriods * fs / max(f_min, eps)));
    VIS_inst = abs(y);
    VIS      = movmean(VIS_inst, winSamps, 'Endpoints','shrink');

    % ----- Max visibility AFTER the change
    maxVIS_post = max(VIS(post_idx));

    % ----- Plot (Figure 1) — update for this condition
    figure(1); clf; set(gcf,'Color','w','Position',[100 100 1000 1200]);

    subplot(6,1,1);
    plot(t, s0, 'g', 'LineWidth', 1); hold on
    xline(T_stim,'--k','Change');
    ylabel('Stimulus s_0(t)');
    if strcmpi(post_mode,'off')
        yline(I_bg,'--','Color',[0.4 0.4 0.4],'DisplayName','I_{bg}');
        ttl1 = sprintf(['INIT: f=%.2f Hz (mean=%.2f, App=%.2f)\n' ...
                        'POST: OFF to I_{bg}=%.2f'], ...
                        f_init, mean0, kelly_tCSF(f_init), I_bg);
    else
        ttl1 = sprintf(['INIT: f=%.2f Hz (mean=%.2f, App=%.2f)\n' ...
                        'POST: f=%.2f Hz (mean=%.2f, App=%.2f)'], ...
                        f_init, mean0, kelly_tCSF(f_init), ...
                        f_post, mean_post, kelly_tCSF(f_post));
    end
    title(ttl1); grid on;

    subplot(6,1,2);
    plot(t, offset, 'b', 'LineWidth', 1.2); hold on
    xline(T_stim,'--k');
    ylabel('offset(t)');
    title(sprintf('ADAPTATION STATE (mean tracker): \\tau_{mean} = %.2f s', tau_mean));
    grid on;

    subplot(6,1,3);
    plot(t, x, 'm', 'LineWidth', 1); hold on
    xline(T_stim,'--k'); yline(0,'k:');
    ylabel('s_0 - offset');
    title('PREDICTED VISUAL SIGNAL FOR ADAPTATION ONLY');
    grid on;

    subplot(6,1,4);
    plot(t, A1, 'g', 'LineWidth', 1.2); hold on
    plot(t, A2, 'b', 'LineWidth', 1.2);
    xline(T_stim,'--k');
    ylabel('A(t)');
    if strcmpi(post_mode,'off')
        title(sprintf('DESENSITIZATION GAINS (recover to 1, \\tau_{rec}=%.1f s): pre f=%.2f Hz', tau_recovery, f_init));
    else
        title(sprintf('DESENSITIZATION GAINS (recover to 1, \\tau_{rec}=%.1f s): pre f=%.2f Hz, post f=%.2f Hz', ...
              tau_recovery, f_init, f_post));
    end
    ylim([-0.05 1.05]); grid on;

    subplot(6,1,5);
    plot(t, y, 'r', 'LineWidth', 1); hold on
    xline(T_stim,'--k'); yline(0,'k:');
    xlabel('Time (s)'); ylabel('Resultant y(t)');
    title('PREDICTED VISUAL SIGNAL WITH ADAPTATION AND DESENSITIZATION');
    grid on;

    subplot(6,1,6);
    plot(t, VIS, 'k', 'LineWidth', 1.5); hold on
    xline(T_stim,'--k');
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
function H = parvo_transfer(f, tau, A)
% PARVO_TRANSFER  First-order low-pass: H(f) = A / (1 + i*2*pi*f*tau)
    w = 2*pi*f(:);
    H = A ./ (1 + 1i*w*tau);
end

function H = magno_transfer(f, tau_hp, tau_lp, A)
% MAGNO_TRANSFER  Band-pass via HP-LP cascade:
%   H(f) = A * [ (i*w*tau_hp)/(1 + i*w*tau_hp) ] * [ 1/(1 + i*w*tau_lp) ]
    w    = 2*pi*f(:);
    H_hp = (1i*w*tau_hp) ./ (1 + 1i*w*tau_hp);  % high-pass
    H_lp = 1 ./ (1 + 1i*w*tau_lp);              % low-pass
    H    = A .* H_hp .* H_lp;
end