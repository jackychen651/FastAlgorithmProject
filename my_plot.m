% Clear workspace and figures
clear;
close all;
clc;

% Define the data
N = [256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 2^24];

% Times corresponding to each N
T1_build = [0.007100, 0.073868, 1.112, 17.46, NaN, NaN, NaN, NaN];
T1_solve = [0.002495, 0.050911, 2.041, 126.9, NaN, NaN, NaN, NaN];
T1 = T1_build+T1_solve;
T2       = [1.310, 1.307, 1.287, 1.360, 1.589, 2.373, 6.846, 25.62, 164.545212];

% Create a new figure
%figure;
%hold on;
%grid on;

% Plot T1_build (exclude NaN values)
valid_T1_build = ~isnan(T1_build);
% loglog(N(valid_T1_build), T1_build(valid_T1_build), 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'T_1^{build}');

% Plot T1_solve (exclude NaN values)
valid_T1_solve = ~isnan(T1_solve);
% loglog(N(valid_T1_solve), T1_solve(valid_T1_solve), 's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'T_1^{solve}');

% Plot T2 (all values available)
% Ensure that T2 has no zero or negative values for log scaling

valid_T1= valid_T1_build;
valid_T2 = T2 > 0;

loglog(N(valid_T1), T1(valid_T1),'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'T_1 of Vanilla Solver');
hold on;
loglog(N(valid_T2), T2(valid_T2),'s-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'T_2 of Iterative Solver');

% Choose the base point for C1 and C2
base_N = N(1);      % N = 256
base_T2 = T2(1);    % T2 = 1.310 seconds

% Calculate C1 and C2 for T2
C1 = 0.00000000004;           % For O(N^3)
C2 = 0.0000005; % For O(N logN)

% Create reference lines using C1 and C2
ref_N = N;
ref_ONlogN = C2 * ref_N .* log(ref_N);    % O(N logN)
ref_ON3    = C1 * ref_N.^3;               % O(N^3)

% Plot reference lines
loglog(ref_N, ref_ONlogN, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference O(NlogN)');
%loglog(ref_N, ref_O_N2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference O(N^2)');
loglog(ref_N, ref_ON3, 'k-.', 'LineWidth', 1.5, 'DisplayName', 'Reference O(N^3)');

legend('Location', 'northwest', 'FontSize', 12);
hold off;

filename = sprintf('Complexity.pdf');
set(gcf, 'PaperUnits', 'centimeters');
paperWidth = 15;
paperHeight = 10;
set(gcf, 'PaperSize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0 0 paperWidth paperHeight]);
saveas(gcf, filename);