clear;
close all;
clc;

N = [256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 2^24];

T1_build = [0.007100, 0.073868, 1.112, 17.46, NaN, NaN, NaN, NaN];
T1_solve = [0.002495, 0.050911, 2.041, 126.9, NaN, NaN, NaN, NaN];
T1 = T1_build+T1_solve;
T2       = [1.310, 1.307, 1.287, 1.360, 1.589, 2.373, 6.846, 25.62, 164.545212];

valid_T1_build = ~isnan(T1_build);
valid_T1_solve = ~isnan(T1_solve);

valid_T1= valid_T1_build;
valid_T2 = T2 > 0;

loglog(N(valid_T1), T1(valid_T1),'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'T_1 of Vanilla Solver');
hold on;
loglog(N(valid_T2), T2(valid_T2),'s-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'T_2 of Iterative Solver');

base_N = N(1);
base_T2 = T2(1);

C1 = 0.00000000004;
C2 = 0.0000005;

ref_N = N;
ref_ONlogN = C2 * ref_N .* log(ref_N);
ref_ON3    = C1 * ref_N.^3;


loglog(ref_N, ref_ONlogN, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference O(NlogN)');
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