clc;
clear;

addpath('./Operator');

% Parameters
delta = 4;
N = 256;
M = round(delta * N);

q = 0.55;

% 1 on; 0 off.
ADC_switch = 0;
bit = 5;

% SNR = 10;
nuw = 0.1;
% nuw = 10 / 255;

% Iteration number
iter_num = 30;
% max_iterations = 2500;
max_iterations = 1e5;

mes = 0.95;

mesB = 0.90;
mesS = 0.90;

% Test number
test_num = 2e1;

outFileName = sprintf(mfilename);
Str = [outFileName, '-', ...
		'M', num2str(M), ...
		'N', num2str(N), ...
		'q', num2str(q), ...
		'A', num2str(ADC_switch), '-', ...
		'B', num2str(bit), '-', ...
		'n', num2str(nuw), '-', ...
		'I', num2str(iter_num), '-', ...
		'm', num2str(mes) ...
	];

% load parameters
Input.M = M;
Input.N = N;

Input.q = q;

Input.ADC_switch = ADC_switch;
Input.bit = bit;

% Input.nuw = 10^(-SNR / 10);
Input.nuw = nuw;

Input.iter_num = iter_num;
Input.max_iterations = max_iterations;

Input.mes = mes;

Input.mesB = mesB;
Input.mesS = mesS;

Input.In_X = In_Real_0_1_Estimation(q);
Input.Out = Out_Real_Quantization_Estimation(ADC_switch, bit);

HSM_X = zeros(iter_num, test_num);
HSM_X_iter = zeros(iter_num, 1);

HSM_X1 = zeros(iter_num, test_num);
HSM_X1_iter = zeros(iter_num, 1);

HSM_SE_X = zeros(iter_num, test_num);
HSM_SE_X_iter = zeros(iter_num, 1);

Solve_X = zeros(max_iterations, test_num);
Solve_X_iter = zeros(max_iterations, 1);

Solve_X1 = zeros(max_iterations, test_num);
Solve_X1_iter = zeros(max_iterations, 1);

Solve_iterations = zeros(1, test_num);
Solve_xt = zeros(N, test_num);

TWF_X = zeros(max_iterations, test_num);
TWF_X_iter = zeros(max_iterations, 1);

TWF_X1 = zeros(max_iterations, test_num);
TWF_X1_iter = zeros(max_iterations, 1);

TWF_iterations = zeros(1, test_num);
TWF_xt = zeros(N, test_num);

% Parfor_Progress(test_num);
% parfor ii = 1 : test_num
for ii = 1 : test_num
	obj = MIMO_System_Real_Gaussian_HSM(Input);

	[HSM_X(:, ii), HSM_X1(:, ii)] = Real_HSM(Input, obj);
	HSM_SE_X(:, ii) = Real_HSM_SE(Input);
	[Solve_X(:, ii), Solve_X1(:, ii), Solve_iterations(:, ii), Solve_xt(:, ii)] = Solve(Input, obj);
	[TWF_X(:, ii), TWF_X1(:, ii), TWF_iterations(:, ii), TWF_xt(:, ii)] = TWF_Solve(Input, obj);
	'----';	
% 	Parfor_Progress;
end
Parfor_Progress(0);

HSM_X_iter = mean(HSM_X, 2);
HSM_X1_iter = mean(HSM_X1, 2);

HSM_SE_X_iter = mean(HSM_SE_X, 2);

Solve_X_iter = mean(Solve_X, 2);
Solve_X1_iter = mean(Solve_X1, 2);

TWF_X_iter = mean(TWF_X, 2);
TWF_X1_iter = mean(TWF_X1, 2);

HSM_X_MI = zeros(max_iterations, 1);
HSM_X1_MI = zeros(max_iterations, 1);
HSM_SE_X_MI = zeros(max_iterations, 1);

HSM_X_MI(1 : iter_num, 1) = HSM_X_iter;
HSM_X_MI(iter_num + 1 : end, 1) = HSM_X_iter(end, 1);

HSM_X1_MI(1 : iter_num, 1) = HSM_X1_iter;
HSM_X1_MI(iter_num + 1 : end, 1) = HSM_X1_iter(end, 1);

HSM_SE_X_MI(1 : iter_num, 1) = HSM_SE_X_iter;
HSM_SE_X_MI(iter_num + 1 : end, 1) = HSM_SE_X_iter(end, 1);

iter = 1 : max_iterations;

figure;
plot(iter, 10 * log10(HSM_X_MI), '-+b');
hold on;
plot(iter, 10 * log10(HSM_X1_MI), '-*b');
hold on;
plot(iter, 10 * log10(HSM_SE_X_MI), '-sr');
hold on;

plot(iter, 10 * log10(Solve_X_iter), '-+g');
hold on;
plot(iter, 10 * log10(Solve_X1_iter), '-*g');
hold on;

plot(iter, 10 * log10(TWF_X_iter), '-+k');
hold on;
plot(iter, 10 * log10(TWF_X1_iter), '-*k');
hold on;
legend( ...
		'HSM-X', ...
		'HSM-X1', ...
		'HSM-SE-X', ...
		'Solve', ...
		'Solve-1', ...
		'TWF', ...
		'TWF-1' ...
	);
hold on;
xlabel('Iter');
hold on;
ylabel('MSE-X');
hold on;
saveas(gcf, [Str, '-Figure.fig']);
saveas(gcf, [Str, '-Figure.png']);
close;

save([Str, '.mat']);