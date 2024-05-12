function Main_DHSM_Iter(index, nuw)

	addpath(genpath('./BM3D'));
	addpath(genpath('./BM3D/BM3D-SAPCA'));

	addpath('../');
	addpath('../Operator');
	addpath('../D-AMP_Toolbox-master/Utils');
	addpath('../matconvnet/matlab');
	addpath('./TestingData');

	% Parameters
	denoiser_BM3D = 'BM3D';
	denoiser_DnCNN = 'DnCNN';

	% Other option is 17
	n_DnCNN_layers = 20;

	Dir = fullfile('../', 'D-AMP_Toolbox-master');
	LoadNetworkWeights(Dir, n_DnCNN_layers);

	File_List = {
			'Akun_16.jpg', ...
			'Alek_16.jpg', ...
			'Cat_16.jpg', ...
			'Cloud_16.jpg', ...
			'Drive_16.jpg', ...
			'Hexagon_16.jpg', ...
			'Moon_16.jpg', ...
			'Pin_16.jpg', ...
			'Shield_16.jpg', ...
			'Star_16.jpg', ...
			'Triangle_16.jpg' ...
		};

	index = min(index, numel(File_List));

	delta = 4; 

	Height = 16;
	Width = 16;
	N = int16(Height * Width);
	M = int16(delta * N);

	% 1 on; 0 off.
	ADC_switch = 0;
	bit = 5;

	% nuw = 1 / 255;

	iter_num = 30;
% 	max_iterations = 2500;
	max_iterations = 1e5;

	mes = 0.95;

	outFileName = sprintf(mfilename);
	Str = [outFileName, '-', ...
			denoiser_BM3D, '-', ...
			denoiser_DnCNN, '-', ...
			'Height', num2str(Height), ...
			'Width', num2str(Width), ...
			'M', num2str(M), ...
			'N', num2str(N), ...
			'A', num2str(ADC_switch), '-', ...
			'B', num2str(bit), '-', ...
			'n', num2str(nuw), '-', ...
			'I', num2str(iter_num), '-', ...
			'm', num2str(mes), '-', ...
			'i', num2str(index) ...
		];

	% load parameters

	Input.File_List = File_List;
	Input.index = index;

	Input.Height = Height;
	Input.Width = Width;

	Input.M = M;
	Input.N = N;

	Input.ADC_switch = ADC_switch;
	Input.bit = bit;

	Input.nuw = nuw;

	Input.iter_num = iter_num;
	Input.max_iterations = max_iterations;

	Input.mes = mes;

	Input.Out = Out_Real_Quantization_Estimation(ADC_switch, bit);

	BM3D_HSM_X = zeros(iter_num, 1);
	BM3D_HSM_X1 = zeros(iter_num, 1);
	BM3D_HSM_Hmx = zeros(N, 1);

	DnCNN_HSM_X = zeros(iter_num, 1);
	DnCNN_HSM_X1 = zeros(iter_num, 1);
	DnCNN_HSM_Hmx = zeros(N, 1);

	DSolve_X = zeros(max_iterations, 1);
	DSolve_X1 = zeros(max_iterations, 1);

	DSolve_iters = zeros(1, 1);
	DSolve_Hmx = zeros(N, 1);

	obj = MIMO_System_RG_DHSM(Input);

	[BM3D_HSM_X(:, 1), BM3D_HSM_X1(:, 1), BM3D_HSM_Hmx(:, 1), BM3D_HSM_Res] = Real_DHSM(Input, obj, denoiser_BM3D);
	[DnCNN_HSM_X(:, 1), DnCNN_HSM_X1(:, 1), DnCNN_HSM_Hmx(:, 1), DnCNN_HSM_Res] = Real_DHSM(Input, obj, denoiser_DnCNN);
	[DSolve_X(:, 1), DSolve_X1(:, 1), DSolve_iters(:, 1), DSolve_Hmx(:, 1)] = DSolve(Input, obj);
	'----------';

	BM3D_HSM_X_MI = zeros(max_iterations, 1);
	BM3D_HSM_X1_MI = zeros(max_iterations, 1);

	DnCNN_HSM_X_MI = zeros(max_iterations, 1);
	DnCNN_HSM_X1_MI = zeros(max_iterations, 1);

	BM3D_HSM_X_MI(1 : iter_num, 1) = BM3D_HSM_X;
	BM3D_HSM_X_MI(iter_num + 1 : end, 1) = BM3D_HSM_X(end, 1);

	BM3D_HSM_X1_MI(1 : iter_num, 1) = BM3D_HSM_X1;
	BM3D_HSM_X1_MI(iter_num + 1 : end, 1) = BM3D_HSM_X1(end, 1);

	DnCNN_HSM_X_MI(1 : iter_num, 1) = DnCNN_HSM_X;
	DnCNN_HSM_X_MI(iter_num + 1 : end, 1) = DnCNN_HSM_X(end, 1);

	DnCNN_HSM_X1_MI(1 : iter_num, 1) = DnCNN_HSM_X1;
	DnCNN_HSM_X1_MI(iter_num + 1 : end, 1) = DnCNN_HSM_X1(end, 1);

	iter = 1 : max_iterations;

	x = obj.x;

	figure;
	plot(iter, 10 * log10(BM3D_HSM_X_MI), '-+r');
	hold on;
	plot(iter, 10 * log10(BM3D_HSM_X1_MI), '-*r');
	hold on;
	plot(iter, 10 * log10(DnCNN_HSM_X_MI), '-+g');
	hold on;
	plot(iter, 10 * log10(DnCNN_HSM_X1_MI), '-*g');
	hold on;
	plot(iter, 10 * log10(DSolve_X), '-*b');
	hold on;
	plot(iter, 10 * log10(DSolve_X1), '-*b');
	hold on;

	legend( ...	
			'BM3D-HSM-X', ...
			'BM3D-HSM-X1', ...
			'DnCNN-HSM-X', ...
			'DnCNN-HSM-X1', ...
			'Solve', ...
			'Solve-1' ...
		);
	hold on;
	xlabel('Iter');
	hold on;
	ylabel('MSE-X');
	hold on;
	saveas(gcf, [Str, '-Figure1.fig']);
	saveas(gcf, [Str, '-Figure1.png']);
	close;

	figure;
	subplot(2, 2, 1); imshow(reshape(x, [Width, Height])); title('Original');
	subplot(2, 2, 2); imshow(reshape(BM3D_HSM_Hmx, [Width, Height])); title('BM3D-Estimated');
	subplot(2, 2, 3); imshow(reshape(DnCNN_HSM_Hmx, [Width, Height])); title('DnCNN-Estimated');
	subplot(2, 2, 4); imshow(reshape(DSolve_Hmx, [Width, Height])); title('DSolve-Estimated');
	saveas(gcf, [Str, '-Figure2.fig']);
	saveas(gcf, [Str, '-Figure2.png']);
	close;

	BM3D_HSM_PSNR = PSNR(x, BM3D_HSM_Hmx)
	DnCNN_HSM_PSNR = PSNR(x, DnCNN_HSM_Hmx)
	DSolve_PSNR = PSNR(x, DSolve_Hmx)

	disp(['BM3D_HSM_PSNR: ', num2str(BM3D_HSM_PSNR)]);
	disp(['DnCNN_HSM_PSNR: ', num2str(DnCNN_HSM_PSNR)]);
	disp(['DSolve_PSNR: ', num2str(DSolve_PSNR)]);
	disp(['DSolve_iters: ', num2str(DSolve_iters)])

	save([Str, '.mat']);

	"---";
end