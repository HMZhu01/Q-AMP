clc;
clear;

P_Index = 7;

P_Output = ['n0.039216-I30-m0.95-i', ...
		num2str(P_Index)];

P_Name = ['Main_DHSM_Iter-BM3D-DnCNN-', ...
		'Height16Width16M1024N256A0-B5-', ...
		P_Output, '.mat' ...
	];

load(P_Name);

P_outFileName = sprintf(mfilename);
P_Str = [P_outFileName, '-', P_Output];

for i = 1 : iter_num

	BM3D_m_STX = BM3D_HSM_Res.m_STX_List(:, i);
	BM3D_Residual = BM3D_HSM_Res.Residual_List(:, i);
	BM3D_Hm_PTX = BM3D_HSM_Res.Hm_PTX_List(:, i);

	DnCNN_m_STX = DnCNN_HSM_Res.m_STX_List(:, i);
	DnCNN_Residual = DnCNN_HSM_Res.Residual_List(:, i);
	DnCNN_Hm_PTX = DnCNN_HSM_Res.Hm_PTX_List(:, i);

	figure;
	subplot(2, 4, 1);
	imshow(reshape(x, [Width, Height]));
	title('x');
	hold on;

	subplot(2, 4, 2);
	imshow(reshape(BM3D_m_STX, [Width, Height]));
	title('BM3D-m-STX');
	hold on;

	subplot(2, 4, 3);
	imshow(reshape(BM3D_Residual, [Width, Height]));
	title('BM3D-Residual');
	hold on;

	subplot(2, 4, 4);
	imshow(reshape(BM3D_Hm_PTX, [Width, Height]));
	title('BM3D-Hm-PTX');
	hold on;

	subplot(2, 4, 6);
	imshow(reshape(DnCNN_m_STX, [Width, Height]));
	title('DnCNN-m-STX');
	hold on;

	subplot(2, 4, 7);
	imshow(reshape(DnCNN_Residual, [Width, Height]));
	title('DnCNN-Residual');
	hold on;

	subplot(2, 4, 8);
	imshow(reshape(DnCNN_Hm_PTX, [Width, Height]));
	title('DnCNN-Hm-PTX');
	hold on;

	saveas(gcf, [P_Str, '-Iter', num2str(i), '.fig']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), '.png']);

	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - x.png']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - BM3D_m_STX.png']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - BM3D_Residual.png']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - BM3D_Hm_PTX.png']);

	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - DnCNN_m_STX.png']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - DnCNN_Residual.png']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - DnCNN_Hm_PTX.png']);

	close;

	figure;
	qqplot(BM3D_Residual)
	saveas(gcf, [P_Str, '-Iter', num2str(i), '-BM3D.fig']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), '-BM3D.png']);

	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - BM3D_QQplot.png']);
	close;

	figure;
	qqplot(DnCNN_Residual)
	saveas(gcf, [P_Str, '-Iter', num2str(i), '-DnCNN.fig']);
	saveas(gcf, [P_Str, '-Iter', num2str(i), '-DnCNN.png']);

	saveas(gcf, [P_Str, '-Iter', num2str(i), ' - DnCNN_QQplot.png']);
	close;
end