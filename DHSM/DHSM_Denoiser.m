function denoised = DHSM_Denoiser(noisy, sigma_hat, width, height, denoiser)
	noisy = reshape(noisy, [width, height]);
	switch denoiser
		case 'BM3D'
			% noisy = noisy * 255;
			[NA, output] = BM3D_1(1, noisy, sigma_hat, 'np', 0);
			% output = 255 * output;
		case 'fast-BM3D'
			noisy = real(noisy);
			[NA, output] = BM3D_1(1, noisy, sigma_hat, 'lc', 0);
			% output =255 * output;
		case 'BM3D-SAPCA'
			% output = 255 * BM3DSAPCA2009((1 / 255) * noisy, (1 / 255) * sigma_hat);
			output = BM3DSAPCA2009(noisy, sigma_hat);
		case 'DnCNN'
			noisy = real(noisy);
			global_vars = who('global');
			if ~any(ismember(global_vars, 'net_300to500'));
				error('You need to run LoadNetworkWeights before you can use the DnCNN denoiser');
			end
			input = gpuArray(single(noisy));
			if sigma_hat > 500
				global net_500to1000
				res = vl_simplenn(net_500to1000, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 300
				global net_300to500
				res = vl_simplenn(net_300to500, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 150
				global net_150to300
				res = vl_simplenn(net_150to300, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 100
				global net_100to150
				res = vl_simplenn(net_100to150, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 80
				global net_80to100
				res = vl_simplenn(net_80to100, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 60
				global net_60to80
				res = vl_simplenn(net_60to80, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 40
				global net_40to60
				res = vl_simplenn(net_40to60, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 20
				global net_20to40
				res = vl_simplenn(net_20to40, input, [], [], 'conserveMemory', true, 'mode', 'test');
			elseif sigma_hat > 10
				global net_10to20
				res = vl_simplenn(net_10to20, input, [], [], 'conserveMemory', true, 'mode', 'test');
			else
				global net_0to10
				res = vl_simplenn(net_0to10, input, [], [], 'conserveMemory', true, 'mode', 'test');
			end
			output = input - res(end).x;
			output = double(gather(output));

		otherwise
			error('Unrecognized Denoiser')
	end
	denoised = output(:);
end