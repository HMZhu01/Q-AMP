function obj = MIMO_System_Real_Gaussian_HSM(Input)
	%% load parameters
	M = Input.M;
	N = Input.N;

	nuw = Input.nuw;
	
	In_X = Input.In_X;
	Out = Input.Out;

	%% Generate x
	x = In_X.Generation([N, 1]);

	A_List = Channel_Real_Gaussian_HSM(M, N);

	%% Noise
	w = Noise_Real_AWGN(M, 1, nuw);

	%% Uncoded system
	[z, z_n] = Process_HSM(x, A_List, w);

	y = Out.Quantization(z_n);

	%% load Inputmeters
	obj.x = x;
	obj.A_List = A_List;
	obj.w = w;
	obj.z = z;
	obj.z_n = z_n;
	obj.y = y;
end