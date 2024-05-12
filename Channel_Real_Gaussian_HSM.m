function A_List = Channel_Real_Gaussian_HSM(M, N)
	A_List = randn(N, N, M) / sqrt(double(N * N));
end