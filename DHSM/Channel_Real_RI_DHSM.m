function A_List = Channel_Real_RI_DHSM(M, N)

	coef = randn(N * N, M); 
	coefs = coef / N;

	ensemble = reshape(coefs, [N, N, M]);
	ensemble = permute(ensemble, [2, 1, 3]);

	coefs = sum(abs(coefs).^(2), 2);
	kappa_vals = (1 ./ (coefs / M)).^(0.5);

	A_List = ensemble * mean(kappa_vals)

end