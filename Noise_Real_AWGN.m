function w = Noise_Real_AWGN(M, L, nuw)
	w = sqrt(nuw) * randn(M, L);
end