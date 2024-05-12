classdef In_Real_Gaussian_Estimation
	
	properties
		TX;
	end
	
	methods
		% Constructor
		function obj = In_Real_Gaussian_Estimation(TX)
			obj.TX = TX;
		end

		function X = Generation(obj, N, L)
			TX = obj.TX;
			X = randn(N, L) * sqrt(TX);
		end

		function [HmX, HvX] = Estimation(obj, mX, vX)
			TX = obj.TX;
			HvX = 1 ./ (1 / TX + 1 ./ vX);
			HmX = HvX .* (mX ./ vX);
			HvX = real(HvX);
		end
	end
end

