classdef Out_Real_AWGN_Estimation 
	
	properties
	end
	
	methods
		function obj = Out_Real_AWGN_Estimation()
		end

		function [HmZ, HvZ] = Estimation(obj, y, nuw, mZ, vZ)

			HvZ = 1 ./ (1 ./ vZ + 1 / nuw);
			HmZ = HvZ .* (mZ ./ vZ + y / nuw);

			% HvZ = real(HvZ);
		end
	
	end
end