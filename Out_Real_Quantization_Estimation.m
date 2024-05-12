classdef Out_Real_Quantization_Estimation 

	properties
		ADC_switch;
		bit;
		quan_step;
		DeltaTh;
		Q_Out;
		Y_Sets;
		MAX_TERM;
	end

	methods
		function Obj = Out_Real_Quantization_Estimation(ADC_switch, bit)

			Obj.MAX_TERM = 1e5;
			Obj.ADC_switch = ADC_switch;
			Obj.bit = bit;

			if bit == 1
				quan_step = 1;
			else
				quan_step = 1 / (2^(bit - 1));    
			end

			DeltaTh = (0 : 1 : 2^(bit - 1) - 1) * quan_step;
			Q_Out = (1 : 2 : 2^(bit) - 1) * quan_step / 2;
			Y_Sets = (- 2^(bit) + 1 : 2 : 2^(bit) - 1) * quan_step / 2;

			Obj.quan_step = quan_step;
			Obj.DeltaTh = DeltaTh;
			Obj.Q_Out = Q_Out;
			Obj.Y_Sets = Y_Sets;
		end

		function Y = Quantization(Obj, Z)
			ADC_switch = Obj.ADC_switch;
			if 0 == ADC_switch
				Y = Z;
			elseif 1 == ADC_switch
				[M, K] = size(Z);
				Z = reshape(Z, M * K, 1);

				quan_step = Obj.quan_step;
				DeltaTh = Obj.DeltaTh;
				Q_Out = Obj.Q_Out;

				Z_R = Q_Out(end) * (abs(Z) >= DeltaTh(end));
				for bIdx = length(DeltaTh) : - 1 : 2
					Z_R = Z_R + Q_Out(bIdx - 1) * ((abs(Z) < DeltaTh(bIdx)) & (abs(Z) >= DeltaTh(bIdx - 1)));
				end
				Z_R = sign(Z) .* Z_R;
				Y = reshape(Z_R, M, K);
			else
				throw(MException('Foo:FatalError', 'Out_Real_Quantization_Estimation.Quantization(...)'));
			end

		end

		function [hatz, hatv] = Estimation(Obj, init_Y, nuw, init_Z, init_V)

			ADC_switch = Obj.ADC_switch;
			MAX_TERM = Obj.MAX_TERM;

			if 0 == ADC_switch
				hatv = 1 ./ (1 ./ init_V + 1 / nuw);
				hatz = hatv .* (init_Z ./ init_V + init_Y / nuw);
				
				hatv = max(hatv, eps);
				%hatv = (init_V .* nuw) ./ (nuw + init_V);
				%hatz = (init_Z .* nuw + init_Y .* init_V) ./ (nuw + init_V);

			elseif 1 == ADC_switch

				quan_step = Obj.quan_step;
				DeltaTh = Obj.DeltaTh;
				[M, L] = size(init_Y);

				Y = reshape(init_Y, M * L, 1);
				Z = reshape(init_Z, M * L, 1);
				V = reshape(init_V, M * L, 1);

				y_up = Y + quan_step / 2;
				y_low = Y - quan_step / 2;
				[pos1, ~] = find(Y > max(DeltaTh));
				[pos2, ~] = find(Y < - max(DeltaTh));
				y_up(pos1) = MAX_TERM;
				y_low(pos2) = - MAX_TERM;

				eta1 = (y_up - Z) ./ sqrt(nuw + V);
				eta2 = (y_low - Z) ./ sqrt(nuw + V);

				tem1 = normpdf(eta1) - normpdf(eta2);
				tem2 = normcdf(eta1) - normcdf(eta2);
				tem2 = tem2 + eps;
				tem3 = eta1 .* normpdf(eta1) - eta2 .* normpdf(eta2);

				z_tem = Z - V ./ sqrt(nuw + V) .* (tem1 ./ tem2);
				v_tem = V - V.^(2) ./ (nuw + V) .* (tem3 ./ tem2 + (tem1 ./ tem2).^(2));

				judge = Judge(z_tem) + Judge(v_tem);
				if judge
					[hatz, hatv] = Obj.Estimation_Damping(init_Y, nuw, init_Z, init_V);
				else
					hatz = z_tem;
					hatv = max(v_tem, eps);

					hatz = reshape(hatz, M, L);
					hatv = reshape(hatv, M, L);
				end
			else
				throw(MException('Foo:FatalError', 'Out_Real_Quantization_Estimation.Estimation(...)'));
			end
		end

		function [hatz, hatv] = Estimation_Damping(Obj, init_Y, nuw, init_Z, init_V)

			ADC_switch = Obj.ADC_switch;
			MAX_TERM = Obj.MAX_TERM;

			if 0 == ADC_switch
				hatv = 1 ./ (1 ./ init_V + 1 / nuw);
				hatz = hatv .* (init_Z ./ init_V + init_Y / nuw);
				%hatv = (init_V .* nuw) ./ (nuw + init_V);
				%hatz = (init_Z .* nuw + init_Y .* init_V) ./ (nuw + init_V);

			elseif 1 == ADC_switch

				quan_step = Obj.quan_step;
				DeltaTh = Obj.DeltaTh;
				[M, L] = size(init_Y);
				Y = reshape(init_Y, M * L, 1);
				Z = reshape(init_Z, M * L, 1);
				V = reshape(init_V, M * L, 1);

				y_up = Y + quan_step / 2;
				y_low = Y - quan_step / 2;
				[pos1, ~] = find(Y > max(DeltaTh));
				[pos2, ~] = find(Y < -max(DeltaTh));
				y_up(pos1) = MAX_TERM;
				y_low(pos2) = -MAX_TERM;

				eta1 = (y_up - Z) ./ sqrt(nuw + V);
				eta2 = (y_low - Z) ./ sqrt(nuw + V);

				tem1 = normpdf(eta1) - normpdf(eta2);
				tem2 = normcdf(eta1) - normcdf(eta2);
				tem2 = tem2 + eps;
				tem3 = eta1 .* normpdf(eta1) - eta2 .* normpdf(eta2);
				
				pos = eta2 < - 100; 
				tem1(pos) = normpdf(eta1(pos));
				tem2(pos) = normcdf(eta1(pos));
				tem2 = tem2 + eps;
				tem3(pos) = eta1(pos) .* normpdf(eta1(pos));

				z_tem = Z - V ./ sqrt(nuw + V) .* (tem1 ./ tem2);
				v_tem = V - V.^(2) ./ (nuw + V) .* (tem3 ./ tem2 + (tem1 ./ tem2).^(2));

				hatz = z_tem;
				hatv = max(v_tem, eps);

				hatz = reshape(hatz, M, L);
				hatv = reshape(hatv, M, L); 
			else
				throw(MException('Foo:FatalError', 'Out_Real_Quantization_Estimation.Estimation_Damping(...)'));
			end
		end

		function vz_sub = MSE_SE(obj, Tz, v1_plus, nuw)
			ADC_switch = obj.ADC_switch;
			MAX_TERM = obj.MAX_TERM;
			
			if 0 == ADC_switch
				vz_sub = (nuw * v1_plus) / (v1_plus + nuw);
			elseif 1 == ADC_switch
				quan_step = obj.quan_step;
				Y_Sets = obj.Y_Sets;
				Y_Up = Y_Sets + quan_step / 2;
				Y_Up(end) = Inf;
				Y_Low = Y_Sets - quan_step / 2;
				Y_Low(1) = -Inf;
				alpha = zeros(length(Y_Sets), 1);
				Gaussian = @(x, m, v) 1 ./ sqrt(2 * pi * v) .* exp(-(x - m).^(2) ./ (2 * v));
				for index = 1 : length(Y_Sets)	  
					alpha(index) = integral(@(u) ...
						( ...
							normpdf( ...
								( ...
									Y_Up(index) - sqrt(Tz - v1_plus) * u ...
								) / sqrt(v1_plus + nuw) ...
							) - normpdf( ...
								( ...
									Y_Low(index) - sqrt(Tz - v1_plus) * u ...
								) / sqrt(v1_plus + nuw) ...
							) ...
						).^(2) ./ ( ...
							normcdf( ...
								( ...
									Y_Up(index) - sqrt(Tz - v1_plus) * u ...
								) / sqrt(v1_plus + nuw) ...
							) - normcdf( ...
								( ...
									Y_Low(index) - sqrt(Tz - v1_plus) * u ...
								) / sqrt(v1_plus + nuw) ...
							) + eps ...
						) .* Gaussian(u, 0, 1), - Inf, Inf ...
					);
				end
				vz_sub = v1_plus - (v1_plus)^(2) / (v1_plus + nuw) * sum(alpha);
			end
		end
	end
end