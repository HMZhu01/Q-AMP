classdef In_Real_BPSK_Estimation

	properties
		All_Symbol;
		TX;
	end

	methods
		% Constructor
		function Obj = In_Real_BPSK_Estimation()
			Obj.All_Symbol = [- 1, 1];
			Obj.TX = 1;
		end

		function X = Generation(Obj, N, K)
			Data = round(rand(N, K));
			X = 2 * Data - 1;
		end

		function [Hmx, Hvx] = Estimation(Obj, mx, vx, Init_Flag, Init_X)

			Flag = false;
			X = 0;

			switch nargin
				case 3
					Flag = false;
				case 4
					Flag = Init_Flag;
					X = zeros(size(mx));
				case 5
					Flag = Init_Flag;
					X = Init_X;
				otherwise
					throw(MException('Foo:FatalError', 'In_Real_BPSK_Estimation.Estimation(...)'));
			end

			if ~ Flag

				[N, K] = size(mx);

				vx = reshape(vx, N * K, 1);
				mx = reshape(mx, N * K, 1);

				tmp = [mx, - mx];
				m_tmp = max(tmp, [], 2);
				n_tmp = tmp - m_tmp;
				E_tmp = exp(n_tmp ./ vx);

				Hmx = (E_tmp(:, 1) - E_tmp(:, 2)) ./ (E_tmp(:, 1) + E_tmp(:, 2));
				Hvx = 1 - Hmx.^(2);

				Hvx = max(Hvx, eps);

				Hmx = reshape(Hmx, N, K);
				Hvx = reshape(Hvx, N, K);
			else
				Hmx = x;
				Hvx = zeros(size(x));
			end

		end

		function [Hmx, Hvx_0, Hvx_1] = Estimation_1RSB(Obj, mx, vx_0, vx_1, L, Init_Flag, Init_X)

			Flag = false;
			X = 0;

			switch nargin
				case 5
					Flag = false;
				case 6
					Flag = Init_Flag;
					X = zeros(size(mx));
				case 7
					Flag = Init_Flag;
					X = Init_X;
				otherwise
					throw(MException('Foo:FatalError', 'In_Real_BPSK_Estimation.Estimation_1RSB(...)'));
			end

			if ~ Flag

				[N, K] = size(mx);

				mx = reshape(mx, N * K, 1);
				vx_0 = reshape(vx_0, N * K, 1);
				vx_1 = reshape(vx_1, N * K, 1);

				tmp1_f = @(bx, m, v0, v1) Real_Gaussian(bx, m, v0) .* Obj.NpN_l(bx, v1, L);
				tmp1 = integral(@(bx) tmp1_f(bx, mx, vx_0, vx_1), - inf, inf, 'ArrayValued', true);

				tmp2_f = @(bx, m, v0, v1) Real_Gaussian(bx, m, v0) .* Obj.NpN_l(bx, v1, L) .* Obj.Est_x(bx, v1);
				tmp2 = integral(@(bx) tmp2_f(bx, mx, vx_0, vx_1), - inf, inf, 'ArrayValued', true);

				tmp3_f = @(bx, m, v0, v1) Real_Gaussian(bx, m, v0) .* Obj.NpN_l(bx, v1, L) .* Obj.Est_xx(bx, v1);
				tmp3 = integral(@(bx) tmp3_f(bx, mx, vx_0, vx_1), - inf, inf, 'ArrayValued', true);

				Hmx = tmp2 ./ tmp1;

				Hvx_0 = tmp3 ./ tmp1 - Hmx.^(2);
				Hvx_1 = 1 - tmp3 ./ tmp1;

				Hvx_0 = max(Hvx_0, eps);
				Hvx_1 = max(Hvx_1, eps);

				Hmx = reshape(Hmx, N, K);
				Hvx_0 = reshape(Hvx_0, N, K);
				Hvx_1 = reshape(Hvx_1, N, K);

			else
				Hmx = X;
				Hvx_0 = zeros(size(x));
				Hvx_1 = zeros(size(x));
			end

		end

		function ret = NpN_l(Obj, bx, v1, l)
			tmp = Real_Gaussian(bx, 1 * ones(size(v1)), v1) + Real_Gaussian(bx, - 1 * ones(size(v1)), v1);
			ret = tmp.^(l);
			stop = Judge(ret);
			if stop
				a = isinf(ret);
				sum(a);
				'----------';
			end
			ret = real(ret);
		end

		%% bx is the integral variable, 1 * 1
		%% v1 is the variances, L * 1
		function ret = Est_x(Obj, bx, v1)

			tmp = [bx, - bx];
			m_tmp = max(tmp, [], 2);
			n_tmp = tmp - m_tmp;
			E_tmp = exp(n_tmp ./ v1);

			ret = (E_tmp(:, 1) - E_tmp(:, 2)) ./ (E_tmp(:, 1) + E_tmp(:, 2));

			stop = Judge(ret);
			if stop
				a = isinf(ret);
				sum(a);
				'----------';
			end
			ret = real(ret);
		end

		%% bx is the integral variable, 1 * 1
		%% v1 is the variances, L * 1
		function ret = Est_xx(Obj, bx, v1)

			tmp = [bx, - bx];
			m_tmp = max(tmp, [], 2);
			n_tmp = tmp - m_tmp;
			E_tmp = exp(n_tmp ./ v1);

			ret1 = (E_tmp(:, 1) - E_tmp(:, 2)) ./ (E_tmp(:, 1) + E_tmp(:, 2));
			ret = ret1.^(2);

			stop = Judge(ret);
			if stop
				a = isinf(ret);
				sum(a);
				'----------';
			end
			ret = real(ret);
		end

	end
end