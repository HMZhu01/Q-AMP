function [MSE_X_List, MSE_X1_List, Hmx, Res] = Real_DHSM(Input, Obj, denoiser)
	% x: N * 1, z: M * 1

	% Load parameters
	Height = Input.Height;
	Width = Input.Width;

	M = Input.M;
	N = Input.N;

	DDenoiser = @(noisy, sigma_hat) DHSM_Denoiser(noisy, sigma_hat, Width, Height, denoiser);

	nuw = Input.nuw;

	iter_num = Input.iter_num;

	mes = Input.mes;

	Out = Input.Out;

	MSE_X_List = zeros(iter_num, 1);
	MSE_X1_List = zeros(iter_num, 1);

	MSE_X_O = 1e3;
	MSE_X1_O = 1e3;

	EPS = 5e-13;

	outFileName = sprintf(mfilename);

	TX = 0.5;

	x = Obj.x;
	A_List = Obj.A_List;
	y = Obj.y;

	A2 = abs(A_List).^(2);
	TZ = TX^(2) * mean(A2(:)) * double(N)^(2);

	% Initialization
	% N * 1
	mX = Get_Initialisation(A_List, y);
	vX = TX * ones(N, 1);

	mZ = zeros(M, 1);
	vZ = TZ * ones(M, 1);

	% Initialize messages
	% M * 1
	v_PZ_inv = 1 ./ vZ;
	v_PZ_inv_O = v_PZ_inv;
	q_PZ = mZ .* v_PZ_inv;
	q_PZ_O = q_PZ;

	% N * 1
	v_PTX_inv = 1 ./ vX;
	v_PTX_inv_O = v_PTX_inv;
	q_PTX = mX .* v_PTX_inv;
	q_PTX_O = q_PTX;

	% N * 1 * M, N * N * M
	m_PX1_List = zeros(N, 1, M);
	C_PX1_List = zeros(N, N, M);
	m_PX1_List = mX;
	C_PX1_List = diag(vX);

	% N * 1 * M, N * N * M
	m_PX2_List = zeros(N, 1, M);
	C_PX2_List = zeros(N, N, M);
	m_PX2_List = mX;
	C_PX2_List = diag(vX);

	Hm_PTX_O = mX;

	m_STX_List = zeros(N, iter_num);
	v_STX_List = zeros(N, iter_num);
	Residual_List = zeros(N, iter_num);
	Hm_PTX_List = zeros(N, iter_num);
	Hv_PTX_List = zeros(N, iter_num);

	% Iteration
	for iter = 1 : iter_num

		% Backward passing
		% The messages of z
		% M * 1
		[Hm_SZ, Hv_SZ] = Out.Estimation(y, nuw, q_PZ ./ v_PZ_inv, 1 ./ v_PZ_inv);
		Hv_SZ = real(Hv_SZ);
		Hv_SZ = max(Hv_SZ, EPS);

		v_SZ = Hv_SZ ./ (1 - Hv_SZ .* v_PZ_inv);
		v_SZ = max(v_SZ, eps);
		m_SZ = v_SZ .* (Hm_SZ ./ Hv_SZ - q_PZ);
		m2_SZ = abs(m_SZ).^(2);

		% Update a_PX1, a_PX2
		A_mPX1 = pagemtimes(A_List, 'none', m_PX1_List, 'none');
		A_CPX1_A = pagemtimes( ...
				pagemtimes(A_List, 'none', C_PX1_List, 'none'), ...
				'none', A_List, 'transpose' ...
			);

		Numerator = pagemtimes( ...
				pagemtimes(A_mPX1, 'transpose', A_CPX1_A, 'none'), ...
				'none', A_mPX1, 'none' ...
			);
		Denominator = (pagemtimes(A_mPX1, 'transpose', A_mPX1, 'none')).^(2);
		a_PX1 = real(Numerator ./ Denominator);

		A_mPX2 = pagemtimes(A_List, 'transpose', m_PX2_List, 'none');
		A_CPX2_A = pagemtimes( ...
				pagemtimes(A_List, 'transpose', C_PX2_List, 'none'), ...
				'none', A_List, 'none' ...
			);

		Numerator = pagemtimes( ...
				pagemtimes(A_mPX2, 'transpose', A_CPX2_A, 'none'), ...
				'none', A_mPX2, 'none' ...
			);
		Denominator = (pagemtimes(A_mPX2, 'transpose', A_mPX2, 'none')).^(2);
		a_PX2 = real(Numerator ./ Denominator);

		% Backward TX
		% The messages of TX
		v_SZ_List = zeros(1, 1, M);
		m_SZ_List = zeros(1, 1, M);
		m2_SZ_List = zeros(1, 1, M);

		v_SZ_List(1, :, :) = v_SZ';
		m_SZ_List(1, :, :) = m_SZ';
		m2_SZ_List(1, :, :) = m2_SZ';

		tmp1 = v_SZ_List + m2_SZ_List .* a_PX1;
		tmp2 = v_SZ_List + m2_SZ_List .* a_PX2;

		A_mPX2 = pagemtimes(A_List, 'transpose', m_PX2_List, 'none');
		A_mPX2_A_mPX2 = pagemtimes(A_mPX2, 'none', A_mPX2, 'transpose');
		L_STX1_List = A_mPX2_A_mPX2 ./ tmp2;
		mSZ_A_mPX2 = pagemtimes(m_SZ_List, 'none', A_mPX2, 'none');
		b_STX1_List = mSZ_A_mPX2 ./ tmp2;

		A_mPX1 = pagemtimes(A_List, 'none', m_PX1_List, 'none');
		A_mPX1_A_mPX1 = pagemtimes(A_mPX1, 'none', A_mPX1, 'transpose');
		L_STX2_List = A_mPX1_A_mPX1 ./ tmp1;
		mSZ_A_mPX1 = pagemtimes(m_SZ_List, 'none', A_mPX1, 'none');
		b_STX2_List = mSZ_A_mPX1 ./ tmp1;

		L_STX = sum(L_STX1_List, 3) + sum(L_STX2_List, 3);
		b_STX = sum(b_STX1_List, 3) + sum(b_STX2_List, 3);

		HC_STX = inv(L_STX + diag(v_PTX_inv));
		Hm_STX = HC_STX * (b_STX + q_PTX);
		Hv_STX = R_D(HC_STX);
		Hv_STX = max(Hv_STX, EPS);

		v_STX = Hv_STX ./ (1 - Hv_STX .* v_PTX_inv);
		v_STX = max(v_STX, eps);
		m_STX = v_STX .* (Hm_STX ./ Hv_STX - q_PTX);

		m_STX_List(:, iter) = m_STX;
		v_STX_List(:, iter) = v_STX;
		Residual_List(:, iter) = (m_STX - x) ./ sqrt(v_STX);

		% Forward passing
		% Forward TX
		% The message of TX
		% N * 1
		% [Hm_PTX, Hv_PTX] = In_X.Estimation(m_STX, v_STX);
		v_tx_s = sqrt(mean(v_STX));

		Hm_PTX = DDenoiser(m_STX, v_tx_s);
		
		epsilon = max( ...
				max(abs(m_STX)) / 1000, 0.00001 ...
			);
		eta = randn(1, N);
		div = eta * ( ...
				( ...
					DDenoiser(m_STX + epsilon .* eta', v_tx_s) - Hm_PTX ...
				) ./ epsilon ...
			) ./ double(N);

		Hv_PTX = div * v_STX;

		Hm_PTX_List(:, iter) = Hm_PTX;
		Hv_PTX_List(:, iter) = Hv_PTX;

		MSE_X = norm(Hm_PTX - x, 'fro').^(2) / norm(x, 'fro').^(2);
		MSE_X1 = Dist(Hm_PTX, x);

		% stop = 1 break; stop = 0 continues
		stop = Judge(Hm_PTX) + Judge(Hv_PTX);

		if stop
			disp([mfilename, ', num = ', num2str(iter), ', X, stop = ', num2str(stop)])
			MSE_X_List(max(1, iter - 1) : iter_num, 1) = MSE_X_O;
			MSE_X1_List(max(1, iter - 1) : iter_num, 1) = MSE_X1_O;

			Hmx = Hm_PTX_O;
			break;
		end
		MSE_X_List(iter, 1) = MSE_X;
		MSE_X_O = MSE_X;

		MSE_X1_List(iter, 1) = MSE_X1;
		MSE_X1_O = MSE_X1;

		Hmx = Hm_PTX;
		Hm_PTX_O = Hm_PTX;

		Hv_PTX = max(Hv_PTX, 0.00001);

		% Hv_PTX = real(Hv_PTX);
		% Hv_PTX = max(Hv_PTX, EPS);

		v_PTX_inv = (v_STX - Hv_PTX) ./ (Hv_PTX .* v_STX);
		q_PTX = (Hm_PTX .* v_STX - m_STX .* Hv_PTX) ./ (Hv_PTX .* v_STX);

		IdxX = v_PTX_inv < eps;
		v_PTX_inv(IdxX) = v_PTX_inv_O(IdxX);
		q_PTX(IdxX) = q_PTX_O(IdxX);

		% Damping
		[v_PTX_inv, v_PTX_inv_O] = Damping(v_PTX_inv, v_PTX_inv_O, mes);
		[q_PTX, q_PTX_O] = Damping(q_PTX, q_PTX_O, mes);

		% Forward X
		% The message of X
		HC_PX = inv(L_STX + diag(v_PTX_inv));
		Hm_PX = HC_PX * (b_STX + q_PTX);
		Hv_PX = R_D(HC_PX);
		Hv_PX = max(Hv_PX, EPS);
		HC_PX(eye(N, 'logical')) = Hv_PX;

		O_m_PX1_List = m_PX1_List;
		O_m_PX2_List = m_PX2_List;

		A_OmPX1 = pagemtimes(A_List, 'none', O_m_PX1_List, 'none');
		A_OmPX1_HCPX_A_OmPX1 = pagemtimes( ...
				pagemtimes(A_OmPX1, 'transpose', HC_PX, 'none') ...
				, 'none', A_OmPX1, 'none' ...
			);
		tmp1 = v_SZ_List + m2_SZ_List .* a_PX1 - A_OmPX1_HCPX_A_OmPX1;

		A_OmPX2 = pagemtimes(A_List, 'transpose', O_m_PX2_List, 'none');
		A_OmPX2_HCPX_A_OmPX2 = pagemtimes( ...
				pagemtimes(A_OmPX2, 'transpose', HC_PX, 'none') ...
				, 'none', A_OmPX2, 'none' ...
			);
		tmp2 = v_SZ_List + m2_SZ_List .* a_PX2 - A_OmPX2_HCPX_A_OmPX2;

		tmp2_0 = pagemtimes(HC_PX, 'none', A_OmPX2, 'none');
		tmp2_1 = pagemtimes(tmp2_0, 'none', tmp2_0, 'transpose');
		tmp2_2 = pagemtimes(A_OmPX2, 'transpose', Hm_PX, 'none');

		C_PX1_List = HC_PX + tmp2_1 ./ tmp2;
		m_PX1_List = Hm_PX - (m_SZ_List - tmp2_2) ./ tmp2 .* tmp2_0;

		tmp1_0 = pagemtimes(HC_PX, 'none', A_OmPX1, 'none');
		tmp1_1 = pagemtimes(tmp1_0, 'none', tmp1_0, 'transpose');
		tmp1_2 = pagemtimes(A_OmPX1, 'transpose', Hm_PX, 'none');

		C_PX2_List = HC_PX + tmp1_1 ./ tmp1;
		m_PX2_List = Hm_PX - (m_SZ_List - tmp1_2) ./ tmp1 .* tmp1_0;

		% The message of Z
		A_mPX2 = pagemtimes(A_List, 'transpose', m_PX2_List, 'none');
		A_CPX2_A = pagemtimes( ...
				pagemtimes(A_List, 'transpose', C_PX2_List, 'none') ...
				, 'none', A_List, 'none' ...
			);

		mPZ = pagemtimes(A_mPX2, 'transpose', m_PX1_List, 'none');
		mPZ = squeeze(mPZ);

		tmp = bsxfun(@times, eye(N), ...
				pagemtimes(A_CPX2_A, 'none', C_PX1_List, 'none') ...
			);
		tmp1 = sum(tmp, [1, 2]);

		vPZ = tmp1 + ...
			pagemtimes( ...
				pagemtimes(m_PX1_List, 'transpose', A_CPX2_A, 'none') ...
				, 'none', m_PX1_List, 'none' ...
			) + ...
			pagemtimes( ...
				pagemtimes(A_mPX2, 'transpose', C_PX1_List, 'none') ...
				, 'none', A_mPX2, 'none' ...
			);
		vPZ = real(squeeze(vPZ));

		v_PZ_inv = 1 ./ vPZ;
		q_PZ = mPZ .* v_PZ_inv;

		IdxZ = v_PZ_inv < 0;
		v_PZ_inv(IdxZ) = v_PZ_inv_O(IdxZ);
		q_PZ(IdxZ) = q_PZ_O(IdxZ);

		% Damping
		[v_PZ_inv, v_PZ_inv_O] = Damping(v_PZ_inv, v_PZ_inv_O, mes);
		[q_PZ, q_PZ_O] = Damping(q_PZ, q_PZ_O, mes);

		'----';
	end
	Res.m_STX_List = m_STX_List;
	Res.v_STX_List = v_STX_List;
	Res.Residual_List = Residual_List;
	Res.Hm_PTX_List = Hm_PTX_List;
	Res.Hv_PTX_List = Hv_PTX_List;
end