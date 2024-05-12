function MSE_X_List = Real_HSM_SE(Input)
	% x: N * 1, z: M * 1
	% Load parameters
	M = Input.M;
	N = Input.N;

	nuw = Input.nuw;

	iter_num = Input.iter_num;

	mesB = Input.mesB;
	mesS = Input.mesS;

	In_X = Input.In_X;

	Out = Input.Out;

	MSE_X_List = zeros(iter_num, 1);

	EPS = eps;
	% EPS = 5e-13;

	outFileName = sprintf(mfilename);

	NUM = 5e0;

	% Generation
	Obj = MIMO_System_Real_Gaussian_HSM(Input);

	A_List = Obj.A_List;
	x = Obj.x;
	y = Obj.y;

	TX = In_X.TX;
	TZ = N^(2) * var(A_List(:)) * TX^(2);

	% Initialization
	% N * 1
	mX = Get_Initialisation(A_List, y);
	vX = TX;

	vZ = TZ - 1e-6;

	% Initialize messages
	v_p_z = vZ;
	v_p_z_O = v_p_z;

	v_p_tx = vX;
	v_p_tx_O = v_p_tx;
	m_p_tx = repmat(mX, [1, 1, 1, NUM]);
	m_p_tx_O = m_p_tx;

	c_p1_x = vX;
	c_p1_x_O = c_p1_x;
	m_p1_x_List = repmat(mX, [1, 1, M, NUM]);
	m_p1_x_List_O = m_p1_x_List;

	c_p2_x = vX;
	c_p2_x_O = c_p2_x;
	m_p2_x_List = repmat(mX, [1, 1, M, NUM]);
	m_p2_x_List_O = m_p2_x_List;

	m2_p1_x_List = pagemtimes(m_p1_x_List, 'none', m_p1_x_List, 'transpose');
	m2_p2_x_List = pagemtimes(m_p2_x_List, 'none', m_p2_x_List, 'transpose');

	M2_p1_x_List = squeeze(sum(m2_p1_x_List, 3));
	M2_p2_x_List = squeeze(sum(m2_p2_x_List, 3));

	m_p1_a_List = pagemtimes(A_List, 'none', m_p1_x_List, 'none');
	m_p2_a_List = pagemtimes(A_List, 'transpose', m_p2_x_List, 'none');

	m2_p1_a_List = pagemtimes(m_p1_a_List, 'none', m_p1_a_List, 'transpose');
	m2_p2_a_List = pagemtimes(m_p2_a_List, 'none', m_p2_a_List, 'transpose');

	M2_p1_a_List = squeeze(sum(m2_p1_a_List, 3));
	M2_p2_a_List = squeeze(sum(m2_p2_a_List, 3));

	m_p1_a_S_List = m_p1_a_List(:, :, 2 : end, :);
	m_p2_a_S_List = m_p2_a_List(:, :, 2 : end, :);

	m2_p1_a_S_List = pagemtimes(m_p1_a_S_List, 'none', m_p1_a_S_List, 'transpose');
	m2_p2_a_S_List = pagemtimes(m_p2_a_S_List, 'none', m_p2_a_S_List, 'transpose');

	M2_p1_a_S_List = squeeze(sum(m2_p1_a_S_List, 3));
	M2_p2_a_S_List = squeeze(sum(m2_p2_a_S_List, 3));

	A2_List = pagemtimes(A_List, 'none', A_List, 'transpose');
	A2_Sum = squeeze(sum(A2_List, 3));

	% Iteration
	for iter = 1 : iter_num

		% Backward passing
		% The messages of Z
		Hv_s_z = Out.MSE_SE(TZ, v_p_z, nuw);
		Hv_s_z = real(Hv_s_z);
		Hv_s_z = max(Hv_s_z, EPS);

		% v_s_z = 1 / (1 / Hv_s_z - 1 / v_p_z);
		v_s_z = (Hv_s_z * v_p_z) / (v_p_z - Hv_s_z);
		v_s_z = max(v_s_z, EPS);

		w_s_z = TZ + v_s_z;

		% Update a_s1_x and a_s2_x
		a_s1_x = c_p1_x / ( ...
				R_T(squeeze(sum(M2_p1_x_List, 3))) / (M * NUM) ...
			);
		a_s1_x = real(a_s1_x);

		a_s2_x = c_p2_x / ( ...
				R_T(squeeze(sum(M2_p2_x_List, 3))) / (M * NUM) ...
			);
		a_s2_x = real(a_s2_x);

		w_s1_x = v_s_z + w_s_z * a_s1_x;
		w_s2_x = v_s_z + w_s_z * a_s2_x;

		res_x = 0.0;
		for i = 1 : NUM
			M2_p1_a = squeeze(M2_p1_a_List(:, :, i));
			M2_p2_a = squeeze(M2_p2_a_List(:, :, i));

			M2 = 1 / w_s2_x * M2_p2_a + 1 / w_s1_x * M2_p1_a;
			M2_eig = real(eig(M2));
			tmp_x = mean(v_p_tx ./ (v_p_tx * M2_eig + 1));
			res_x = res_x + tmp_x;
		end
		Hv_s_tx = res_x / NUM;
		Hv_s_tx = max(Hv_s_tx, EPS);

		% v_s_tx = 1 / (1 / Hv_s_tx - 1 / v_p_tx);
		v_s_tx = (Hv_s_tx * v_p_tx) / (v_p_tx - Hv_s_tx);
		v_s_tx = max(v_s_tx, EPS);

		% Forward TX
		% The message of TX
		Hv_p_tx = In_X.MSE_SE(v_s_tx);
		Hv_p_tx = real(Hv_p_tx);
		Hv_p_tx = max(Hv_p_tx, EPS);
		MSE_X_List(iter, 1) = Hv_p_tx / TX;

		% v_p_tx = 1 / (1 / Hv_p_tx - 1 / v_s_tx);
		v_p_tx = (Hv_p_tx * v_s_tx) / (v_s_tx - Hv_p_tx);

		v_p_tx_inv = 1 / v_p_tx;
		Idx_x = v_p_tx_inv < 1;

		v_p_tx(Idx_x) = v_p_tx_O(Idx_x);

		% Damping
		[v_p_tx, v_p_tx_O] = Damping(v_p_tx, v_p_tx_O, mesS);

		% Generate m_p_tx
		Noise_x = randn(N, 1, 1, NUM) * sqrt(v_s_tx);
		m_s_tx = x + Noise_x;

		[HM_PTX, HV_PTX] = In_X.Estimation(m_s_tx, v_s_tx * ones(N, 1, 1, NUM));
		HV_PTX = real(HV_PTX);
		HV_PTX = max(HV_PTX, EPS);

		V_PTX_inv = (v_s_tx - HV_PTX) ./ (HV_PTX * v_s_tx);
		m_p_tx = (HM_PTX * v_s_tx - m_s_tx .* HV_PTX) ./ (v_s_tx - HV_PTX);

		IdxX = V_PTX_inv < 1;

		m_p_tx(IdxX) = m_p_tx_O(IdxX);

		% Damping m_p_tx
		[m_p_tx, m_p_tx_O] = Damping(m_p_tx, m_p_tx_O, mesB);

		% Forward x
		% The message of x
		% Update c_p1_x
		% Update c_p2_x
		rex_x1 = 0.0;
		rex_x2 = 0.0;
		for i = 1 : NUM
			M2_p1_a = squeeze(M2_p1_a_List(:, :, i));
			M2_p2_a = squeeze(M2_p2_a_List(:, :, i));

			M2_p1_a_S = squeeze(M2_p1_a_S_List(:, :, i));
			M2_p2_a_S = squeeze(M2_p2_a_S_List(:, :, i));

			M2_1 = 1 / w_s2_x * M2_p2_a_S + 1 / w_s1_x * M2_p1_a;
			M2_1_eig = real(eig(M2_1));
			tmp_x1 = mean(v_p_tx ./ (v_p_tx * M2_1_eig + 1));
			rex_x1 = rex_x1 + tmp_x1;

			M2_2 = 1 / w_s2_x * M2_p2_a + 1 / w_s1_x * M2_p1_a_S;
			M2_2_eig = real(eig(M2_2));
			tmp_x2 = mean(v_p_tx ./ (v_p_tx * M2_2_eig + 1));
			rex_x2 = rex_x2 + tmp_x2;

		end
		c_p1_x = rex_x1 / NUM;
		c_p2_x = rex_x2 / NUM;

		c_p1_x_inv = 1 / c_p1_x;
		Idx_x1 = c_p1_x_inv < 1;
		c_p1_x(Idx_x1) = c_p1_x_O(Idx_x1);
		% Damping
		[c_p1_x, c_p1_x_O] = Damping(c_p1_x, c_p1_x_O, mesS);

		c_p2_x_inv = 1 / c_p2_x;
		Idx_x2 = c_p2_x_inv < 1;
		c_p2_x(Idx_x2) = c_p2_x_O(Idx_x2);
		% Damping
		[c_p2_x, c_p2_x_O] = Damping(c_p2_x, c_p2_x_O, mesS);

		cx_1 = v_p_tx - c_p1_x;
		cx_2 = v_p_tx - c_p2_x;
		cx_1 = max(cx_1, 1e-3);
		cx_2 = max(cx_2, 1e-3);

		Noise_x1 = randn(N, 1, M, NUM) * sqrt(cx_1);
		Noise_x2 = randn(N, 1, M, NUM) * sqrt(cx_2);
		m_p1_x_List = m_p_tx + Noise_x1;
		m_p2_x_List = m_p_tx + Noise_x2;

		[m_p1_x_List, m_p1_x_List_O] = Damping(m_p1_x_List, m_p1_x_List_O, mesB);
		[m_p2_x_List, m_p2_x_List_O] = Damping(m_p2_x_List, m_p2_x_List_O, mesB);

		m2_p1_x_List = pagemtimes(m_p1_x_List, 'none', m_p1_x_List, 'transpose');
		m2_p2_x_List = pagemtimes(m_p2_x_List, 'none', m_p2_x_List, 'transpose');

		M2_p1_x_List = squeeze(sum(m2_p1_x_List, 3));
		M2_p2_x_List = squeeze(sum(m2_p2_x_List, 3));

		m_p1_a_List = pagemtimes(A_List, 'none', m_p1_x_List, 'none');
		m_p2_a_List = pagemtimes(A_List, 'transpose', m_p2_x_List, 'none');

		m2_p1_a_List = pagemtimes(m_p1_a_List, 'none', m_p1_a_List, 'transpose');
		m2_p2_a_List = pagemtimes(m_p2_a_List, 'none', m_p2_a_List, 'transpose');

		M2_p1_a_List = squeeze(sum(m2_p1_a_List, 3));
		M2_p2_a_List = squeeze(sum(m2_p2_a_List, 3));

		m_p1_a_S_List = m_p1_a_List(:, :, 2 : end, :);
		m_p2_a_S_List = m_p2_a_List(:, :, 2 : end, :);

		m2_p1_a_S_List = pagemtimes(m_p1_a_S_List, 'none', m_p1_a_S_List, 'transpose');
		m2_p2_a_S_List = pagemtimes(m_p2_a_S_List, 'none', m_p2_a_S_List, 'transpose');

		M2_p1_a_S_List = squeeze(sum(m2_p1_a_S_List, 3));
		M2_p2_a_S_List = squeeze(sum(m2_p2_a_S_List, 3));

		% Forward Z
		v_p_z = 1 / M * c_p1_x * c_p2_x * R_T(A2_Sum) + ...
			1 / (M * NUM) * c_p1_x * R_T(squeeze(sum(M2_p2_a_List, 3))) + ...
			1 / (M * NUM) * c_p2_x * R_T(squeeze(sum(M2_p1_a_List, 3)));

		v_p_z_inv = 1 / v_p_z;
		Idx_z = v_p_z_inv < 1;

		v_p_z(Idx_z) = v_p_z_O(Idx_z);

		% Damping
		[v_p_z, v_p_z_O] = Damping(v_p_z, v_p_z_O, mesB);

		'----';

	end
end