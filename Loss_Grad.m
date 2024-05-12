function grad = Loss_Grad(x_t, A_List, y)
	% computes the gradient

	M = size(A_List, 3);
	N = size(x_t, 1);

	% equation (4)
	% grad = 0;
	% for i = 1 : M

	% 	yi = y(i, 1);
	% 	Ai = squeeze(A_List(:, :, i));

	% 	grad = grad + ( ...
	% 			conj(yi) - x_t.' * conj(Ai) * conj(x_t) ...
	% 		) * ( ...
	% 			-1 * Ai.' * conj(x_t) ...
	% 		) + ( ...
	% 			yi - x_t' * Ai * x_t ...
	% 		) * ( ...
	% 			-1 * conj(Ai) * conj(x_t) ...
	% 		);

	% end

	y_List = zeros(1, 1, M);
	y_List(1, :, :) = y';

	A_m = pagemtimes(A_List, 'none', x_t, 'none');
	At_m = pagemtimes(A_List, 'transpose', x_t, 'none');

	m_A_m = pagemtimes(x_t, 'transpose', A_m, 'none');

	tmp2_1 = (m_A_m - y_List) .* A_m;
	tmp2_2 = (m_A_m - y_List) .* At_m;

	tmp2 = sum(tmp2_1, 3) + sum(tmp2_2, 3);

	grad = squeeze(tmp2)  / double(2 * M);

end