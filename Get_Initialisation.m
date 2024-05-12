function Res = Get_Initialisation(A_List, y)

	M = size(A_List, 3);

	y_List = zeros(1, 1, M);
	y_List(1, :, :) = y';

	List = y_List .* A_List;

	S = squeeze(sum(List, 3))  / double(M);

	[U, L, V] = svd(S);
	v = V(:, 1);

	Res = (y' * y / (2 * M))^(1 / 4) * v;
end