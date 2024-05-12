function Res = Dist(x_t, x)
	a_t = x_t(:);
	a = x(:);

	Res = C_Dist(x_t, x) / norm(a);
end