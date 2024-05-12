function Res = C_Dist(x_t, x)
	a_t = x_t(:);
	a = x(:);

	Res = norm( ...
			a_t - exp(-j * angle(a_t' * a)) * a ...
		);
end