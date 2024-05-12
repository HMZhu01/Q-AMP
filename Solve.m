function [MSE_X_List, MSE_X1_List, iterations, x_t] = Solve(Input, Obj)
	% x: N * 1, z: M * 1

	M = Input.M;
	N = Input.N;

	% arbitrary initial distance to between consecutive iterates
	diff = 1e6; 
	% iterations counter
	% maximum gradient descent iterations
	max_iterations = Input.max_iterations;
	% max_iterations = 2500;
	% max_iterations = 200;
	% maximum allowable distance between consecutive iterates
	diff_tol = 1e-6;

	MSE_X_List = zeros(max_iterations, 1);
	MSE_X1_List = zeros(max_iterations, 1);

	x = Obj.x;
	A_List = Obj.A_List;	
	y = Obj.y;

	% Initialization
	% N * 1
	x_0 = Get_Initialisation(A_List, y);

	MSE_X_List(1, 1) = norm(x_0 - x, 'fro').^(2) / norm(x, 'fro').^(2);
	MSE_X1_List(1, 1) = Dist(x_0, x);

	% step size
	lr = 1e-1;
	norm_estimate = real( ...
			( ...
				y' * y / double(2 * M) ...
			)^(1 / 4) ...
		);
	% scale the step size by the inverse of the norm squared as in the mathematical proofs
	lr = lr / norm_estimate^(2);

	x_t = x_0;

	iterations = 2;

	% gradient descent
	while (diff > diff_tol) && (iterations <= max_iterations)
		grad = Loss_Grad(x_t, A_List, y);
		Ox_t = x_t;
		x_t = x_t - lr * grad;

		% calculate error and print information
		MSE_X_List(iterations, 1) = norm(x_t - x, 'fro').^(2) / norm(x, 'fro').^(2);
		MSE_X1_List(iterations, 1) = Dist(x_t, x);

		diff = C_Dist(x_t, Ox_t) / norm_estimate;
		iterations = iterations + 1;
	end

	Start = min(iterations, max_iterations);

	MSE_X_List(Start : end, 1) = norm(x_t - x, 'fro').^(2) / norm(x, 'fro').^(2);
	MSE_X1_List(Start : end, 1) = Dist(x_t, x);
end