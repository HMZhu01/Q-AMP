classdef In_Real_0_1_Estimation

	properties
		q;
		All_Symbol;
		TX;
	end

	methods
		% Constructor
		% q is the probability of 1.
		function Obj = In_Real_0_1_Estimation(q)
			Obj.q = q;
			Obj.All_Symbol = [0, 1];
			Obj.TX = q;
		end

		function X = Generation(Obj, size)
			q = Obj.q;

			X = double(rand(size) > 1 - q);
		end

		function [Hmx, Hvx] = Estimation(Obj, mx, vx, Init_Flag, Init_X)
			q = Obj.q;

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
					throw(MException('Foo:FatalError', 'In_Real_0_1_Estimation.Estimation(...)'));
			end

			if ~ Flag

				Size = size(mx);

				vx = reshape(vx, numel(mx), 1);
				mx = reshape(mx, numel(mx), 1);

				tmp = q * Real_Gaussian(mx, 1, vx);
				C = (1 - q) * Real_Gaussian(mx, 0, vx) + tmp;

				Hmx = tmp ./ C;
				Hvx = Hmx - Hmx.^(2);

				Hvx = max(Hvx, eps);

				Hmx = reshape(Hmx, Size);
				Hvx = reshape(Hvx, Size);
			else
				Hmx = x;
				Hvx = zeros(size(x));
			end

		end

		function Hvx = MSE_SE(Obj, vx)

			q = Obj.q;
			TX = Obj.TX;

			q_tx = integral(@(u) Func1(u, q, vx), - Inf, Inf);

			Hvx = TX - q_tx;

		end
	end
end

% clc;
% clear;

% addpath('../../Operator');

% n = 1000;

% q_List = rand(1, n);
% v_List = rand(1, n);

% q_x_List = zeros(1, n);
% q_x2_List = zeros(1, n);

% for i = 1 : n

% 	q = q_List(1, i);
% 	v = v_List(1, i);

% 	q_x1 = integral(@(u) ...
% 		Func1(u, q, v), - Inf, Inf ...
% 	);
% 	q_x1_List(1, i) = q_x1;

% 	q_x2 = integral(@(u) ...
% 		Func2(u, q, v), - Inf, Inf ...
% 	);
% 	q_x2_List(1, i) = q_x2;

% end

% Error = norm(q_x1_List - q_x2_List, 'fro').^(2) / norm(q_x2_List, 'fro').^(2);
% disp(['Error: ', num2str(Error)]);

% q_x_List = [q_x1_List; q_x2_List];

function Res = Func1(u, q, v)
	% disp('Func1: ');

	Res = ( ...
			q^(2) * exp(u / sqrt(v) - 1 / (2 * v)) ...
		) ./ ( ...
			(1 - q) * exp(1 / (2 * v) - u / sqrt(v)) + q ...
		) .* Real_Gaussian(u, 0, 1);

	Index1 = isinf(Res);
	Res(Index1) = 0.0;
	% disp(['Index1:', num2str(sum(double(Index1)))]);

	Index2 = isnan(Res);
	Res(Index2) = 0.0;
	% disp(['Index2:', num2str(sum(double(Index2)))]);

	'----';
end

% function Res = Func2(u, q, v)
% 	% disp('Func2: ');

% 	Res = ( ...
% 			q * Real_Gaussian(u, 1, v) ...
% 		).^(2) ./ ( ...
% 			(1 - q) * Real_Gaussian(u, 0, v) + q * Real_Gaussian(u, 1, v) ...
% 		);

% 	Index1 = isinf(Res);
% 	Res(Index1) = 0;
% 	% disp(['Index1:', num2str(sum(double(Index1)))]);

% 	Index2 = isnan(Res);
% 	Res(Index2) = 0;
% 	% disp(['Index2:', num2str(sum(double(Index2)))]);

% 	'----';
% end