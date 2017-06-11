%% showResult: function description
function [outputs] = showResult(dx, HS, Ec, J, dU, Tr)

	% figure('units', 'normalized', 'outerposition', [0 0 1 1]);

	% Axes = {
	% 	subplot(1, 2, 1);
	% 	subplot(1, 2, 2);
	% };

	% plot(Axes{1}, (0:HS-1)*dx, Ec);
	% plot(Axes{2}, log(Tr)', (ones(length(Ec(:, 1)), 1)*linspace(0, max(max(Ec)), length(Tr(1, :))))' );

	Axes = {
		subplot(1, 1, 1);
	};

	plot(Axes{1}, dU, J(1, :));

	outputs = 1;