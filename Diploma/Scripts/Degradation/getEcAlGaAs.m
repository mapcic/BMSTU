%% getEc_AlGaAs: function description
function [Ec, meff] = getEcAlGaAs(C_Al)
	n_Atoms = 4.42*1e28;
	n_Al = n_Atoms/2;
	x = C_Al./n_Al;

	more45 = double(x >= 0.45);
	more45 = 0.232*more45 - 0.259*more45.*x + 1.147*more45.*x.^2;

	less45 = double(x < 0.45);
	less45 = 0.773*less45.*x;

	Ec = (more45 + less45)./2;
	Ec = Ec - min(Ec);

	meff = 0.067 + 0.083*x;
end