%% getBandPropAlGaAs: function description
function [Ec, Eg, me_eff, mp_eff] = getBandPropAlGaAs(C_Al)
	n_Atoms = 4.42*1e28;
	n_Al = n_Atoms/2;
	x = C_Al./n_Al;

	more45 = double(x >= 0.45);
	more45Ec = 0.232*more45 - 0.259*more45.*x + 1.147*more45.*x.^2;
	more45Eg = 1.656*more45 + 0.215*more45.*x + 0.143*more45.*x.^2;

	less45 = double(x < 0.45);
	less45Ec = 0.773*less45.*x;
	less45Eg = 1.424*less45 + 1.247*less45.*x;

	Ec = more45Ec + less45Ec;
	Ec = Ec - min(Ec);

	Eg = more45Eg + less45Eg;

	me_eff = 0.067 + 0.083*x;
	mp_eff = 0.45*ones(size(x));
end