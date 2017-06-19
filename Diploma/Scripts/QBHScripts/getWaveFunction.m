% getWaveFunction: function description
function [waveLeft, waveRigth] = getWaveFunction(delta, meff, U, Ez)
	hbar = 1.0551*1e-34;

	EzLen = length(Ez);
	ULen = length(U);

	waveLeft = zeros(EzLen, ULen);
	waveRigth = zeros(EzLen, ULen);

	for j = 1 : EzLen
		kLeft = sqrt( 2*meff(1)*(Ez(j) - U(1)) )/hbar;
		kRight = sqrt( 2*meff(end)*(Ez(j) - U(end)) )/hbar;

		d1 = ones(ULen-1, 1);
		d2 = 2*delta^2*meff(2:end-1).*( Ez(j)-U(2:end-1) )./hbar^2 - meff(2:end-1)./meff(3:end) - 1;
		d2 = [1i*kLeft*delta - 1, d2, 1i*kRight*delta - 1];
		d3 = [1, meff(2:end-1)./meff(3:end)];

		H = diag(d1, -1) + diag(d2) + diag(d3, +1);

		fLeft = [2*1i*kLeft*delta; zeros(ULen-1, 1)];
		fRight = [zeros(ULen-1, 1); 2*1i*kRight*delta];

		waveLeft(j, :) = (inv(H)*fLeft)';
		waveRigth(j, :) = (inv(H)*fRight)';
	end
	
end