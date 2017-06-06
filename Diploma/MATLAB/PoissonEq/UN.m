function [Vnew,nold]=UN(epsv,Ec,mi,Ni,eps,delta,V,i1,i2)
	qe=1.602E-19;
	
	n=length(Ec);
	
	Vnew=[zeros(1,i1-1),linspace(0,V,i2-i1+1),V*ones(1,n-i2)];
	Vold=Vnew+10;
	
	while (max(abs((Vnew-Vold)))>epsv)
		Vold=Vnew;
		Ui=Ec-qe*Vold;
		nold=nz(Ui,mi,delta,i1,i2);
        Vnew=pois(V,Vold,nold,eps,Ni,delta);
	end

	Ui = Ec-qe*Vnew;
	nold = nz(Ui,mi,delta,i1,i2);
end