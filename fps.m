
L = 80 ; N =256 ; dt = 0.02 ; tmax = 20 ; nmax = round(tmax/dt);

dx = L/N ; x = [-L/2 : dx : L/2-dx]; k = [0:N/2-1 -N/2:-1]*2*pi/L ; k2 = k.^2;
u = 1.2*sech(1.2*(x+20)).* exp(i * x) + 0.8 * sech(0.8 * x);
udata = u ; tdata = 0;
for nn = 1:nmax
	du1 = i*(ifft(-k2.*fft(u))+2*u.*u.*conj(u)); v = u+0.5*du1*dt;
	
	du2 = i*(ifft(-k2.*fft(v))+2*v.*v.*conj(v)); v = u+0.5*du2*dt;
	du3 = i*(ifft(-k2.*fft(v))+2*v.*v.*conj(v)); v = u+du3*dt;
	du4 = i*(ifft(-k2.*fft(v))+2*v.*v.*conj(v));
	u = u+(du1+2*du2+2*du3+du4)*dt/6;
	size(udata,1)
	if (mod(nn,round(nmax/25)) == 0),
		udata = [udata ; u] ; 
		tdata = [tdata nn*dt];
	end;
end;

waterfall(x,tdata,abs(udata));
colormap([0 0 0]) ; view(10,60);


