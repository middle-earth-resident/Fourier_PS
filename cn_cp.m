N = 256;
N2 = N/2;
gamma = 1;
nu = 1;
dx = 0.05;
tmax = 20;
dt = 0.02;
nmax = round(tmax/dt);
for i=1:N,
	for j = 1:N,
		x(i) = (i-N2)*dx;
		y(i) = (i-N2)*dx;
	end;
end;
k = [0:N/2-1 -N/2:-1]*2*pi/(N*dx) ; 
k2 = k.^2;
x2 = x.*x;
y2 = y.*y;
qr_al = sqrt(sqrt(gamma*nu)/pi);
ci = 0+1i;
for J = 1:N,
	for I = 1:N,
				
		v(I,J) = (gamma*x2(I) + nu*y2(J))/2;
		tmp = (gamma*x2(I) + nu*y2(J))/2;
		cp(I,J) = qr_al * exp(-tmp) * (x(I)+ci*y(J)) * exp(2*ci*pi);
	end;
end;


p2 = abs(cp);
phase = atan(imag(cp)/real(cp));

[dtmpx,dtmpy] = gradient(phase);
dtmp = dtmpx+dtmpy;

con = sqrt(p2)*dtmp;
for i=1:N,
		divc(i,:) = divergence(con(i,:),con(:,i)');
end;

pcolor(x,y,abs(divc));
hold on;
shading("flat");
