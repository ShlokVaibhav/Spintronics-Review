clear all
close all
%% MTJ

% constants
hbar = 1.06e-34;
q = 1.6e-19;
m = 9.1e-31;
kt = 0.025;
a = 1e-10;
sigma_z = [1 0; 0 -1];
sigma_y = [0 -i; i 0];
sigma_x = [0 1; 1 0];

Ub =3;
% fixed magnet magnetization
Mx = 0;
My = 0;
Mz = 1;
kkk = 0*2*pi*1e7*(1:5);
delta = 2.15;
Ef = 2.25;
m_FM = 0.73*m;
m_ox = 0.2*m;
t_FM = hbar^2/2/m_FM/a^2/q;
t_ox = hbar^2/2/m_ox/a^2/q;
Vapp = linspace(-1, 1, 50);
I = zeros(size(Vapp));
J_L = I;
J_R = I;
taup = I;
tauperp = I;
Nc = 15;
NL = 25;
NR= 25;
H = zeros(2*(NR+NL+Nc+2), 2*(NR+NL+Nc+2));
zplus = 1e-12;

mz = 1;
my = 0;
mx = 0;


for i=1:length(Vapp)
k = kkk(1);
a_FM = 2*t_FM+hbar^2*k^2/2/m_FM/q;
a_ox = 2*t_ox+hbar^2*k^2/2/m_ox/q;
a_int = 0.5*(a_ox+a_FM);

V = Vapp(i);
E=(Ef-abs(V)/2-4*kt):0.01:(Ef+abs(V)/2+4*kt+Ub);
f1 = 1./(1+exp((E - Ef-V/2)./kt));
f2 = 1./(1+exp((E - Ef+V/2)./kt));

for ii = 1:NL
H = placemat(H, eye(2)*(a_FM+V/2)+(eye(2)-Mz*sigma_z-Mx*sigma_x-My*sigma_y)*delta/2, ii, ii);
H = placemat(H, -t_FM*eye(2), ii, ii+1);
H = placemat(H, -t_FM*eye(2), ii+1, ii);
end

H = placemat(H, eye(2)*(a_int+Ub/2+V/2)+(eye(2)-Mz*sigma_z-Mx*sigma_x-My*sigma_y)*delta/4, NL+1, NL+1);

for ii = 1:Nc
H = placemat(H, eye(2)*(a_ox+Ub+V*(1/2-ii/(Nc+1))), ii+NL+1, ii+NL+1);
H = placemat(H, -t_ox*eye(2), ii+NL, ii+NL+1);
H = placemat(H, -t_ox*eye(2), ii+NL+1, ii+NL);
end

H = placemat(H, -t_ox*eye(2), Nc+NL+1, Nc+NL+2);
H = placemat(H, -t_ox*eye(2), Nc+NL+2, Nc+NL+1);

H = placemat(H, eye(2)*(a_int+Ub/2-V/2)+(eye(2)-mz*sigma_z-mx*sigma_x-my*sigma_y)*delta/4, NL+Nc+2, NL+Nc+2);

for ii = 1:NR
H = placemat(H, eye(2)*(a_FM-V/2)+(eye(2)-mz*sigma_z-mx*sigma_x-my*sigma_y)*delta/2, NL+Nc+2+ii, NL+Nc+2+ii);
H = placemat(H, -t_FM*eye(2), NL+Nc+1+ii, NL+Nc+2+ii);
H = placemat(H, -t_FM*eye(2), NL+Nc+2+ii, NL+Nc+1+ii);
end


k_up_La = acos(1-(E-V/2-hbar^2*k^2/2/m_FM/q)/2/t_FM);
k_up_Ra = acos(1-(E+V/2-hbar^2*k^2/2/m_FM/q)/2/t_FM);  
k_down_La = acos(1-(E-V/2-hbar^2*k^2/2/m_FM-delta)/2/t_FM); 
k_down_Ra = acos(1-(E+V/2-hbar^2*k^2/2/m_FM-delta)/2/t_FM);  

for iii = 1:length(E)
sigma_L= zeros(size(H));
sigma_L= placemat(sigma_L,[-t_FM*exp(j*k_up_La(iii)) 0; 0 -t_FM*exp(j*k_down_La(iii))],1,1);
sigma_R= zeros(size(H));
sigma_R= placemat(sigma_R,[-t_FM*exp(j*k_up_Ra(iii)) 0; 0 -t_FM*exp(j*k_down_Ra(iii))],NR+NL+Nc+2,NR+NL+Nc+2);

G = inv((E(iii)+j*zplus)*eye(length(H))-H - sigma_L-sigma_R);
A = j*(G-G');
Gamma1 = j*(sigma_L - sigma_L');
Gamma2 = j*(sigma_R - sigma_R');
T = real(trace(Gamma1*G*Gamma2*G'));
I(i) = I(i)+2*q^2/2/pi/hbar*T*(f1(iii)-f2(iii))*(E(2)-E(1));


Gn = G*Gamma1*G'*f1(iii)+G*Gamma2*G'*f2(iii);
endi = NR+NL+Nc+2;
endii = Nc+NL+1;
Gnd = Gn';
Hd = H';
J_R(i) = J_R(i)+j*(E(2)-E(1))*(trace(sigma_z*(H(endi-2:endi-1, endi:endi+1)*Gn(endi:endi+1, endi-2:endi-1)-Gnd(endi-2:endi-1, endi:endi+1)*Hd(endi:endi+1, endi-2:endi-1))));
J_L(i) = J_L(i)+j*(E(2)-E(1))*(trace(sigma_z*(H(endii-2:endii-1, endii:endii+1)*Gn(endii:endii+1, endii-2:endii-1)-Gnd(endii-2:endii-1, endii:endii+1)*Hd(endii:endii+1, endii-2:endii-1))));
end
temp = cross([mx my mz],cross([mx my mz], [0 0 J_L(i)]));
taup(i) = abs(temp(1)^2+temp(3)^2);
tauperp(i) = temp(2);
i
J_L(i)

end
plot(Vapp, I, 'LineWidth',1.2)
xlabel('V_{app} (V)')
ylabel('I (Arb. units)')
set(gca, "linewidth", 1, "fontsize", 18);
grid on
figure
plot(Vapp(1:end-1),diff(I)*1000/(Vapp(2)-Vapp(1)))
xlabel('V_{app} (V)')
ylabel('dI/dV (Arb. units)')
set(gca, "linewidth", 1, "fontsize", 18);
grid on
figure
plot(Vapp(1:end-1), diff(tauperp)*1000/(Vapp(2)-Vapp(1)))
figure
plot(Vapp(1:end-1), diff(taup)*1000/(Vapp(2)-Vapp(1)))

tempo = I;

clearvars -except Ub tempo
hbar = 1.06e-34;
q = 1.6e-19;
m = 9.1e-31;
kt = 0.025;
a = 1e-10;
sigma_z = [1 0; 0 -1];
sigma_y = [0 -i; i 0];
sigma_x = [0 1; 1 0];
% fixed magnet magnetization
Mx = 0;
My = 0;
Mz = 1;
kkk = 0*2*pi*1e7*(1:5);
delta = 2.15;
Ef = 2.25;
m_FM = 0.73*m;
m_ox = 0.2*m;
t_FM = hbar^2/2/m_FM/a^2/q;
t_ox = hbar^2/2/m_ox/a^2/q;
Vapp = linspace(-1, 1, 50);
I = zeros(size(Vapp));
J_L_x = I;
J_L_y = I;
J_L_z = I;
J_R = I;
tau = I;

Nc =15;
NL = 25;
NR = 25;
H = zeros(2*(NR+NL+Nc+2), 2*(NR+NL+Nc+2));
zplus = 1e-12;

mz = 0;
my = 0;
mx = 1;

for i=1:length(Vapp)
k = kkk(1);
a_FM = 2*t_FM+hbar^2*k^2/2/m_FM/q;
a_ox = 2*t_ox+hbar^2*k^2/2/m_ox/q;
a_int = 0.5*(a_ox+a_FM);

V = Vapp(i);
E=(Ef-abs(V)/2-4*kt):0.01:(Ef+abs(V)/2+4*kt+Ub);
f1 = 1./(1+exp((E - Ef-V/2)./kt));
f2 = 1./(1+exp((E - Ef+V/2)./kt));


for ii = 1:NL
H = placemat(H, eye(2)*(a_FM+V/2)+(eye(2)-Mz*sigma_z-Mx*sigma_x-My*sigma_y)*delta/2, ii, ii);
H = placemat(H, -t_FM*eye(2), ii, ii+1);
H = placemat(H, -t_FM*eye(2), ii+1, ii);
end

H = placemat(H, eye(2)*(a_int+Ub/2+V/2)+(eye(2)-Mz*sigma_z-Mx*sigma_x-My*sigma_y)*delta/4, NL+1, NL+1);

for ii = 1:Nc
H = placemat(H, eye(2)*(a_ox+Ub+V*(1/2-ii/(Nc+1))), ii+NL+1, ii+NL+1);
H = placemat(H, -t_ox*eye(2), ii+NL, ii+NL+1);
H = placemat(H, -t_ox*eye(2), ii+NL+1, ii+NL);
end

H = placemat(H, -t_ox*eye(2), Nc+NL+1, Nc+NL+2);
H = placemat(H, -t_ox*eye(2), Nc+NL+2, Nc+NL+1);

H = placemat(H, eye(2)*(a_int+Ub/2-V/2)+(eye(2)-mz*sigma_z-mx*sigma_x-my*sigma_y)*delta/4, NL+Nc+2, NL+Nc+2);

for ii = 1:NR
H = placemat(H, eye(2)*(a_FM-V/2)+(eye(2)-mz*sigma_z-mx*sigma_x-my*sigma_y)*delta/2, NL+Nc+2+ii, NL+Nc+2+ii);
H = placemat(H, -t_FM*eye(2), NL+Nc+1+ii, NL+Nc+2+ii);
H = placemat(H, -t_FM*eye(2), NL+Nc+2+ii, NL+Nc+1+ii);
end


k_up_La = acos(1-(E-V/2-hbar^2*k^2/2/m_FM/q)/2/t_FM);
k_up_Ra = acos(1-(E+V/2-hbar^2*k^2/2/m_FM/q)/2/t_FM);  
k_down_La = acos(1-(E-V/2-hbar^2*k^2/2/m_FM-delta)/2/t_FM); 
k_down_Ra = acos(1-(E+V/2-hbar^2*k^2/2/m_FM-delta)/2/t_FM);  

for iii = 1:length(E)
sigma_L= zeros(size(H));
sigma_L= placemat(sigma_L,[-t_FM*exp(j*k_up_La(iii)) 0; 0 -t_FM*exp(j*k_down_La(iii))],1,1);
alpha = -t_FM*exp(j*k_up_Ra(iii));
beta = -t_FM*exp(j*k_down_Ra(iii));

sigma_R= zeros(size(H));
randomm = 1/2*[(alpha+beta+(alpha-beta)*mz), ((alpha-beta)*(mx-j*my)); ((alpha-beta)*(mx+j*my)), (alpha+beta-(alpha-beta)*mz)];
sigma_R= placemat(sigma_R,randomm,NR+NL+Nc+2,NR+NL+Nc+2);


%sigma_R= placemat(sigma_R,[-t_FM*exp(j*k_up_Ra(iii)) 0; 0 -t_FM*exp(j*k_down_Ra(iii))],NR+NL+Nc+2,NR+NL+Nc+2);

G = inv((E(iii)+j*zplus)*eye(length(H))-H - sigma_L-sigma_R);
A = j*(G-G');
Gamma1 = j*(sigma_L - sigma_L');
Gamma2 = j*(sigma_R - sigma_R');
T = real(trace(Gamma1*G*Gamma2*G'));
I(i) = I(i)+2*q^2/2/pi/hbar*T*(f1(iii)-f2(iii))*(E(2)-E(1));

Gn = G*Gamma1*G'*f1(iii)+G*Gamma2*G'*f2(iii);
endi = NR+NL+Nc+2;
endii = Nc+NL+1;
Gnd = Gn';
Hd = H';
J_R(i) = J_R(i)+j*(E(2)-E(1))*(trace(sigma_z*(H(endi-2:endi-1, endi:endi+1)*Gn(endi:endi+1, endi-2:endi-1)-Gnd(endi-2:endi-1, endi:endi+1)*Hd(endi:endi+1, endi-2:endi-1))));
J_L_z(i) = J_L_z(i)+j*(E(2)-E(1))*(trace(sigma_z*(H(endii-2:endii-1, endii:endii+1)*Gn(endii:endii+1, endii-2:endii-1)-Gnd(endii-2:endii-1, endii:endii+1)*Hd(endii:endii+1, endii-2:endii-1))));
J_L_y(i) = J_L_y(i)+j*(E(2)-E(1))*(trace(sigma_y*(H(endii-2:endii-1, endii:endii+1)*Gn(endii:endii+1, endii-2:endii-1)-Gnd(endii-2:endii-1, endii:endii+1)*Hd(endii:endii+1, endii-2:endii-1))));
J_L_x(i) = J_L_x(i)+j*(E(2)-E(1))*(trace(sigma_x*(H(endii-2:endii-1, endii:endii+1)*Gn(endii:endii+1, endii-2:endii-1)-Gnd(endii-2:endii-1, endii:endii+1)*Hd(endii:endii+1, endii-2:endii-1))));

end
temp = cross([mx my mz],cross([mx my mz], [J_L_x(i) J_L_y(i) J_L_z(i)]));
i
[J_L_x(i) J_L_y(i) J_L_z(i)];
taup(i) = abs(temp(1)^2+temp(3)^2);
tauperp(i) = temp(2);

end
figure(1)
hold on
plot(Vapp, I)
figure(2)
hold on
plot(Vapp(1:end-1),diff(I)*1000/(Vapp(2)-Vapp(1)))
figure(3)
hold on
plot(Vapp(1:end-1), diff(tauperp)*1000/(Vapp(2)-Vapp(1)))
figure(4)
hold on
plot(Vapp(1:end-1), diff(taup)*1000/(Vapp(2)-Vapp(1)))



figure(1)
legend('parallel', 'antiparallel')
figure(2)
legend('parallel', 'antiparallel')
figure(3)
legend('parallel', 'antiparallel')

figure
plot(Vapp(1:end-1),100*diff(tempo)./diff(I)-100)
xlabel('V_{app} (V)')
ylabel('TMR')
set(gca, "linewidth", 1, "fontsize", 18);
grid on




function Hnew = placemat(H, block, n,m)
Hnew = H;
Hnew(2*(n-1)+1, 2*(m-1)+1) = block(1,1);
Hnew(2*(n-1)+2, 2*(m-1)+1) = block(2,1);
Hnew(2*(n-1)+1, 2*(m-1)+2) = block(1,2);
Hnew(2*(n-1)+2, 2*(m-1)+2) = block(2,2);
end
