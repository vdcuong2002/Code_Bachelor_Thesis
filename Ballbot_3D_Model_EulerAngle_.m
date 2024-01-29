clear
tic
%% Parameters
syms mb Jb rb ma Jax Jay Jaz l rw Jw g real
l = 0.37;
rw = 0.05;
rb = 0.12;
alpha = pi/4;
beta = 2*pi/3;

% Variables
syms x(t) y(t) phix(t) phiy(t) phiz(t)
syms dx(t) dy(t) dphix(t) dphiy(t) dphiz(t)

dx = diff(x);
dy = diff(y);
dphix = diff(phix);
dphiy = diff(phiy);
dphiz = diff(phiz);

ddx = diff(dx);
ddy = diff(dy);
ddphix = diff(dphix);
ddphiy = diff(dphiy);
ddphiz = diff(dphiz);

% q = [phix; phiy; phiz; q4];
% dq = [dphix; dphiy; dphiz; dq4];

PSI(t) = [x(t); y(t); phix(t); phiy(t); phiz(t)];
dPSI(t) = [dx(t); dy(t); dphix(t); dphiy(t); dphiz(t)];
ddPSI(t) = [ddx(t); ddy(t); ddphix(t); ddphiy(t); ddphiz(t)];

assumeAlso(PSI(t), 'real')
assumeAlso(dPSI(t), 'real')
assumeAlso(ddPSI(t), 'real')


% Inverse Kinetics
dC1 = [-sin(0*beta);
    cos(0*beta);
    0];
dC2 = [-sin(1*beta);
    cos(1*beta);
    0];
dC3 = [-sin(2*beta);
    cos(2*beta);
    0];

pC1 = [cos(0*beta)*sin(alpha)*rb;
    sin(0*beta)*sin(alpha)*rb;
    cos(alpha)*rb];
pC2 = [cos(1*beta)*sin(alpha)*rb;
    sin(1*beta)*sin(alpha)*rb;
    cos(alpha)*rb];
pC3 = [cos(2*beta)*sin(alpha)*rb;
    sin(2*beta)*sin(alpha)*rb;
    cos(alpha)*rb];


RX = [1 0 0; 
    0 cos(phix) -sin(phix); 
    0 sin(phix) cos(phix)];

RY = [cos(phiy) 0 sin(phiy); 
    0 1 0; 
    -sin(phiy) 0 cos(phiy)];
RZ = [cos(phiz) -sin(phiz) 0; 
    sin(phiz) cos(phiz) 0; 
    0 0 1];

R_a_to_b = RZ*RY*RX;

wb = 1/rb*[-dy; dx; 0];
wa_b = [dphix; dphiy; dphiz]; 
wa_a = R_a_to_b^(-1)*wa_b;
wb_a = R_a_to_b^(-1)*wb;

w_b_a_a = wb_a - wa_a;
vC1b = cross(w_b_a_a, pC1);
vC2b = cross(w_b_a_a, pC2);
vC3b = cross(w_b_a_a, pC3);


dtheta1 = 1/rw * dC1' * vC1b;
dtheta2 = 1/rw * dC2' * vC2b;
dtheta3 = 1/rw * dC3' * vC3b;


% Ball energy
vb = [dx; dy; 0];
Jb_ = diag([Jb, Jb, Jb]);
% Tb = 1/2 * mb * (vb'*vb) + 1/2 * wb' * Jb_ * wb;
Tb = 1/2*(mb + Jb/rb^2)*(dx^2 + dy^2);
Vb = 0;


% Body Energy
% derivative of rotation matrix
dR_a_to_b =  [        - dphiy*sin(phiy) - dphiz*cos(phiy)*sin(phiz),                                   dphiy*cos(phiy)*sin(phix) - dphiz*(cos(phix)*cos(phiz) + sin(phix)*sin(phiy)*sin(phiz)),                                     dphiz*(cos(phiz)*sin(phix) - cos(phix)*sin(phiy)*sin(phiz)) + dphiy*cos(phix)*cos(phiy);
          dphix*sin(phiy) + dphiz*cos(phiy)*cos(phiz),                                 - dphiz*(cos(phix)*sin(phiz) - cos(phiz)*sin(phix)*sin(phiy)) - dphix*cos(phiy)*sin(phix),                                     dphiz*(sin(phix)*sin(phiz) + cos(phix)*cos(phiz)*sin(phiy)) - dphix*cos(phix)*cos(phiy);
dphix*cos(phiy)*sin(phiz) - dphiy*cos(phiy)*cos(phiz), dphix*(cos(phix)*cos(phiz) + sin(phix)*sin(phiy)*sin(phiz)) + dphiy*(cos(phix)*sin(phiz) - cos(phiz)*sin(phix)*sin(phiy)), - dphix*(cos(phiz)*sin(phix) - cos(phix)*sin(phiy)*sin(phiz)) - dphiy*(sin(phix)*sin(phiz) + cos(phix)*cos(phiz)*sin(phiy))];
 

pb = [x; y; 0];
vb = [dx; dy; 0];
pCOM_a = [0; 0; l];
pCOM_b = R_a_to_b * pCOM_a;
vCOM_b = dR_a_to_b * pCOM_a;

pCOM = pb + pCOM_b;
vCOM = vb + vCOM_b;
h = [0 0 1] * pCOM;
Ja = diag([Jax, Jay, Jaz]);

Ta = 1/2 * ma * (vCOM' * vCOM) + 1/2 * wa_a' * Ja * wa_a;
% Ta = 1/2 * wa_a' * Ja * wa_a;
% Ta = 1/2 * ma * (vCOM' * vCOM);
% Ta = 0;
Va = ma*g*h;

% Wheels Energy
Jw_ = diag([Jw, Jw, Jw]);
dtheta = [dtheta1; dtheta2; dtheta3];
Tw = 1/2 * dtheta' * Jw_ * dtheta;
Tw = simplify(Tw);
Vw = 0;

% Tw = 0;
% Tb = 0;

% Ta = 0;
% Va = 0;
L = Tb + Ta + Tw - Va;
TTT = Tb + Ta + Tw + Va;
L = simplify(L);
TTT = simplify(TTT);
%% Calculate derivate Lagrange function
% Calculate dL/dpsi
syms dx(t) dy(t) dphix(t) dphiy(t) dphiz(t)
% define dpsi
dPSI_ = [dx(t); dy(t); dphix(t); dphiy(t); dphiz(t)];
ddPSI_ = diff(dPSI_, t);
assumeAlso(dPSI_, 'real')       % contraint dPSI_
assumeAlso(ddPSI_, 'real') 
% replace diff() in L by d_(t) in dPSI_
LL = subs(L, dPSI, dPSI_);
eqn = diff(functionalDerivative(LL, dPSI_), t) - functionalDerivative(LL, PSI);


%% Re-define variable
syms x y phix phiy phiz real
syms dx dy dphix dphiy dphiz real
syms ddx ddy ddphix ddphiy ddphiz real

psi= [x; y; phix; phiy; phiz];
dpsi= [dx; dy; dphix; dphiy; dphiz];
ddpsi= [ddx; ddy; ddphix; ddphiy; ddphiz];


eqn_sym = subs(eqn, [PSI; dPSI; ddPSI; dPSI_; ddPSI_], [psi; dpsi; ddpsi; dpsi; ddpsi]);
% L = subs(L, [PSI; dPSI; ddPSI; dPSI_; ddPSI_], [psi; dpsi; ddpsi; dpsi; ddpsi]);
% matlabFunction(L(0), "File", "output_E/Lagranger_function", "Vars", [x,y,phix,phiy,phiz,dx,dy,dphix,dphiy,dphiz,rb,rw,l,Jb,Jax,Jay,Jaz,Jw,ma,mb,g]);
% TTT = subs(TTT, [PSI; dPSI; ddPSI; dPSI_; ddPSI_], [psi; dpsi; ddpsi; dpsi; ddpsi]);
% matlabFunction(TTT(0), "File", "output_E/E_function", "Vars", [x,y,phix,phiy,phiz,dx,dy,dphix,dphiy,dphiz,rb,rw,l,Jb,Jax,Jay,Jaz,Jw,ma,mb,g]);

eq = simplify(eqn_sym);
% eq = eq(0);
% eq = collect(eq, ddpsi);

%% acsadf
M = sym([]);
% Row 1
eqM = [1 0 0 0 0 ]*eq;
M1 = [diff(eqM, ddx)    diff(eqM, ddy)    diff(eqM, ddphix)   ...
    diff(eqM, ddphiy)    diff(eqM, ddphiz)];
% Row 2
eqM = [0 1 0 0 0 ]*eq;
M2 = [diff(eqM, ddx)    diff(eqM, ddy)    diff(eqM, ddphix)   ...
    diff(eqM, ddphiy)    diff(eqM, ddphiz)];
% Row 3
eqM = [0 0 1 0 0 ]*eq;
M3 = [diff(eqM, ddx)    diff(eqM, ddy)    diff(eqM, ddphix)   ...
    diff(eqM, ddphiy)    diff(eqM, ddphiz)];
% Row 4
eqM = [0 0 0 1 0 ]*eq;
M4 = [diff(eqM, ddx)    diff(eqM, ddy)    diff(eqM, ddphix)   ...
    diff(eqM, ddphiy)    diff(eqM, ddphiz)];
% Row 5
eqM = [0 0 0 0 1 ]*eq;
M5 = [diff(eqM, ddx)    diff(eqM, ddy)    diff(eqM, ddphix)   ...
    diff(eqM, ddphiy)    diff(eqM, ddphiz)];
% Row 6

% Combine
M = [M1; M2; M3; M4; M5];
M = simplify(M);

% np = 5;
% M = sym([]);
% for (m = 1:np)
%     for (n = 1:np)
%         z = zeros(1,np);
%         z(m) = 1;
%         tmp = coeffs(z*collect(eq, ddpsi), ddpsi(n), 'all');
%         if (length(tmp) < 2)
%             tmp = zeros(1,2);
%         end       
%         M(m,n) = [1 0]*tmp'; % take end element
%     end
% end
% 
% M
% Gravity matrix, size = 6x1
% replace all derivate symbolic by zero
G = subs(eq, [dpsi; ddpsi], zeros(10, 1));
Cq = simplify(eq - M*ddpsi - G);
%%
dtheta1 = subs(dtheta1, [PSI; dPSI; ddPSI; dPSI_; ddPSI_], [psi; dpsi; ddpsi; dpsi; ddpsi]);
dtheta2 = subs(dtheta2, [PSI; dPSI; ddPSI; dPSI_; ddPSI_], [psi; dpsi; ddpsi; dpsi; ddpsi]);
dtheta3 = subs(dtheta3, [PSI; dPSI; ddPSI; dPSI_; ddPSI_], [psi; dpsi; ddpsi; dpsi; ddpsi]);
dtheta1 = simplify(dtheta1(0)); 
dtheta2 = simplify(dtheta2(0)); 
dtheta3 = simplify(dtheta3(0));

J = [diff(dtheta1, dx)      diff(dtheta2, dx)       diff(dtheta3, dx);
    diff(dtheta1, dy)      diff(dtheta2, dy)       diff(dtheta3, dy);
    diff(dtheta1, dphix)      diff(dtheta2, dphix)       diff(dtheta3, dphix);
    diff(dtheta1, dphiy)      diff(dtheta2, dphiy)       diff(dtheta3, dphiy);
    diff(dtheta1, dphiz)      diff(dtheta2, dphiz)       diff(dtheta3, dphiz)];

dtheta = [dtheta1; dtheta2; dtheta3];
%% Export
% matlabFunction(M, "File", "output\EulerAngle_m_matrix", "Vars", [x,y,phix,phiy,phiz,dx,dy,dphix,dphiy,dphiz,rb,rw,l,Jb,Jax,Jay,Jaz,Jw,ma,mb,g]);
% matlabFunction(Cq, "File", "output\EulerAngle_cq_matrix", "Vars", [x,y,phix,phiy,phiz,dx,dy,dphix,dphiy,dphiz,rb,rw,l,Jb,Jax,Jay,Jaz,Jw,ma,mb,g]);
% matlabFunction(G, "File", "output\EulerAngle_g_matrix", "Vars", [x,y,phix,phiy,phiz,dx,dy,dphix,dphiy,dphiz,rb,rw,l,Jb,Jax,Jay,Jaz,Jw,ma,mb,g]);
% matlabFunction(J, "File", "output\EulerAngle_j_matrix",  "Vars", [x,y,phix,phiy,phiz,dx,dy,dphix,dphiy,dphiz,rb,rw,l,Jb,Jax,Jay,Jaz,Jw,ma,mb,g]);
% matlabFunction(dtheta, "File", "output_E/InverseKineticEulerAngle_v6");
toc