%------------------------------------------------------------------
% Inner estimates of the (uncertain) region of attraction (ERA) with Invariant Sets
%   Compute estimates of the region of attraction using the Invariant Set formulation. Threee different
%   iterative algorithms are available (nomenclature from the journal paper
%   "Robust estimations of the region of attraction using invariant sets" (The Journal Of Franklin Institute)
% 
%   NOTE: the suite of libraries "SOSAnalysis" downloadable at 
%   http://www.aem.umn.edu/~AerospaceControl/
%   must be loaded before starting the analyses.
%   e.g.  addpath(genpath('...\Software\SOSAnalysis'))
%   Specifically: 
%   - the library "multipoly" is used for manipulation of polynomials;
%   - the library "sosopt" is used for SOS optimization 
%   
% Author: Andrea Iannelli (University of Bristol)
%------------------------------------------------------------------

clear all; clc; close all;






% % % % System selection
% 1 -> Van der Pol
% 2 -> Short Period

system = 1;


% Create vector field

switch system
     
     case 1 
% Van Der Pol Oscillators
         
pvar x1 x2 delta;

x = [x1;x2];
x_delta = [x1;x2;delta];
x1dot = -x2 + 0.2 * delta *x2;
x2dot = x1+(x1^2-1)*x2;
f = [x1dot; x2dot];

f_nom=subs ( f ,delta,0 );


     case 2
         
% Short Period
         
pvar x1 x2 x3 x4 x5 delta_1 delta_2;
x = [x1;x2;x3;x4;x5]; 
delta = [ delta_1; delta_2 ];
x_delta = [x;delta];

% [z1,z2,eta1,eta2,z3]; % paper nomenclature
% z1= pitch rate; z2 = angle of attack; z3 = pitch angle; \eta1,2= controller states

f1 =  3.0000*x(4)+1.3489*x(3)-.56250*x(5)-1.3489*x(2)-3.0000*x(1)-.36992e-1*x(2)*x(3)-.82272e-1*x(2)*x(4) + (1+delta_1) *(.82272e-1*x(1)*x(2)+.43981*x(2)^2+.15426e-1*x(2)*x(5)+.22335*x(2)^3);
    
f2 = .91136*x(1)-.64450*x(2)-.16621e-1*x(5)+.39857e-1*x(3)+.88643e-1*x(4) + (1+delta_2) *(-.54444e-1*x(2)^2-.54444e-1*x(5)^2+.10889*x(2)*x(5));
    
f3 = -.60401e-1*x(1)-.16621e-1*x(5)-.60464*x(3)+.88643e-1*x(4);
f4 = -.75000*x(1)-.28125*x(5);
f5 = x(1);
f = [f1;f2;f3;f4;f5];


f_nom = subs ( f ,delta,[0;0] );

 


end

n_x = size(x,1);


% Lyapunov Function of the linearized equilibrium
Q = eye(n_x);

[Vlin,~,~]=linstab(f_nom,x,Q);

R_0 = Vlin;


% Initialization for R



%%%%% Create level set functions variables (with multipoly utilities)

% V_n
VnDeg=4;
para_V = 2;
switch para_V
     
     case 1 
% Two options:
% para_V=1: V independent of uncertainties (cheaper, but possibly more
% conservative)
zVn = monomials(x,1:VnDeg);
[V_n,~] = polydecvar('Vn',zVn); % generic definition of poly

     case 2 

% para_V=2: V the monomial basis of V_n is such that contains, for the defined VnDeg:
% - all the possible monomials in x (gathered initially in zVn_1)

% -  all the possible monomials x*delta (added later in zVn_b)


zVn_1 = monomials(x,1:VnDeg);

% -  all the possible monomials x*delta
zVn_b=zVn_1;
for jj = 1 : VnDeg-1

zVn_1d = monomials(delta,jj);
zVn_1x = monomials(x,1:VnDeg-jj);

size_d = max(size(zVn_1d));

for ii = 1:size_d
    zVn_21 = zVn_1x*zVn_1d(ii);

  zVn_b = [zVn_b;zVn_21];  
    
end

end
 zVn=zVn_b;   


[V_n_2_2,~] = polydecvar('Vn',zVn); % generic definition of poly

%    keyboard
V_n = V_n_2_2;

end

% R_n
RnDeg = VnDeg;
zRn = monomials(x,1:RnDeg);
[R_n,~] = polydecvar('Rn',zRn); % generic definition of poly







% Definition of the polynomials used to describe the uncertainty set as
% semialgebraic set

switch system
     
     case 1 
% Van Der Pol Oscillators
         


% uncertainty bound : [-1;1]
d_l = -1; d_u = 1; 

h = 1/2*(d_l+d_u);
k =  (d_u-h)^2;


m_c = -(delta-h).^2+k; % (s_2 in the paper (Eq. 22))



     case 2

% Short Period

%  The "cheap" option which employs a single multiplier
%  m_c(\delta_1,\delta_2) is implemented here. 

% uncertainties bounds : [-0.1;0.1]
d_l = -.1; d_u = .1; 

h = 1/2*(d_l+d_u);
k =  (d_u-h)^2;


m_d1 = -(delta_1-h).^2+k; % (s_2 in the paper (Eq. 22))
m_d2 = -(delta_2-h).^2+k; % (m_2 in the paper (Eq. 22))

m_c = m_d1 + m_d2; %  (m_c in the paper (Eq. 23))
         
 
end





% % Multipliers parameterization (results might be sensititive; these are
% recommented options, but it is interesting to study the effect of
% changing the parameterization)

switch system
     
    
     case 1 
         
         % Van Der Pol Oscillators



% s_1,s_2,s_3 in the paper
zmDeg=1;
zm = monomials(x, 0:zmDeg );
s_1 = sosdecvar('m0',zm);

s_2 = sosdecvar('m1',zm);

s_3 = sosdecvar('m3',zm);

% s_0 in the paper
pDeg=4; % 
zp = monomials(x,0:pDeg);
[s_0,~] = polydecvar('s_0',zp); %


%%%%% Create \delta-multipliers variables (they will feature in the
%%%%% polynomial \Gamma Eq. 18d)

zmdDeg=1;
zmd = monomials(x_delta, 0:zmdDeg );

s_01 = sosdecvar('md2',zmd);
s_11 = sosdecvar('md1',zmd);
s_21 = sosdecvar('md3',zmd);



     case 2
      
% Short Period


% s_2,m_2,s_3 (i.e. s_1,s_2,s_3 in the paper)
zmDeg=2;
zm = monomials(x, 0:zmDeg );
s_1 = sosdecvar('m0',zm);

s_2 = sosdecvar('m1',zm);

s_3 = sosdecvar('m3',zm);

% s_0 (i.e. s_0 in the paper)
pDeg=4; % 
zp = monomials(x,0:pDeg);
[s_0,~] = polydecvar('p',zp); %
         


%%%%% Create \delta-multipliers variables (they will feature in the
%%%%% polynomial \Gamma Eq. 18d)

zmdDeg=1;
zmd = monomials(x_delta, 0:zmdDeg );

s_01 = sosdecvar('md2',zmd);
s_11 = sosdecvar('md1',zmd);
s_21 = sosdecvar('md3',zmd);


end



gopts = gsosoptions; % default definition of gsosopt
% Definition of some of the options (see documentatio of sosopt)
gopts.maxobj=0;  
gopts.minobj=-80;  
gopts.checkfeas='both'; 
gopts.solver = 'sedumi'; 
gopts.display = 'on';



sopts = sosoptions; % default definition of sosopt
% Definition of some of the options (see documentatio of sosopt)
sopts.checkfeas='both'; 
sopts.solver = 'sedumi'; % 'sedumi'
 
 

iter_max = 50;
n_iter=0;



% % % % Algorithm selection
% 1 -> Extension of Algorithm 1 to the robust case 
% 2 -> Extension of Algorithm 2 to the robust case; this is Algorithm 6  (IS_R-3 Steps) in the paper
% 3 -> Extension of Algorithm 3 to the robust case 

Alg_ID = 3;



 
 % initializing the cell where level set functions and size are written at each step of the iteration
c0 = cell(iter_max,1);

switch Alg_ID
     
     case 1 
         
unfeas_A2=1;
unfeas_A1=0;
iter_h= struct('R_h_A1',c0,'V_h_A1',c0,'gamma_A1',c0);

     case 2
         
unfeas_A2=0;
unfeas_A1=1;
iter_h= struct('R_h_A2',c0,'V_h_A2',c0,'gamma_A2',c0);
V_0=Vlin;
     
    case 3
         
unfeas_A2=0;
unfeas_A1=0;

iter_h= struct('R_h_A1',c0,'V_h_A1',c0,'gamma_A1',c0,'R_h_A2',c0,'V_h_A2',c0,'gamma_A2',c0);

 

end


 
% rERA Algorithm

 while n_iter <= iter_max
     
 n_iter = n_iter+1
  
 
   if Alg_ID == 1 || Alg_ID == 3 
       
%        Algorithm 1 (IS-2 Steps)
       
  [unfeas_A1,iter_h] = JFI_rERA_Alg1(V_n,R_n,s_3,s_1,s_2,s_0,R_0,x,f,gopts,sopts,iter_h,n_iter,x_delta,s_11,s_01,s_21,m_c);
 

 
if unfeas_A1 
    
    if unfeas_A2

    fprintf('\n Iteration stopped for unfeasibility \n');
 
    break
    
    end
    
else
    

%     This is done to speed up the next bisection
 gopts.minobj=-iter_h(n_iter).gamma_A1*2;
 
%  Initialization of R and V for the next iteration (note, V_0 needed here
%  only if the hybrid approach is used, because  Algorithm 1 does not need
%  an initialization for V)

 R_0 = iter_h(n_iter).R_h_A1;
 V_0=iter_h(n_iter).V_h_A1;
     
 gamma=iter_h(n_iter).gamma_A1; % this has only plotting purposes
end
 
 
 
 
     
   end
   
   
    
if Alg_ID == 2 || Alg_ID == 3 


   [unfeas_A2,iter_h] = JFI_rERA_Alg2(V_0,s_1,s_2,s_0,R_n,R_0,s_3,x,f,gopts,sopts,V_n,iter_h,n_iter,s_11,s_01,s_21,m_c,x_delta);




  
 
 if unfeas_A2 
    
    if unfeas_A1

    fprintf('\n Iteration stopped for unfeasibility \n');
 
    break
    
    end
    
else
    
 %     This is done to speed up the next bisection

 gopts.minobj=-iter_h(n_iter).gamma_A2*2;
 
 %  Initialization of R and V for the next iteration 
 
 R_0 = iter_h(n_iter).R_h_A2;
 V_0=iter_h(n_iter).V_h_A2;
 
 
 gamma=iter_h(n_iter).gamma_A2; % this has only plotting purposes
   
 end

 
 

  end

 
 
 end
 
 
% % %  Example of Plotting (see documentation of multipoly for utilities to
% % % plot the results)

switch system
     
     case 1 
         
  domain = [-3 3 -3 3];

 figure(1)                                                   
[C,ph(1)]=pcontour(R_0,gamma,domain,'r');
hold on


     case 2
         
 domain = 9*[-3 3 -3 3];

R_0_draw = subs ( R_0 ,[x(3:5)], [0; 0; 0] ); % projection onto z1-z2

figure(1)                                                   
[C,ph(1)]=pcontour(R_0_draw,gamma,domain,'g');
         
end

