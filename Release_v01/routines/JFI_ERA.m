%------------------------------------------------------------------
% Inner estimates of the (nominal) region of attraction (ERA) with Invariant Sets
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
         
pvar x1 x2;
x = [x1;x2];

x1dot = -x2;
x2dot = x1+(x1^2-1)*x2;

f = [x1dot; x2dot];




     case 2
% Short Period
         
pvar x1 x2 x3 x4 x5;
x=[x1;x2;x3;x4;x5]; 
% [z1,z2,eta1,eta2,z3]; % paper nomenclature
% z1= pitch rate; z2 = angle of attack; z3 = pitch angle; \eta1,2= controller states
%



f1 =  3.0000*x(4)+1.3489*x(3)-.56250*x(5)-1.3489*x(2)-3.0000*x(1)+.82272e-1*x(1)*x(2)+.43981*x(2)^2-.36992e-1*x(2)*x(3)...
    -.82272e-1*x(2)*x(4)+.15426e-1*x(2)*x(5)+.22335*x(2)^3;
f2 = .91136*x(1)-.64450*x(2)-.54444e-1*x(2)^2-.16621e-1*x(5)+.10889*x(2)*x(5)-.54444e-1*x(5)^2+.39857e-1*x(3)...
    +.88643e-1*x(4);
f3 = -.60401e-1*x(1)-.16621e-1*x(5)-.60464*x(3)+.88643e-1*x(4);
f4 = -.75000*x(1)-.28125*x(5);
f5 = x(1);
f = [f1;f2;f3;f4;f5];


 


end

n_x = size(x,1);


% Lyapunov Function of the linearized equilibrium
Q = eye(n_x);

[Vlin,~,~]=linstab(f,x,Q);


% Initialization for R
R_0 = Vlin;


%%%%% Create level set functions variables (with multipoly utilities)

% V_n
VnDeg=4;
zVn = monomials(x,1:VnDeg);
[V_n,~] = polydecvar('Vn',zVn); % generic definition of poly

% R_n
RnDeg = VnDeg;
zRn = monomials(x,1:RnDeg);
[R_n,~] = polydecvar('Rn',zRn); % generic definition of poly


%%%%% Create multipliers variables

% s_1,s_2,s_3 in the paper
zmDeg=1;
zm = monomials(x, 0:zmDeg );
s_1 = sosdecvar('m0',zm);

s_2 = sosdecvar('m1',zm);

s_3 = sosdecvar('m3',zm);

% s_0 in the paper
pDeg=2; % 
zp = monomials(x,0:pDeg);
[s_0,~] = polydecvar('p',zp); %

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
% 1 -> Algorithm 1 (IS-2 Steps): Step A1-1 & Step A1-2
% 2 -> Algorithm 2 (IS-3 Steps): Step A2-1 & Step A2-2 & Step A2-3
% 3 -> Algorithm 3 (IS-Hyb): Stage 1, Stage 2

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


% ERA Algorithm

 while n_iter <= iter_max
     
 n_iter = n_iter+1
  
 
   if Alg_ID == 1 || Alg_ID == 3 
       
%        Algorithm 1 (IS-2 Steps)
       
  [unfeas_A1,iter_h] = JFI_ERA_Alg1(V_n,R_n,s_3,s_1,s_2,s_0,R_0,x,f,gopts,sopts,iter_h,n_iter);
 

 
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


   [unfeas_A2,iter_h] = JFI_ERA_Alg2(V_0,s_1,s_2,s_0,R_n,R_0,s_3,x,f,gopts,sopts,V_n,iter_h,n_iter);




  
 
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
[C,ph(1)]=pcontour(R_0,gamma,domain,'g');
hold on


     case 2
         
 domain = 9*[-3 3 -3 3];

R_0_draw = subs ( R_0 ,[x(3:5)], [0; 0; 0] ); % projection onto z1-z2

figure(1)                                                   
[C,ph(1)]=pcontour(R_0_draw,gamma,domain,'g');
         
end

