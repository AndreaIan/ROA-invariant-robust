function [info,gamma_2,V_n_M,m_d1_fN,m_d2_fN,m_d3_fN] = ERA_V_Step_R(V_n,R_n,m_0_f,m_1_f,p_f,x,f,sopts,gamma1_LB,m_d1,m_d2,m_d3,m_d1_f,m_d2_f,m_d3_f,N_d,x_delta)

  L2 =  0*1e-8*(x'*x);

  % Definition of SOS constraints and objective function (-\gamma, because the solver implements a minimization)

  pvar t2                                                                 
sosconstr2 = polyconstr;

%%% Note that the polynomials m_d# associated with the uncertainties can be
%%% in principle kept fixed to the value at previous steps (as done for the
%%% other multipliers). Therefore, first the problem is solved by keeping
%%% them fixed and, in case of unfeasibility, they will be employed as
%%% optimization variables (note that this will not determine bilinearities
%%% in the problem). This strategy is devised to limit the computational
%%% cost, but the code can be easily modified to optimize directly over the
%%% m_d# polynomials.

 sosconstr2(1) = V_n - L2 - m_0_f * (-t2-R_n) - m_d2_f*N_d>=0;
 sosconstr2(2) = - jacobian(V_n,x)*f - m_1_f * (-t2-R_n) - m_d3_f*N_d>=0;
   sosconstr2(3) = -t2-gamma1_LB >=0;
    sosconstr2(4) = - jacobian(R_n,x)*f - p_f * (-t2-R_n) - m_d1_f*N_d>=0;





objM = t2 ; 
[info1,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts); 

if info1.feas~=0
    fprintf('\n Step rA2-2 OK - optimal gamma found \n');

t2_M = info1.obj;
gamma_2 = -t2_M;
V_n_M = subs(V_n,doptM);

m_d1_fN = [];
m_d2_fN = [];
m_d3_fN = [];
info=info1;
else
    fprintf('\n Step rA2-2 optimization failed - Let us try with just feasibility \n');
    objM = [] ; 

    [info2,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts);
    
    if info2.feas~=0
        fprintf(' \n Step rA2-2 OK - feasible solution found \n');

 t2_M = subs(t2,doptM);
gamma_2 = -t2_M.coefficient;
V_n_M = subs(V_n,doptM);
m_d1_fN = [];
m_d2_fN = [];
m_d3_fN = [];
info=info2;

    else
        
        fprintf(' \n Step rA2-2 unfeasible - m_di polynomials added as decision variables \n');
  sosconstr2(1) = V_n - L2 - m_0_f * (-t2-R_n) - m_d2*N_d>=0;
 sosconstr2(2) = - jacobian(V_n,x)*f - m_1_f * (-t2-R_n) - m_d3*N_d>=0;
    sosconstr2(4) = - jacobian(R_n,x)*f - p_f * (-t2-R_n) - m_d1*N_d>=0;

     sosconstr2(5) = m_d1>=0; 
sosconstr2(6) = m_d2>=0; 
sosconstr2(7) = m_d3>=0; 
    
   objM = t2 ; 
[info3,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts); 



if info3.feas~=0
        fprintf(' \n Step rA2-2 OK - optimal gamma found  \n');

t2_M = info3.obj;
gamma_2 = -t2_M;
V_n_M = subs(V_n,doptM);
m_d1_fN = subs(m_d1,doptM);
m_d2_fN = subs(m_d2,doptM);
m_d3_fN = subs(m_d3,doptM);
info=info3;

else
                fprintf('\n Step rA2-2 optimization failed - Let us try with just feasibility \n');

    objM = [] ; 

    [info4,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts); 
    info=info4;

    if info4.feas~=0
        
        fprintf(' \n Step rA2-2 OK - feasible solution found \n');

 t2_M = subs(t2,doptM);
gamma_2 = -t2_M.coefficient;
V_n_M = subs(V_n,doptM);
m_d1_fN = subs(m_d1,doptM);
m_d2_fN = subs(m_d2,doptM);
m_d3_fN = subs(m_d3,doptM);
    else

                fprintf('\n Step rA1-2 unfeasible \n');
  
    gamma_2 = [];
V_n_M = [];
m_d1_fN = [];
m_d2_fN = [];
m_d3_fN = [];
    end
    
end

    end
end


 
 
 