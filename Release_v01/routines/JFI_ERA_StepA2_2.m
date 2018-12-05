function [infoM,gamma_2,V_n_M] = JFI_ERA_StepA2_2(V_n,R_n,m_0_f,m_1_f,p_f,x,f,sopts,gamma1_LB)

L2 =  0*1e-8*(x'*x);

% Definition of SOS constraints and objective function (-\gamma, because the solver implements a minimization)


pvar t2
sosconstr2 = polyconstr;

sosconstr2(1) = V_n - L2 - m_0_f * (-t2-R_n) >=0;
sosconstr2(2) = - jacobian(V_n,x)*f - m_1_f * (-t2-R_n) >=0;
sosconstr2(3) = -t2-gamma1_LB >=0;
sosconstr2(4) = - jacobian(R_n,x)*f - p_f * (-t2-R_n) >=0;




% Minimization (via bisection). This time I use sosopt instead of gsosopt
% because there are not bilinear terms, hence there is no need for
% bisection (see documentation on difference between gsosopt/sosopt for
% further info)



objM = t2 ;
[infoM,doptM,sossolM] = sosopt(sosconstr2,x,objM,sopts);


if infoM.feas~=0
    fprintf('\n Step A2-2 OK - optimal gamma found \n');
    
    t2_M = infoM.obj;
    gamma_2 = -t2_M;
    V_n_M = subs(V_n,doptM);
    
    
else
    fprintf('\n Step A2-2 optimization failed - Let us try with just feasibility \n');
    
    objM = [] ;
    
    [infoM,doptM,sossolM] = sosopt(sosconstr2,x,objM,sopts);
    
    if infoM.feas~=0
        
        fprintf(' \n Step A2-2 OK - feasible solution found \n');
        t2_M = subs(t2,doptM);
        gamma_2 = -t2_M.coefficient;
        V_n_M = subs(V_n,doptM);
        
    else
        fprintf('\n Step A2-2 unfeasible \n');
        
        gamma_2 = [];
        V_n_M = [];
        
    end

end

 
 
 