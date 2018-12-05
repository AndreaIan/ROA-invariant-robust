function [infoM,gamma_2,R_n_M] = JFI_ERA_StepA2_3(V_n,R_n,m_3,m_0_f,m_1_f,p_f,R_0,t12,x,f,sopts)

L2 =  0*1e-8*(x'*x);

% Definition of SOS constraints and objective function (-\gamma, because the solver implements a minimization)

pvar t2
sosconstr2 = polyconstr;

sosconstr2(1) = V_n - L2 - m_0_f * (-t2-R_n) >=0;
sosconstr2(2) = - jacobian(V_n,x)*f - m_1_f * (-t2-R_n) >=0;
sosconstr2(3) = - jacobian(R_n,x)*f - p_f * (-t2-R_n) >=0;
sosconstr2(4) = m_3 >=0;
sosconstr2(5) =  (-t2-R_n) -m_3*(-t12-R_0) >=0;



% Minimization (via bisection). This time I use sosopt instead of gsosopt
% because there are not bilinear terms, hence there is no need for
% bisection (see documentation on difference between gsosopt/sosopt for
%     further info)
    
    
    
    objM = t2 ;
    [infoM,doptM,sossolM] = sosopt(sosconstr2,x,objM,sopts);
    
    if infoM.feas~=0
        fprintf('\n Step A2-3 OK - optimal gamma found \n');
        
        t2_M = infoM.obj;
        gamma_2 = -t2_M;
        R_n_M = subs(R_n,doptM);
        
        
    else
        fprintf('\n Step A2-3 optimization failed - Let us try with just feasibility \n');
        objM = [] ;
        
        [infoM,doptM,sossolM] = sosopt(sosconstr2,x,objM,sopts);
        
        if infoM.feas~=0
            fprintf(' \n Step A2-3 OK - feasible solution found \n');
            
            t2_M = subs(t2,doptM);
            gamma_2 = -t2_M.coefficient;
            R_n_M = subs(R_n,doptM);
            
        else
            
            fprintf('\n Step A2-3 unfeasible \n');
            
            gamma_2=[];
            R_n_M=[];
        end

end
