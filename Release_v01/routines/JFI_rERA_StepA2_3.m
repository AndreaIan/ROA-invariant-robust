function [infoM,gamma_2,R_n_M] = ERA_Rg_Step_R(V_n,R_n,m_0_f,m_1_f,p_f,x,f,sopts,t12,m_d1,m_d2,m_d3,m_d1_f,m_d2_f,m_d3_f,N_d,x_delta,m_3,R_0)
                                               

L2 =  0*1e-8*(x'*x);

pvar t2                                                                 
sosconstr2 = polyconstr;


sosconstr2(1) = V_n - L2 - m_0_f * (-t2-R_n) - m_d2_f*N_d>=0;
 sosconstr2(2) = - jacobian(V_n,x)*f - m_1_f * (-t2-R_n) - m_d3_f*N_d>=0;
 sosconstr2(3) = - jacobian(R_n,x)*f - p_f * (-t2-R_n) - m_d1_f*N_d>=0;
 
 sosconstr2(4) = m_3 >=0; 
 sosconstr2(5) =  (-t2-R_n) -m_3*(-t12-R_0) >=0; 


  
  


objM = t2 ; 
[infoM,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts); 


if infoM.feas~=0
    fprintf('\n Step rA2-3 OK - optimal gamma found \n');


    t2_M = infoM.obj;
gamma_2 = -t2_M;
R_n_M = subs(R_n,doptM);



else
    fprintf('\n Step rA2-3 optimization failed - Let us try with just feasibility \n');
    objM = [] ; 

    [infoM,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts); 
    
    if infoM.feas~=0
        fprintf(' \n Step rA2-3 OK - feasible solution found \n');

 t2_M = subs(t2,doptM);
gamma_2 = -t2_M.coefficient;
R_n_M = subs(R_n,doptM);

    else
        fprintf(' \n Step rA2-3 unfeasible - m_di polynomials added as decision variables \n');
  sosconstr2(1) = V_n - L2 - m_0_f * (-t2-R_n) - m_d2*N_d>=0;
 sosconstr2(2) = - jacobian(V_n,x)*f - m_1_f * (-t2-R_n) - m_d3*N_d>=0;
    sosconstr2(3) = - jacobian(R_n,x)*f - p_f * (-t2-R_n) - m_d1*N_d>=0;

sosconstr2(6) = m_d1>=0; 
sosconstr2(7) = m_d2>=0;
sosconstr2(8) = m_d3>=0; 
    
   objM = t2 ;
[infoM,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts); 


if infoM.feas~=0
        fprintf(' \n Step rA2-3 OK - optimal gamma found  \n');

t2_M = infoM.obj;
gamma_2 = -t2_M;
R_n_M = subs(R_n,doptM);


else
                fprintf('\n Step rA2-3 optimization failed - Let us try with just feasibility \n');
    objM = [] ; 

    [infoM,doptM,sossolM] = sosopt(sosconstr2,x_delta,objM,sopts);
    
    if infoM.feas~=0
        fprintf(' \n Step rA2-3 OK - feasible solution found \n');

 t2_M = subs(t2,doptM);
gamma_2 = -t2_M.coefficient;
R_n_M = subs(R_n,doptM);

    else

                fprintf('\n Step rA2-3 unfeasible \n');
  
   gamma_2 = [];
R_n_M = [];
    end
    
end

    end
end


