function [info,g_LB,m_0_f,m_1_f,p_f,m_d1_f,m_d2_f,m_d3_f] = JFI_rERA_StepA1_1(V_n,m_0,m_1,p,R_0,x,f,gopts,m_d1,m_d2,m_d3,N_d,x_delta)

    L2 =  0*1e-8*(x'*x);
% Definition of SOS constraints and objective function (-\gamma, because the solver implements a minimization)

pvar t                                                                 
sosconstr = polyconstr;
sosconstr(1) = m_0>=0; 
sosconstr(2) = m_1>=0; 
sosconstr(3) = m_d1>=0; 
sosconstr(4) = m_d2>=0; 
sosconstr(5) = m_d3>=0;

sosconstr(6) = - jacobian(R_0,x)*f - p * (-t-R_0) - m_d1*N_d>=0;
sosconstr(7) = V_n - L2 - m_0 * (-t-R_0)- m_d2*N_d>=0;
sosconstr(8) = - jacobian(V_n,x)*f - m_1 * (-t-R_0) - m_d3*N_d>=0;

% Minimization (via bisection)
[info,dopt] = gsosopt(sosconstr,x_delta,t,gopts); 


if ~isempty(info.tbnds)
        fprintf('\n Step rA1-1 OK - optimal gamma found with bisection \n');

g_LB = -info.tbnds(2); 
m_0_f = subs(m_0,dopt);
m_1_f = subs(m_1,dopt);
p_f = subs(p,dopt); 

m_d1_f = subs(m_d1,dopt);
m_d2_f = subs(m_d2,dopt);
m_d3_f = subs(m_d3,dopt);


else
    
    fprintf('\n Step rA1-1 unfeasible \n');
m_0_f = [];
m_1_f = [];
p_f = [];
 g_LB = []; % 
m_d1_f =[];
m_d2_f = [];
m_d3_f = [];

end



end