function [info,g_LB,m_0_f,m_1_f,p_f] = JFI_ERA_StepA2_1(V_n_0,m_0,m_1,p,R_0,x,f,gopts,V_n)

L2 =  0*1e-8*(x'*x);

% Definition of SOS constraints and objective function (-\gamma, because the solver implements a minimization)

pvar t
sosconstr = polyconstr;
sosconstr(1) = m_0>=0;
sosconstr(2) = m_1>=0;


sosconstr(3) = - jacobian(R_0,x)*f - p * (-t-R_0) >=0;
sosconstr(4) = V_n_0 - L2 - m_0 * (-t-R_0) >=0;
sosconstr(5) = - jacobian(V_n_0,x)*f - m_1 * (-t-R_0) >=0;

% Minimization (via bisection)

[info,dopt] = gsosopt(sosconstr,x,t,gopts);

if ~isempty(info.tbnds)
    fprintf('\n Step A2-1 OK - optimal gamma found with bisection \n');
    
    
    g_LB = -info.tbnds(2);
    m_0_f = subs(m_0,dopt);
    m_1_f = subs(m_1,dopt);
    p_f = subs(p,dopt);
    
else
    
    
    fprintf('\n Step A2-1 unfeasible \n');
    g_LB=[];
    m_0_f=[];
    m_1_f=[];
    p_f=[];
    
end







end