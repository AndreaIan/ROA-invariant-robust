function [unfeas,iter] = JFI_ERA_Alg2(V_0,m_0,m_1,p,R_n,R_0,m_3,x,f,gopts,sopts,V_n,iter,n_iter)

unfeas = 0;
% keyboard
[info_s1,gamma1_LB,m_0_g,m_1_g,p_g] = JFI_ERA_StepA2_1(V_0,m_0,m_1,p,R_0,x,f,gopts);

if ~(info_s1.feas==1)
    
    
    unfeas=1;
    
    iter(n_iter).R_h_A2=[];
iter(n_iter).gamma_A2=[];
iter(n_iter).V_h_A2=[];


    return
    
    
end

[info_s2,gamma2,V_n_M] = JFI_ERA_StepA2_2(V_n,R_0,m_0_g,m_1_g,p_g,x,f,sopts,gamma1_LB);

if ~(info_s2.feas==1)
    
    unfeas=1;
    
      iter(n_iter).R_h_A2=[];
iter(n_iter).gamma_A2=[];
iter(n_iter).V_h_A2=[];

    return
    
    
end

[info_s3,gamma3,R_n_M] = JFI_ERA_StepA2_3(V_n_M,R_n,m_3,m_0_g,m_1_g,p_g,R_0,-gamma2,x,f,sopts);

if ~(info_s3.feas==1)

unfeas=1;


end

iter(n_iter).R_h_A2=R_n_M;
iter(n_iter).gamma_A2=gamma3;
iter(n_iter).V_h_A2=V_n_M;



end
