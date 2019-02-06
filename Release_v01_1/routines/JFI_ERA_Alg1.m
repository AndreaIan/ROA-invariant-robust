function [unfeas,iter] = JFI_ERA_Alg1(V_n,R_n,m_3,m_0,m_1,p,R_0,x,f,gopts,sopts,iter,n_iter)

unfeas = 0;

% : Step A1-1 
[info_s1,gamma1_LB,m_0_g,m_1_g,p_g] = JFI_ERA_StepA1_1(V_n,m_0,m_1,p,R_0,x,f,gopts);


if ~(info_s1.feas==1)

unfeas=1;

iter(n_iter).R_h_A1=[];
iter(n_iter).gamma_A1=[];
iter(n_iter).V_h_A1=[];


return

    
end    



% & Step A1-2

[info_s2,gamma2,V_n_M,R_n_M,~] = JFI_ERA_StepA1_2(V_n,R_n,m_3,m_0_g,m_1_g,p_g,R_0,-gamma1_LB,x,f,sopts);


if ~(info_s2.feas==1)

unfeas=1;
    
end 




iter(n_iter).R_h_A1=R_n_M;
iter(n_iter).gamma_A1=gamma2;
iter(n_iter).V_h_A1=V_n_M;


   
 




end
  
  
  
  
  
  
  
  
  
  
