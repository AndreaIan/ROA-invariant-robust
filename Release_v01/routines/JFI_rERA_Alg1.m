function [unfeas,iter] = JFI_rERA_Alg1(V_n,R_n,m_3,m_0,m_1,p,R_0,x,f,gopts,sopts,iter,n_iter,x_delta,m_d1,m_d2,m_d3,N_d)

unfeas = 0;

[info_s1,gamma1_LB,m_0_g,m_1_g,p_g,m_d1_f,m_d2_f,m_d3_f] = JFI_rERA_StepA1_1(V_n,m_0,m_1,p,R_0,x,f,gopts,m_d1,m_d2,m_d3,N_d,x_delta);


if ~(info_s1.feas==1)

unfeas=1;

iter(n_iter).R_h_A1=[];
iter(n_iter).gamma_A1=[];
iter(n_iter).V_h_A1=[];


return

    
end    





% Step A1-2

[info_s2,gamma2,V_n_M,R_n_M,~] = JFI_rERA_StepA1_2(V_n,R_n,m_3,m_0_g,m_1_g,p_g,R_0,-gamma1_LB,x,f,sopts,m_d1,m_d2,m_d3,N_d,x_delta,m_d1_f,m_d2_f,m_d3_f); 



if ~(info_s2.feas==1)

unfeas=1;
    
end 




iter(n_iter).R_h_A1=R_n_M;
iter(n_iter).gamma_A1=gamma2;
iter(n_iter).V_h_A1=V_n_M;


   
 




end
  
  
  
  
  
  
  
  
  
  
