function [unfeas,iter] = JFI_rERA_Alg2(V_0,m_0,m_1,p,R_n,R_0,m_3,x,f,gopts,sopts,V_n,iter,n_iter,m_d1,m_d2,m_d3,N_d,x_delta)
   
unfeas = 0;


[info_s1,gamma1_LB,m_0_g,m_1_g,p_g,m_d1_f,m_d2_f,m_d3_f] = JFI_rERA_StepA2_1(V_0,m_0,m_1,p,R_0,x,f,gopts,m_d1,m_d2,m_d3,N_d,x_delta);

 if ~(info_s1.feas==1)

unfeas=1;

  iter(n_iter).R_h_A2=[];
iter(n_iter).gamma_A2=[];
iter(n_iter).V_h_A2=[];

return

    
end    


        [info_s2,gamma2,V_n_M,m_d1_fv,m_d2_fv,m_d3_fv] = JFI_rERA_StepA2_2(V_n,R_0,m_0_g,m_1_g,p_g,x,f,sopts,gamma1_LB,m_d1,m_d2,m_d3,m_d1_f,m_d2_f,m_d3_f,N_d,x_delta);

        
  if ~(info_s2.feas==1)

unfeas=1;

  iter(n_iter).R_h_A2=[];
iter(n_iter).gamma_A2=[];
iter(n_iter).V_h_A2=[];


return

    
if ~isempty(m_d1_fv)
    m_d1_f=m_d1_fv;
    m_d2_f=m_d2_fv;
    m_d3_f=m_d3_fv;
end


end    


    [info_s3,gamma3,R_n_M] = JFI_rERA_StepA2_3(V_n_M,R_n,m_0_g,m_1_g,p_g,x,f,sopts,-gamma2,m_d1,m_d2,m_d3,m_d1_f,m_d2_f,m_d3_f,N_d,x_delta,m_3,R_0);

if ~(info_s3.feas==1)

unfeas=1;


end

iter(n_iter).R_h_A2=R_n_M;
iter(n_iter).gamma_A2=gamma3;
iter(n_iter).V_h_A2=V_n_M;



end
