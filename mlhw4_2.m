x_i = [0 0; 0 1; 1 0; 1 1]
test = [0 1 1 0]';
i = 2;
j = 4;
k = 1;
a = 0.5;
b = -0.5;
v_oj= 0.4;
w_ok = 0.3;
max = 50;

v_ij = [0.2710 , -0.4133 ; 0.3476, -0.3258 ; -0.0383,-0.4961 ;0.4218,0.0573]

w_jk=[-0.4597;0.1596;-0.0853; 0.2762]


for iterations = 1:max
        Z_inj = sum(x_i*v_ij')' + v_oj;
        z_j = sigmf(Z_inj,[1,1]);
        Y_ink = w_jk.*z_j+w_ok;
        y_k = sigmf(Y_ink,[1,1]);
        error = test - y_k;
        alpha = 0.9;
        
        f_dash_Y_ink = y_k.*(1-y_k);
        delta_k = error.*f_dash_Y_ink;
        delta_w_jk = alpha * delta_k .* z_j;
        delta_wok = alpha * delta_k;
        
       f_dash_z_inj = z_j' * (1-z_j);
       delta_j = sum(delta_k .*w_jk) * f_dash_z_inj;
       delta_v_ij = alpha * delta_j*x_i;
       delta_v_oj = alpha * delta_j;
       
       mu = 0.5;
       w_jk = w_jk + delta_w_jk + mu * (w_jk - (w_jk -delta_w_jk));
       w_ok = w_ok + delta_wok + mu * (w_ok - (w_ok - delta_wok));
       v_ij = v_ij + delta_v_oj + mu * (v_oj - (v_oj - delta_v_oj));
       v_oj = v_oj + delta_v_oj + mu * (v_oj - (v_oj - delta_v_oj));
       
       s(iterations)= 0.5 * (sum(error))^2;
end

figure;plot(0:49,s(1:50))