x_i = [0 0; 0 1; 1 0; 1 1]
test = [0 1 1 0]';
i = 2;
j = 4;
k = 1;
a = 0.5;
b = -0.5;
max = 50;
n = 2;
h = 4;

v_ij = (b-a).*rand(i,j) + a;
beta = 0.7 * ((h)^(1/n));
v_ij = beta * (v_ij/norm(v_ij));
v_ij = v_ij';
w_jk = (b-a).*rand(j,k) + a;
beta_min = -1.4;
beta_max = 1.4;
bias1 = (beta_max - beta_min) .*rand(1) + beta_min;
bias2 = (beta_max - beta_min) .*rand(1) + beta_min;
v_oj = bias1;
w_ok = bias2;

for iterations = 1:max
    
    Z_inj = sum(x_i*v_ij')' + v_oj;
    z_j = tansig(Z_inj,[1,1]);
    Y_ink = w_jk.*z_j+w_ok;
    y_k = tansig(Y_ink,[1,1]);
    error = test - y_k;
    alpha = 0.9;
    
    f_dash_Y_ink = 0.5 * ((1+y_k) .* (1-y_k));
    delta_k = error.*f_dash_Y_ink;
    delta_w_jk = alpha * delta_k .* z_j;
    delta_wok = alpha * delta_k;
    
    f_dash_z_inj = 0.5 * ((1+z_j)' * (1-z_j));
    delta_j = sum(delta_k .* w_jk) * f_dash_z_inj;
    delta_v_ij = alpha * delta_j*x_i;
    delta_v_oj = alpha * delta_j;
    
    w_jk = w_jk + delta_w_jk;
    w_ok = w_ok + delta_wok;
    v_ij = v_ij + delta_v_ij;
    v_oj = v_oj + delta_v_oj;
    
    s(iterations) = 0.5 * (sum(error))^2;
end
figure;plot(0:49,s(1:50))