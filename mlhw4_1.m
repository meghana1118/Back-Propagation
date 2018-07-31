 %xi = [0 0; 0 1; 1 0; 1 1]
 %test = [0 1 1 0]';
xi = [-1 -1 -1;-1 -1 1;-1 1 -1;-1 1 1;1 -1 -1;1 -1 1;1 1 -1;1 1 1]
test = [0 1 1 0 1 0 0 1]';
i = 3;
j = 8;
k = 1;
a = 0.5;
b = -0.5;
voj= 0.4;
wok = 0.3;
max = 50;
vij = (b-a).*rand(i,j) + a;
vinitial = vij;
vij = vij';
wjk = (b-a).*rand(j,k) + a;
winitial = wjk;

for iterations = 1:max
        Zinj = sum(xi*vij')' + voj;
        zj = sigmf(Zinj,[1,1]);
        Yink = wjk.*zj+wok;
        yk = sigmf(Yink,[1,1]);
        error = test - yk;
        alpha = 0.9;
        fdashYink = yk.*(1-yk);
        deltak = error.*fdashYink;
        deltawjk = alpha * deltak .* zj;
        deltawok = alpha * deltak;
        
        wjk = wjk + deltawjk;
        wok = wok + deltawok;
        
        fdashzinj = zj' * (1-zj);
        deltaj = sum(deltak .* wjk) * fdashzinj;
        deltavij = alpha * deltaj*xi;
        deltavoj = alpha * deltaj;
        
        vij = vij + deltavij;
        voj = voj + deltavoj;
        
        s(iterations) = 0.5 * (sum(error))^2;
end

figure;plot(0:49,s(1:50))
yk