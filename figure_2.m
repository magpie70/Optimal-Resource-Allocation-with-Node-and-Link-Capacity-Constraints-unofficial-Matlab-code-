clc
clear all
%%% Initialization
N = 200;
L = 600;
S = 1000;
M = 5000;
%%% BA-network model
ba_network = scalefree(N,10,L/N);
G1 = graph(ba_network,'upper');
shortest_path = distances(G1);


i = randi([1,N],1);
j = randi([1,N],1);
coord_S = [[i j]];
len = 1;
while len<S
    i = randi([1,N],1);
    j = randi([1,N],1);
    if ~any(ismember(coord_S,[i j],'rows')) && ~any(ismember(coord_S,[j i],'rows')) && shortest_path(i,j)>0
        coord_S = [coord_S; [i j]];
        len=len+1;
    end
end

ds = [];
for i=1:S
    ds = [ds shortest_path(coord_S(i,1),coord_S(i,2))];
end
ds = ds+1;

V_s = {};
for i=1:S
    V_s=[V_s; (shortestpath(G1,coord_S(i,1),coord_S(i,2)))];
end

F_n = cell(N,1);
for i=1:S
    temp = V_s{i};
    for j=1:size(temp,2)
       F_n{temp(j)} = [F_n{temp(j)} i];
    end
end

F_l = cell(N,N);
for i=1:S
    temp = V_s{i};
    for j=1:size(temp,2)-1
        F_l{temp(j),temp(j+1)} = [F_l{temp(j),temp(j+1)} i];
        F_l{temp(j+1),temp(j)} = [F_l{temp(j+1),temp(j)} i];
    end
end

%%% BA
BA_utility_alpha = [];
T1 = ba_network;
p = (1./ds.^2)';

for alpha = [linspace(0.01,0.1,10),linspace(0.1,1,10)]
    alpha
    cvx_begin sdp quiet
        variable x(S)
        variable D(N)
        maximize sum(-p.*inv_pos(x))
        subject to

    x >= 0;
    
    for i=1:N
        sum1 = 0;
        for j=F_n{i}
            sum1 = sum1 + x(j);
        end
        sum1<=D(i);
    end

    for m=1:N
        for n=1:N
            sum2 = 0;
            for j=F_l{m,n}
                sum2 = sum2 + x(j);
            end
            sum2 <= T1(m,n)*alpha*(D(m)+D(n));
        end
    end

    sum(D) == M ;
    D>=0;
    
    cvx_end
    
    res = -1./((ds.^2).*x');
    
    BA_utility_alpha = [BA_utility_alpha sum(res)];
    
end

%%% ER network model
er_network = randomGraph(N,0.5,L);
G2 = graph(er_network,'upper');
shortest_path = distances(G2);

i = randi([1,N],1);
j = randi([1,N],1);
coord_S = [[i j]];
len = 1;
while len<S
    i = randi([1,N],1);
    j = randi([1,N],1);
    if ~any(ismember(coord_S,[i j],'rows')) && ~any(ismember(coord_S,[j i],'rows')) && shortest_path(i,j)>0
        coord_S = [coord_S; [i j]];
        len=len+1;
    end
end

ds = [];
for i=1:S
    ds = [ds shortest_path(coord_S(i,1),coord_S(i,2))];
end
ds = ds+1;

V_s = {};
for i=1:S
    V_s=[V_s; (shortestpath(G2,coord_S(i,1),coord_S(i,2)))];
end

F_n = cell(N,1);
for i=1:S
    temp = V_s{i};
    for j=1:size(temp,2)
       F_n{temp(j)} = [F_n{temp(j)} i];
    end
end

F_l = cell(N,N);
for i=1:S
    temp = V_s{i};
    for j=1:size(temp,2)-1
        F_l{temp(j),temp(j+1)} = [F_l{temp(j),temp(j+1)} i];
        F_l{temp(j+1),temp(j)} = [F_l{temp(j+1),temp(j)} i];
    end
end
%%%
ER_utility_alpha = [];
T2 = er_network;
p = (1./ds.^2)';
for alpha = [linspace(0.01,0.1,10),linspace(0.1,1,10)]
    alpha
    cvx_begin sdp quiet
        variable x(S)
        variable D(N)
        maximize sum(-p.*inv_pos(x))
        subject to

    x >= 0;
    
    for i=1:N
        sum1 = 0;
        for j=F_n{i}
            sum1 = sum1 + x(j);
        end
        sum1<=D(i);
    end

    for m=1:N
        for n=1:N
            sum2 = 0;
            for j=F_l{m,n}
                sum2 = sum2 + x(j);
            end
            sum2 <= T2(m,n)*alpha*(D(m)+D(n));
        end
    end

    sum(D) == M ;
    D>=0;
    
    cvx_end
    
    res = -1./((ds.^2).*x');
    
    ER_utility_alpha = [ER_utility_alpha sum(res)];
    
end


%%%
values_alpha = [linspace(0.01,0.1,10),linspace(0.1,1,10)];
values_M = 0:1000:10000;
figure(1)

plot(values_alpha, BA_utility_alpha, 'b-^', 'LineWidth',1,'MarkerSize',10);
hold on;
plot(values_alpha, ER_utility_alpha, 'r--*', 'LineWidth',1,'MarkerSize',8);

set(gca,'xtick',M);
xlim([values_alpha(1), values_alpha(end)]);
xlabel('Parameter \alpha'),ylabel('Utility'); 
grid on;
handle=legend('BA network, optimal resource allocation','ER network, optimal resource allocation');
set(gca,'Xscale','log');
set(handle, 'Interpreter', 'LaTex');
set(handle, 'FontSize', 13, 'Location', 'best');
set(gca, 'FontSize', 12);