
function plot_my_network(A,thres)

A=A-diag(diag(A));
A(abs(A)<thres)=0;

g=digraph(A');

plot(g,'layout','layered','LineWidth',g.Edges.Weight*20);

%plot(digraph(A'),'layout','layered');

