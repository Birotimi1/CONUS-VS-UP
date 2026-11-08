function plot_layer_model(matrix, color, width, label)
% Step profile including an extended half-space for plotting.
edges=[0;cumsum(matrix(1:end-1,1));max(100,sum(matrix(:,1))+max(10,.1*sum(matrix(:,1))))];
x=reshape([matrix(:,3) matrix(:,3)]',[],1);
z=reshape([edges(1:end-1) edges(2:end)]',[],1);
plot(x,z,'Color',color,'LineWidth',width,'DisplayName',label);
end
