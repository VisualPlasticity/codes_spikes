figure;
x = [-5000:10:5000-10];
y = [-5000:10:5000-10];
imagesc(x,y,ones(800,800));

% ML; AP
Center = [4557,-2973];
Radius = 500;
points = [4673,-3230;4693,-2940;4474,-2681;4390,-3189];

%draw circle
 drawcircle('Center',Center,'Radius',Radius,'FaceAlpha',0.1,'Color',[1 0 0],'LineWidth',2)
hold on
cols = lines(length(points));
legendname={};
for pid=1:length(points)
    h(pid)=plot(points(pid,1),points(pid,2),'*','color',cols(pid,:));
    legendname={legendname{:} ['Insertion ' num2str(pid)]};
    Distance(pid) = sqrt((points(pid,1)-Center(1))^2+(points(pid,2)-Center(2))^2);
end
legend([h(:)],legendname)
set(gca,'XDir','normal','YDir','normal')
%Zoom in on craniotomy
xlim([Center(1)-Radius*1.1 Center(1)+Radius*1.1])
ylim([Center(2)-Radius*1.1 Center(2)+Radius*1.1])

ylabel('AnteriorPosterior')
xlabel('MedialLateral')
box off

Distance