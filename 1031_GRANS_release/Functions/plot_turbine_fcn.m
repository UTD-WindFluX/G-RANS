function plot_turbine_fcn(x,y,z,angle,diameter)
    load turbine 
    angle=angle+180;
    part1=[-turbine.hub{1,1}';turbine.hub{1,2}';zeros(1,length(turbine.hub{1,2}))];
    part2=[-turbine.hub{2,1}';turbine.hub{2,2}';zeros(1,length(turbine.hub{2,2}))];
    part3=[-turbine.down{1,1}';turbine.down{1,2}';zeros(1,length(turbine.down{1,2}))];
    part4=[-turbine.up{1,1}';turbine.up{1,2}';zeros(1,length(turbine.up{1,2}))];
    M_rot=[cos(angle*pi/180) -sin(angle*pi/180) 0;sin(angle*pi/180) cos(angle*pi/180) 0; 0 0 1];
    part1_rot=M_rot*part1*diameter+repmat([x,y,z]',1,length(part1(1,:)));
    part2_rot=M_rot*part2*diameter+repmat([x,y,z]',1,length(part2(1,:)));
    part3_rot=M_rot*part3*diameter+repmat([x,y,z]',1,length(part3(1,:)));
    part4_rot=M_rot*part4*diameter+repmat([x,y,z]',1,length(part4(1,:)));
    hold on
    patch('XData',part1_rot(1,:),'YData',part1_rot(2,:),'ZData',part1_rot(3,:),'FaceColor', [.7 .7 .7])
    hold on
    patch('XData',part2_rot(1,:),'YData',part2_rot(2,:),'ZData',part2_rot(3,:),'FaceColor', [.7 .7 .7])
    hold on
    patch('XData',part3_rot(1,:),'YData',part3_rot(2,:),'ZData',part3_rot(3,:),'FaceColor', [.7 .7 .7]);
    hold on
    patch('XData',part4_rot(1,:),'YData',part4_rot(2,:),'ZData',part4_rot(3,:),'FaceColor', [.7 .7 .7]);
end