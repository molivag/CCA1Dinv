function [ F6 ] = Fig6( i, RMS )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

F6 =figure(6); 
set( F6, 'Position', [200 10 600 300], 'Name', 'RMS');
bar(i,RMS)
xlabel('Iteraciones')
ylabel('RMS','FontSize', 14,'Rotation',90,'interpreter','latex')
end

