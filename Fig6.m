function [ F6 ] = Fig6
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

F6 =figure(3); 
set( F6, 'Position', [150 10 540 310], 'Name', 'RMS');
xlabel('Iteraciones')
ylabel('RMS','FontSize', 14,'Rotation',90,'interpreter','latex')
end

