% MATLAB program to convert video into slow motion
clc;clear;close all;
  
 % load the video.
obj = VideoReader('/Users/youssof/Downloads/seismic_modelling_folder/Tina_seismic_modeling/testmovienew.mp4');  
  
% Write in new variable
obj2= VideoWriter('xyz.mp4');    
  
% decrease framerate 
obj2.FrameRate = 5;              
open(obj2);
  
% for reading frames one by one
while hasFrame(obj)              
    k = readFrame(obj); 
  
    % write the frames in obj2.         
    obj2.writeVideo(k);          
end
  
close(obj2);