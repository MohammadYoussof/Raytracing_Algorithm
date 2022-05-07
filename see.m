vid=videoinput('winvideo',1);                                            
figure(3);preview(vid);                    
set(vid,'ReturnedColorspace','rgb')
img = getsnapshot(vid);
figure (2); imagesc(img);
red = img(:,:,1); % Canal Rojo
green = img(:,:,2); % Canal Verde
blue = img(:,:,3); % Canal Azul
a = zeros(size(img, 1), size(img, 2));
just_red = cat(3, red, a, a);
just_green = cat(3, a, green, a);
just_blue = cat(3, a, a, blue);
I= just_red - just_green; 
I2=im2bw(I, graythresh(I));  
I3= bwareaopen(I2, 100);  
se=strel('disk',6);
I4=imclose(I3,se); 
Ic=imfill(I4,'holes');