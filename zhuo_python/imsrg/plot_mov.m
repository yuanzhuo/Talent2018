clc;
clear;
N = 30;
data=zeros(N,256,256);

F = ['data/srg_conv_' , num2str(0) , '.txt'];
D_t = importdata(F);
size(D_t)
x = 1:size(D_t);
y = x;
D_t = importdata(F);
contour(x,y,D_t)

M = moviein(N);
for i = 1:N
    F = ['data/srg_conv_' , num2str(i) , '.txt'];
    D_t = importdata(F);
    size(D_t)
    x = 1:size(D_t);
    y = x;
    
    contour(x,y,D_t)
    %set(gcf, 'doublebuffer', 'on')
    %shading interp
    
    %colormap(hot)
    M(i) = getframe(gcf);
end
size(data)
movie(M,2);




for i = 1:1
    x = 1:size(data(i));
    y = x;
    
    
    %surf(x,y,z)
   %M(i) = getframe;
end

%movie(M,2);
%data

