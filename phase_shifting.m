%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Zhengzhong Huang, 2024
% The version of Matlab for this code is R2020a
% Reference:Quantitative phase imaging based on holography: Trends and new perspectives
% 《Light: Science & Applications》
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
%% Amplitude & phase of object

load('amplitude_object','amplitude_object')% Amplitude of object
load('phase_object','phase_object')% Phase of object

sample = amplitude_object.*exp(1i.*phase_object);

%% Using small region (Non-essential)

image_size=500;
starn_m=640;
starn_n=200;
sample1=sample(starn_m+1:starn_m+image_size,starn_n+1:starn_n+image_size);
figure(2),imshow(abs(sample1),[])

%% Defined CTF

pupil_radius=round(image_size/8)-1;% Radius of CTF
[x_array,y_array]=meshgrid(1:image_size,1:image_size);
x_array=x_array-floor(max(x_array(:))/2+1); % center of image to be zero
y_array=y_array-floor(max(y_array(:))/2+1); % center of image to be zero
pupil_function=(x_array./pupil_radius).^2+(y_array./pupil_radius).^2 <= 1;% CTF
figure(3),imshow(abs(pupil_function),[])

%% Defined band-limit object wave

Sample1=ifftshift(fft2(fftshift(sample1)));
object=ifftshift(ifft2(fftshift(Sample1.*pupil_function)));% Band-limit object wave

figure(4),imshow(angle(object),[])

%% Defined reference wave: 4 step phase shifting

ref_field_Fourier=zeros(image_size,image_size);
ref_spatial_freq=0;% Carried spatial frequency
ref_field_Fourier(floor(image_size/2+1),floor(image_size/2+1)-ref_spatial_freq) = ...
max(abs(object(:)))*length(object(:))*2;% The spatial frequency of reference wave
ref_field = fftshift(ifft2(ifftshift(ref_field_Fourier)));% Reference wave
ref_field1=ref_field.*exp(1i.*0);% Reference wave with phase shifting step 1
ref_field2=ref_field.*exp(1i.*pi./2);% Reference wave with phase shifting step 2
ref_field3=ref_field.*exp(1i.*pi);% Reference wave with phase shifting step 3
ref_field4=ref_field.*exp(1i.*3.*pi./2);% Reference wave with phase shifting step 4

figure(5),imshow(angle(ref_field),[])

%% Hologram

hologram1=abs(object+ref_field1).^2;% Hologram
hologram2=abs(object+ref_field2).^2;% Hologram
hologram3=abs(object+ref_field3).^2;% Hologram
hologram4=abs(object+ref_field4).^2;% Hologram
Hologram=ifftshift(fft2(fftshift(hologram1)));% Spatial frequency of hologram

figure(6),imshow((abs(hologram4)),[])

%% Reconstruction by phase shifting

reconstruction = (hologram4-hologram2+1i.*(hologram1-hologram3))./(4.*ref_field);
figure(7),imshow((angle(reconstruction)),[])


