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

%% Using small regions of two object wave (Non-essential)

image_size=500;
starn_m1=640;
starn_n1=100;
starn_m2=640;
starn_n2=600;

sample1=sample(starn_m1+1:starn_m1+image_size,starn_n1+1:starn_n1+image_size);
sample2=sample(starn_m2+1:starn_m2+image_size,starn_n2+1:starn_n2+image_size);
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
Sample2=ifftshift(fft2(fftshift(sample2)));
object1=ifftshift(ifft2(fftshift(Sample1.*pupil_function)));% Band-limit object wave
object2=ifftshift(ifft2(fftshift(Sample2.*pupil_function)));% Band-limit object wave

Object1=ifftshift(ifft2(fftshift(object1)));
Object1=circshift(Object1,[-round(3*image_size/8),-round(3*image_size/8)]);% Create tilt object wave 1
object1=ifftshift(fft2(fftshift(Object1)));
Object2=ifftshift(ifft2(fftshift(object2)));
Object2=circshift(Object2,[round(3*image_size/8),-round(3*image_size/8)]);% Create tilt object wave 2
object2=ifftshift(fft2(fftshift(Object2)));

figure(4),imshow(log(abs(Object1)),[])

%% Defined reference wave

ref_field_Fourier=zeros(image_size,image_size);
ref_spatial_freq=0;% Carried spatial frequency
ref_field_Fourier(floor(image_size/2+1),floor(image_size/2+1)-ref_spatial_freq) = ...
max(abs(object1(:)))*length(object1(:))*2;% The spatial frequency of reference wave
ref_field = fftshift(ifft2(ifftshift(ref_field_Fourier)));% Reference wave

figure(5),imshow(angle(ref_field),[])

%% Hologram

hologram=abs(object1+object2+ref_field).^2;% Multipleced hologram
Hologram=ifftshift(fft2(fftshift(hologram)));% Spatial frequency of hologram

figure(6),imshow(log(abs(Hologram)),[])

%% Reconstruction of off-axis hologram

Filter_radius=round(image_size/8)-1;% Radius of filtered function
Filter_function1=(x_array./Filter_radius).^2+(y_array./Filter_radius).^2 <= 1;
Filter_function1=circshift(Filter_function1,[round(3*image_size/8),round(3*image_size/8)]);% Defined filtered function 1
Filter_function2=(x_array./Filter_radius).^2+(y_array./Filter_radius).^2 <= 1;
Filter_function2=circshift(Filter_function2,[-round(3*image_size/8),round(3*image_size/8)]);% Defined filtered function 2

Reconstruction1=Hologram.*Filter_function1;
Reconstruction1=circshift(Reconstruction1,[-round(3*image_size/8),-round(3*image_size/8)]);% Reconstructed spatial frequency of object wave 1
reconstruction1 = fftshift(ifft2(ifftshift(Reconstruction1)));% Reconstructed object wave 1

Reconstruction2=Hologram.*Filter_function2;
Reconstruction2=circshift(Reconstruction2,[round(3*image_size/8),-round(3*image_size/8)]);% Reconstructed spatial frequency of object wave 2
reconstruction2 = fftshift(ifft2(ifftshift(Reconstruction2)));% Reconstructed object wave 2

figure(7),imshow((angle(reconstruction2)),[])



