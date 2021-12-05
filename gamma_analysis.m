function varargout = gamma_analysis(varargin)
% GAMMA_ANALYSIS MATLAB code for gamma_analysis.fig
%      GAMMA_ANALYSIS, by itself, creates a new GAMMA_ANALYSIS or raises the existing
%      singleton*.
%
%      H = GAMMA_ANALYSIS returns the handle to a new GAMMA_ANALYSIS or the handle to
%      the existing singleton*.
%
%      GAMMA_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAMMA_ANALYSIS.M with the given input arguments.
%
%      GAMMA_ANALYSIS('Property','Value',...) creates a new GAMMA_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gamma_analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gamma_analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gamma_analysis

% Last Modified by GUIDE v2.5 01-Aug-2018 14:27:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gamma_analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @gamma_analysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gamma_analysis is made visible.
function gamma_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gamma_analysis (see VARARGIN)

% Choose default command line output for gamma_analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gamma_analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gamma_analysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

clear all
% --- Executes on button press in input_measure.
function input_measure_Callback(hObject, eventdata, handles)
% hObject    handle to input_measure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in input_simulation.
global flag_pdd;    global flag_ocr;
global z_measure; global dose_measure;
global dose_measure_5;    global dose_measure_10; global dose_measure_20;

if flag_pdd==1
[FileName,PathName] = uigetfile('D:\DONG NAI\*.txt');
set(handles.measure,'string',FileName);

formatSpec='%f';
fid=fopen([PathName,FileName],'r');
store=fscanf(fid,formatSpec,inf);
fclose(fid);

store_z_measure=[];
chay=1;
for i=1:2:length(store)
    store_z_measure(chay,1)=store(i,1);
    chay=chay+1;
end
chay=1;
store_dose_measure=[];
for i=2:2:length(store)
    store_dose_measure(chay,1)=store(i,1);
    chay=chay+1;
end
z_measure=store_z_measure;  dose_measure=store_dose_measure;
%display('flag_pdd==1');
end



if flag_ocr==1
[FileName,PathName] = uigetfile('D:\DONG NAI\*.txt');
set(handles.measure,'string',FileName);

formatSpec='%f';
fid=fopen([PathName,FileName],'r');
store=fscanf(fid,formatSpec,inf);
fclose(fid);

                                %5%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
store_z_measure=[];
chay=1;
for i=1:6:length(store)
    store_z_measure(chay,1)=store(i,1);
    chay=chay+1;
end

chay=1;
store_dose_measure_5=[];
for i=2:6:length(store)
    store_dose_measure_5(chay,1)=store(i,1);
    chay=chay+1;
end
z_measure=store_z_measure;  
dose_measure_5=100.*store_dose_measure_5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %10%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chay=1;
store_dose_measure_10=[];
for i=4:6:length(store)
    store_dose_measure_10(chay,1)=store(i,1);
    chay=chay+1;
end
dose_measure_10=100.*store_dose_measure_10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %20%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chay=1;
store_dose_measure_20=[];
for i=6:6:length(store)
    store_dose_measure_20(chay,1)=store(i,1);
    chay=chay+1;
end
dose_measure_20=100.*store_dose_measure_20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%axes(handles.axes1); %title('PDD ban dau');
%a=plot(z_measure,dose_measure,'.k');   hold on
%title('PDD');
%legend(a,'Dose measure');   grid on
%xlabel('Depth'); ylabel('Dose');

function input_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to input_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flag_pdd;    global flag_ocr;
global z_simulation; global dose_simulation;
global dose_simulation_5;  global dose_simulation_10;   global dose_simulation_20;

if flag_pdd==1 && flag_ocr==0
%[FileName,PathName] = uigetfile('D:\DONG NAI\*.3ddose');
[FileName,PathName] = uigetfile({'*.jpg';'*.dcm';'*.jpeg';'*.*';'*.ima'},'File Selector');

set(handles.simul,'string',FileName);

display(FileName);
%save_name = ([PathName,FileName]);
%save_name = strrep(save_name, '.3ddose', '_Results');

fid = fopen ([PathName,FileName], 'r');

nx = cell2mat(textscan(fid, '%f', 1));
ny = cell2mat(textscan(fid, '%f', 1));
nz = cell2mat(textscan(fid, '%f', 1));

nx_2 = round(nx/2);
ny_2 = round(ny/2);
nz_2 = round(nz/2);

x_dim = nx+1;
y_dim = ny+1;
z_dim = nz+1;
no_of_voxels = nx*ny*nz;

x_cor_bndry = cell2mat(textscan(fid, '%f', x_dim));
y_cor_bndry = cell2mat(textscan(fid, '%f', y_dim));
z_cor_bndry = cell2mat(textscan(fid, '%f', z_dim));

Dose_Val = cell2mat(textscan(fid, '%f', no_of_voxels));
Error_Val = cell2mat(textscan(fid, '%f', no_of_voxels));

end_1 = ftell(fid);
eofstat_1 = feof(fid);
rest_of_file = cell2mat(textscan(fid, '%f'));
end_2 = ftell(fid);
eofstat_2 = feof(fid);
fclose(fid);

Dose_3D_Mat = reshape(Dose_Val, nx,ny,nz);
Error_3D_Mat = reshape(Error_Val, nx,ny,nz);

clear Dose_Val Error_Val

PDD=1:nz;
PDD=PDD*0;
Dose_3D_Mat(Dose_3D_Mat<0)=0;

CAX_D_Max = max(max(max(Dose_3D_Mat)));
temp_3D = 100 .* Dose_3D_Mat ./ CAX_D_Max;
PDD(:)=temp_3D(nx_2,ny_2,:);

A(:,:) = temp_3D (nx_2,:, :);
B(:,:) = temp_3D (:,ny_2, :);
C(:,:) = temp_3D (:, :, 1);
A_e(:,:) = Error_3D_Mat (nx_2,:, :);
B_e(:,:) = Error_3D_Mat (:,ny_2, :);

x_simulation(1:nx)=0;
for m=1:nx;
    x_simulation(m)= (x_cor_bndry(m)+x_cor_bndry(m+1))/2;
end

y_simulation(1:ny)=0;
for m=1:ny;
    y_simulation(m)= (y_cor_bndry(m)+y_cor_bndry(m+1))/2;
end

y_simulation=transpose(y_simulation);
x_simulation=transpose(x_simulation);

z_simulation(1:nz)=0;
for m=1:nz;
    z_simulation(m)= (z_cor_bndry(m)+z_cor_bndry(m+1))/2;
end
dose_simulation = 100 .* PDD ./ max(PDD);
z_simulation=transpose(z_simulation);   dose_simulation=transpose(dose_simulation);
%display('flag_pdd==1');
end


if flag_ocr==1 && flag_pdd==0
    [FileName,PathName] = uigetfile('D:\DONG NAI\*.3ddose');
set(handles.simul,'string',FileName);

%save_name = ([PathName,FileName]);
%save_name = strrep(save_name, '.3ddose', '_Results');

fid = fopen ([PathName,FileName], 'r');

nx = cell2mat(textscan(fid, '%f', 1));
ny = cell2mat(textscan(fid, '%f', 1));
nz = cell2mat(textscan(fid, '%f', 1));

nx_2 = round(nx/2);
ny_2 = round(ny/2);
nz_2 = round(nz/2);

x_dim = nx+1;
y_dim = ny+1;
z_dim = nz+1;
no_of_voxels = nx*ny*nz;

x_cor_bndry = cell2mat(textscan(fid, '%f', x_dim));
y_cor_bndry = cell2mat(textscan(fid, '%f', y_dim));
z_cor_bndry = cell2mat(textscan(fid, '%f', z_dim));

Dose_Val = cell2mat(textscan(fid, '%f', no_of_voxels));
Error_Val = cell2mat(textscan(fid, '%f', no_of_voxels));

end_1 = ftell(fid);
eofstat_1 = feof(fid);
rest_of_file = cell2mat(textscan(fid, '%f'));
end_2 = ftell(fid);
eofstat_2 = feof(fid);
fclose(fid);

Dose_3D_Mat = reshape(Dose_Val, nx,ny,nz);
Error_3D_Mat = reshape(Error_Val, nx,ny,nz);

clear Dose_Val Error_Val

PDD=1:nz;
PDD=PDD*0;
Dose_3D_Mat(Dose_3D_Mat<0)=0;
%display(size(Dose_3D_Mat));

CAX_D_Max = max(max(max(Dose_3D_Mat)));
temp_3D = 100 .* Dose_3D_Mat ./ CAX_D_Max;
PDD(:)=temp_3D(nx_2,ny_2,:);

A(:,:) = temp_3D (nx_2,:, :);
B(:,:) = temp_3D (:,ny_2, :);
C(:,:) = temp_3D (:, :, 1);
A_e(:,:) = Error_3D_Mat (nx_2,:, :);
B_e(:,:) = Error_3D_Mat (:,ny_2, :);


x_simulation(1:nx)=0;
for m=1:nx;
    x_simulation(m)= (x_cor_bndry(m)+x_cor_bndry(m+1))/2;
end
x_simulation=transpose(x_simulation);

y_simulation(1:ny)=0;
for m=1:ny;
    y_simulation(m)= (y_cor_bndry(m)+y_cor_bndry(m+1))/2;
end
y_simulation=transpose(y_simulation);

%z_simulation(1:nz)=0;
%for m=1:nz;
%    z_simulation(m)= (z_cor_bndry(m)+z_cor_bndry(m+1))/2;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OCR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dose_simulation = 100 .* PDD ./ max(PDD);
dose=100.*Dose_3D_Mat ./ CAX_D_Max;
%ocr=[];
%display('size(dose)');    display(size(dose));
%ocr=dose(nx_2,ny_2,:);
dose_simulation_5 = dose(2,:,16);
dose_simulation_5=transpose(dose_simulation_5);

dose_simulation_10 = dose(2,:,26);
dose_simulation_10=transpose(dose_simulation_10);

dose_simulation_20 = dose(2,:,46);
dose_simulation_20=transpose(dose_simulation_20);
%x_simulation=transpose(x_simulation);
%y_simulation=transpose(y_simulation);
z_simulation=y_simulation;
    %display('flag_ocr==1');
    

%display(size(dose_measure)); 
end
%%%%%%%% KET THUC PHAN DOC FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%% SAVE FILE TINH TOAN %%%%%%%%%%%%%%%
%save(save_name, 'x_cor', 'y_cor', 'z_cor', 'temp_PDD', 'temp_3D', 'Error_3D_Mat');
%savename=strrep(FileName,'.3ddose','-dose-simulation.txt');
%fid=fopen(savename,'w');
%fprintf(fid,'%f\r\n',dose_simulation);
%fclose(fid);

%savename1=strrep(FileName,'.3ddose','-z-simulation.txt');
%fid=fopen(savename1,'w');
%fprintf(fid,'%f\r\n',z_simulation);
%fclose(fid);


% --- Executes on button press in fit.
function fit_Callback(hObject, eventdata, handles)
% hObject    handle to fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global z_measure; global dose_measure;
global dose_measure_5;    global dose_measure_10;  global dose_measure_20;
global z_simulation; global dose_simulation;
global dose_simulation_5;  global dose_simulation_10;   global dose_simulation_20;

global fit_first; global fit_step;  global fit_last;
global flag_pdd;    global flag_ocr;

global z_meas; global dose_meas;
global z_meas_5;    global z_meas_10;    global z_meas_20;
global dose_meas_5;  global dose_meas_10; global dose_meas_20;
global z_simul; global dose_simul;
global z_simul_5;   global z_simul_10;   global z_simul_20;
global dose_simul_5;  global dose_simul_10;  global dose_simul_20;

first_fit=fit_first;  last_fit=fit_last;

if flag_pdd==1 

%axes(handles.axes1); %title('PDD ban dau');
%plot(z_measure,dose_measure,'.k');   hold on
%plot(z_simulation,dose_simulation,'.b');    title('PDD');
%legend('Dose measure','Dose simulation');   grid on
%xlabel('Depth'); ylabel('Dose');

%axes(handles.axes1); %title('PDD ban dau');
%plot(z_measure,dose_measure,'.k');   hold on
%plot(z_simulation,dose_simulation,'.b');    title('PDD');
%legend('Dose measure','Dose simulation');   grid on
%xlabel('Depth'); ylabel('Dose');

if first_fit < min(z_measure)
    first_fit=min(z_measure);
end
if last_fit > max(z_measure)
    last_fit=max(z_measure);
end

[z_measure_new,dose_measure_new]=fit_data(z_measure,dose_measure,first_fit,fit_step,last_fit);
  
%%%%%%%%% CAP NHAT GIA TRI CHO Z_MEASURE BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
z_measure=z_measure_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% CAP NHAT GIA TRI CHO DOSE_MEASURE BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
dose_measure=dose_measure_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

first_fit=fit_first;    last_fit=fit_last;
if first_fit < min(z_simulation)
    first_fit=min(z_simulation);
end
if last_fit > max(z_simulation)
    last_fit=max(z_simulation);
end

[z_simulation_new,dose_simulation_new]=fit_data(z_simulation,dose_simulation,first_fit,fit_step,last_fit);

%%%%%%%%% CAP NHAT GIA TRI CHO Z_simulation BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
z_simulation=z_simulation_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% CAP NHAT GIA TRI CHO DOSE_simulation BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
dose_simulation=dose_simulation_new;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display(max(store_dose_simulation));
%axes(handles.axes1);
%plot(z_measure,dose_measure,'k'); hold on
%plot(z_simulation,dose_simulation,'b');
%title('PDD AFTER FIT');   %legend('Dose measure', 'Dose simulation'); grid on
%xlabel('Depth');    ylabel('Dose');
%display('flag_pdd==1');

axes(handles.axes1); %title('PDD ban dau');
plot(z_measure,dose_measure,'.k');   hold on
plot(z_simulation,dose_simulation,'.b');    title('PDD');
legend('Dose measure','Dose simulation');   grid on
xlabel('Depth'); ylabel('Dose');

%%%% dua cac gia tri nay vao cau lenh tinh chi so gamma %%%%%%%%%
z_meas=z_measure;    
z_simul=z_simulation; 
dose_meas=dose_measure;
dose_simul=dose_simulation;
end




if flag_ocr==1 && flag_pdd==0 

%%%%%%%%%%%%% loai nhanh 0 cho dose_measure %%%%%%%%%%%%%%%%%%%%%%%
chay=1;
z_measure_5=z_measure;  z_measure_10=z_measure;  z_measure_20=z_measure;


store_z_measure_5=[];  store_z_measure_10=[];    store_z_measure_20=[];
store_dose_measure_5=[];    store_dose_measure_10=[];    store_dose_measure_20=[];
for i=1:floor(length(z_measure_5)/2)
    if dose_measure_5(i,1)>min(dose_simulation_5(1:floor(length(dose_simulation_5)/2)))%0.2
        %dose_measure(i,1)=0;
        store_dose_measure_5(chay,1)=dose_measure_5(i,1);
        store_z_measure_5(chay,1)=z_measure_5(i,1);
        chay=chay+1;
    end
end
for i=floor(length(z_measure_5)/2)+1:length(z_measure_5)
    if dose_measure_5(i,1)>min(dose_simulation_5(floor(length(dose_simulation_5)/2)+1:length(dose_simulation_5),1))%0.2
        %dose_measure(i,1)=0;
        store_dose_measure_5(chay,1)=dose_measure_5(i,1);
        store_z_measure_5(chay,1)=z_measure_5(i,1);
        chay=chay+1;
    end
end
chay=1;

for i=1:floor(length(z_measure_10)/2)
    if dose_measure_10(i,1)>min(dose_simulation_10(1:floor(length(dose_simulation_10)/2)))%0.2
        %dose_measure(i,1)=0;
        store_dose_measure_10(chay,1)=dose_measure_10(i,1);
        store_z_measure_10(chay,1)=z_measure_10(i,1);
        chay=chay+1;
    end
end
for i=floor(length(z_measure_10)/2)+1:length(z_measure_10)
    if dose_measure_10(i,1)>min(dose_simulation_10(floor(length(dose_simulation_10)/2)+1:length(dose_simulation_10),1))%0.2
        %dose_measure(i,1)=0;
        store_dose_measure_10(chay,1)=dose_measure_10(i,1);
        store_z_measure_10(chay,1)=z_measure_10(i,1);
        chay=chay+1;
    end
end

chay=1;
for i=1:floor(length(z_measure_20)/2)
    if dose_measure_20(i,1)>min(dose_simulation_20(1:floor(length(dose_simulation_20)/2)))%0.2
        %dose_measure(i,1)=0;
        store_dose_measure_20(chay,1)=dose_measure_20(i,1);
        store_z_measure_20(chay,1)=z_measure_20(i,1);
        chay=chay+1;
    end
end
for i=floor(length(z_measure_20)/2)+1:length(z_measure_20)
    if dose_measure_20(i,1)>min(dose_simulation_20(floor(length(dose_simulation_20)/2)+1:length(dose_simulation_20),1))%0.2
        %dose_measure(i,1)=0;
        store_dose_measure_20(chay,1)=dose_measure_20(i,1);
        store_z_measure_20(chay,1)=z_measure_20(i,1);
        chay=chay+1;
    end
end

dose_measure_5=store_dose_measure_5;    dose_measure_10=store_dose_measure_10;
dose_measure_20=store_dose_measure_20;
z_measure_5=store_z_measure_5;  z_measure_10=store_z_measure_10;  z_measure_20=store_z_measure_20;


%display(size(dose_measure));
%%%%%%%%%%%%%%%%%%%%% FIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_fit=fit_first;    last_fit=fit_last;
if first_fit < min(z_measure_5)
    first_fit=-7;%min(z_measure);
end
if last_fit > max(z_measure_5)
    last_fit=7;%max(z_measure);
end
[z_measure_new_5,dose_measure_new_5]=fit_data(z_measure_5,dose_measure_5,first_fit,fit_step,last_fit);

first_fit=fit_first;    last_fit=fit_last;
if first_fit < min(z_measure_10)
    first_fit=-7;%min(z_measure);
end
if last_fit > max(z_measure_10)
    last_fit=7;%max(z_measure);
end
[z_measure_new_10,dose_measure_new_10]=fit_data(z_measure_10,dose_measure_10,first_fit,fit_step,last_fit);

first_fit=fit_first;    last_fit=fit_last;
if first_fit < min(z_measure_20)
    first_fit=-7;%min(z_measure);
end
if last_fit > max(z_measure_20)
    last_fit=7;%max(z_measure);
end
[z_measure_new_20,dose_measure_new_20]=fit_data(z_measure_20,dose_measure_20,first_fit,fit_step,last_fit);
  
%%%%%%%%% CAP NHAT GIA TRI CHO Z_MEASURE BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
z_measure_5=z_measure_new_5;    z_measure_10=z_measure_new_10;    z_measure_20=z_measure_new_20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% CAP NHAT GIA TRI CHO DOSE_MEASURE BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
dose_measure_5=dose_measure_new_5;  
dose_measure_10=dose_measure_new_10;
dose_measure_20=dose_measure_new_20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% FIT simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% loai nhanh 0 cho dose_simulation 5, 10, 20 %%%%%%%%%%%%%%%%%%%%%%%
z_simulation_5=z_simulation;    z_simulation_10=z_simulation;    z_simulation_20=z_simulation;

chay=1;
store_z_simulation_5=[];  store_z_simulation_10=[];  store_z_simulation_20=[];
store_dose_simulation_5=[]; store_dose_simulation_10=[]; store_dose_simulation_20=[];

for i=1:floor(length(z_simulation_5)/2)
    if dose_simulation_5(i,1)>min(dose_simulation_5(1:floor(length(dose_simulation_5)/2),1))%0.2
        %dose_simulation(i,1)=0;
        store_dose_simulation_5(chay,1)=dose_simulation_5(i,1);
        store_z_simulation_5(chay,1)=z_simulation(i,1);
        chay=chay+1;
    end
end
for i=floor(length(z_simulation_5)/2)+1:length(z_simulation_5)
    if dose_simulation_5(i,1)>min(dose_simulation_5(floor(length(dose_simulation_5)/2)+1:length(dose_simulation_5),1))%0.2
        %dose_simulation(i,1)=0;
        store_dose_simulation_5(chay,1)=dose_simulation_5(i,1);
        store_z_simulation_5(chay,1)=z_simulation(i,1);
        chay=chay+1;
    end
end
chay=1;

for i=1:floor(length(z_simulation_10)/2)
    if dose_simulation_10(i,1)>min(dose_simulation_10(1:floor(length(dose_simulation_10)/2),1))%0.2
        %dose_simulation(i,1)=0;
        store_dose_simulation_10(chay,1)=dose_simulation_10(i,1);
        store_z_simulation_10(chay,1)=z_simulation(i,1);
        chay=chay+1;
    end
end
for i=floor(length(z_simulation_10)/2)+1:length(z_simulation_10)
    if dose_simulation_10(i,1)>min(dose_simulation_10(floor(length(dose_simulation_10)/2)+1:length(dose_simulation_10),1))%0.2
        %dose_simulation(i,1)=0;
        store_dose_simulation_10(chay,1)=dose_simulation_10(i,1);
        store_z_simulation_10(chay,1)=z_simulation(i,1);
        chay=chay+1;
    end
end
chay=1;

for i=1:floor(length(z_simulation_20)/2)
    if dose_simulation_20(i,1)>min(dose_simulation_20(1:floor(length(dose_simulation_20)/2),1))%0.2
        %dose_simulation(i,1)=0;
        store_dose_simulation_20(chay,1)=dose_simulation_20(i,1);
        store_z_simulation_20(chay,1)=z_simulation(i,1);
        chay=chay+1;
    end
end
for i=floor(length(z_simulation_20)/2)+1:length(z_simulation_20)
    if dose_simulation_20(i,1)>min(dose_simulation_20(floor(length(dose_simulation_20)/2)+1:length(dose_simulation_20),1))%0.2
        %dose_simulation(i,1)=0;
        store_dose_simulation_20(chay,1)=dose_simulation_20(i,1);
        store_z_simulation_20(chay,1)=z_simulation(i,1);
        chay=chay+1;
    end
end

z_simulation_5=store_z_simulation_5;    z_simulation_10=store_z_simulation_10;    
z_simulation_20=store_z_simulation_20;
dose_simulation_5=store_dose_simulation_5;  dose_simulation_10=store_dose_simulation_10;
dose_simulation_20=store_dose_simulation_20;

first_fit=fit_first;    last_fit=fit_last;
if first_fit <= min(z_simulation_5)
    first_fit=-7;%min(z_simulation);
end
if last_fit >= max(z_simulation_5)
    last_fit=7;%max(z_simulation);
end
[z_simulation_new_5,dose_simulation_new_5]=fit_data(z_simulation_5,dose_simulation_5,first_fit,fit_step,last_fit);

first_fit=fit_first;    last_fit=fit_last;
if first_fit <= min(z_simulation_10)
    first_fit=-7;%min(z_simulation);
end
if last_fit >= max(z_simulation_10)
    last_fit=7;%max(z_simulation);
end
[z_simulation_new_10,dose_simulation_new_10]=fit_data(z_simulation_10,dose_simulation_10,first_fit,fit_step,last_fit);

first_fit=fit_first;    last_fit=fit_last;
if first_fit <= min(z_simulation_20)
    first_fit=-7;%min(z_simulation);
end
if last_fit >= max(z_simulation_20)
    last_fit=7;%max(z_simulation);
end
[z_simulation_new_20,dose_simulation_new_20]=fit_data(z_simulation_20,dose_simulation_20,first_fit,fit_step,last_fit);

%%%%%%%%% CAP NHAT GIA TRI CHO Z_simulation BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
z_simulation_5=z_simulation_new_5;  z_simulation_10=z_simulation_new_10;
z_simulation_20=z_simulation_new_20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% CAP NHAT GIA TRI CHO DOSE_simulation BAN DAU THEO CAC GIA TRI FIT %%%%%%%%%%%%%%%%%
dose_simulation_5=dose_simulation_new_5;
dose_simulation_10=dose_simulation_new_10;
dose_simulation_20=dose_simulation_new_20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% dua cac gia tri nay vao cau lenh tinh chi so gamma %%%%%%%%%
z_meas=z_measure;
z_meas_5=z_measure_5;   z_meas_10=z_measure_10;   z_meas_20=z_measure_20;
z_simul_5=z_simulation_5;  z_simul_10=z_simulation_10; z_simul_20=z_simulation_20; 
dose_meas_5=dose_measure_5; dose_meas_10=dose_measure_10; dose_meas_20=dose_measure_20;
dose_simul_5=dose_simulation_5; dose_simul_10=dose_simulation_10; dose_simul_20=dose_simulation_20;
display('DONE!');
beep
set(handles.pop,'visible','on');    set(handles.text27,'visible','on');
end


%clear all % xoa cac gia tri da luu trong bo nho de lan sau thuc hien khong bi sai
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fit_value_first_Callback(hObject, eventdata, handles)
% hObject    handle to fit_value_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_value_first as text
%        str2double(get(hObject,'String')) returns contents of fit_value_first as a double


% --- Executes during object creation, after setting all properties.

global fit_first;
fit_first=get(hObject,'string');  fit_first=str2num(fit_first);

function fit_value_first_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_value_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function step_Callback(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step as text
%        str2double(get(hObject,'String')) returns contents of step as a double


% --- Executes during object creation, after setting all properties.

global fit_step;
fit_step=get(hObject,'string');  fit_step=str2num(fit_step);

function step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fit_value_last_Callback(hObject, eventdata, handles)
% hObject    handle to fit_value_last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fit_value_last as text
%        str2double(get(hObject,'String')) returns contents of fit_value_last as a double


% --- Executes during object creation, after setting all properties.

global fit_last;
fit_last=get(hObject,'string');  fit_last=str2num(fit_last);

function fit_value_last_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_value_last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in input_dose_measure.

function gamma_index_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global dta; global dose_diff;
global z_meas; global dose_meas;
global z_meas_5;    global z_meas_10;    global z_meas_20;
global dose_meas_5;  global dose_meas_10; global dose_meas_20;
global z_simul; global dose_simul;
global z_simul_5; global z_simul_10; global z_simul_20;
global dose_simul_5;  global dose_simul_10;  global dose_simul_20;
global flag_pdd;    global flag_ocr;
global flag;    %global fit_step

%z_measure=z_meas;   dose_measure=dose_meas; z_simulation=z_simul;   dose_simulation=dose_simul;

z_simulation=z_simul;
z_simulation_5=z_simul_5;   z_simulation_10=z_simul_10;   z_simulation_20=z_simul_20;
z_measure=z_meas;
z_measure_5=z_meas_5;   z_measure_10=z_meas_10;   z_measure_20=z_meas_20;
dose_measure=dose_meas; dose_simulation=dose_simul;
dose_simulation_5=dose_simul_5;  dose_measure_5=dose_meas_5; 
dose_simulation_10=dose_simul_10;  dose_measure_10=dose_meas_10;
dose_simulation_20=dose_simul_20;  dose_measure_20=dose_meas_20;

if flag_pdd==1
    [G,passage]=func_gamma_index(z_measure,dose_measure,z_simulation,dose_simulation,dta,dose_diff);
    cla(handles.axes1,'reset');
    axes(handles.axes1);
    z=plot(z_measure,dose_measure,'.k'); hold on
    plot(z_simulation,dose_simulation,'.b'); hold on
    [t,z1,z2]=plotyy(z_simulation,dose_simulation,z_simulation,G);
    b=[t,z1,z2];
    set(z1,'LineStyle','.','Color','b');
    set(z2,'Color','r');
    legend([z,z1,z2],'Dose measure','Dose simulation',strcat('Gamma Index',' %pass= ',num2str(passage),'%'))
    title('PDD AFTER FIT');
    xlabel('Depth');    ylabel('Dose');
    ylabel(b(2),'Gamma index');
    grid on

%axes(handles.axes3);
%plot(z_simulation,G,'r');
%legend(strcat('Gamma Index',' %pass= ',num2str(passage),'%'))
end

if flag_ocr==1
    [G_5,passage_5]=func_gamma_index(z_measure_5,dose_measure_5,z_simulation_5,dose_simulation_5,dta,dose_diff);
    [G_10,passage_10]=func_gamma_index(z_measure_10,dose_measure_10,z_simulation_10,dose_simulation_10,dta,dose_diff);
    [G_20,passage_20]=func_gamma_index(z_measure_20,dose_measure_20,z_simulation_20,dose_simulation_20,dta,dose_diff);
    
    %%%%%%%%%%%%%%%%%%%%%% VE GAMMA INDEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag==1
        axes(handles.axes1);
        y=plot(z_measure_5,dose_measure_5,'.k'); hold on
        plot(z_simulation_5,dose_simulation_5,'.b'); hold on
        [x,y1,y2]=plotyy(z_simulation_5,dose_simulation_5,z_simulation_5,G_5);
        a=[x,y1,y2];
        set(y1,'LineStyle','.','Color','b');
        set(y2,'Color','r');
        legend([y,y1,y2],'Dose measure','Dose simulation',strcat('Gamma Index',' %pass= ',num2str(passage_5),'%'))
        title('OCR 5 AFTER FIT');
        xlabel('Depth');    ylabel('Dose');
        ylabel(a(2),'Gamma index');
        grid on
    end

%cla(handles.axes2,'reset');
    if flag==2
        axes(handles.axes1);
        z=plot(z_measure_10,dose_measure_10,'.k'); hold on
        plot(z_simulation_10,dose_simulation_10,'.b'); hold on
        [t,z1,z2]=plotyy(z_simulation_10,dose_simulation_10,z_simulation_10,G_10);
        b=[t,z1,z2];
        set(z1,'LineStyle','.','Color','b');
        set(z2,'Color','r');
        legend([z,z1,z2],'Dose measure','Dose simulation',strcat('Gamma Index',' %pass= ',num2str(passage_10),'%'))
        title('OCR 10 AFTER FIT');
        xlabel('Depth');    ylabel('Dose');
        ylabel(b(2),'Gamma index');
        grid on
    end

%cla(handles.axes3,'reset');
    if flag==3
        axes(handles.axes1);
        l=plot(z_measure_20,dose_measure_20,'.k'); hold on
        plot(z_simulation_20,dose_simulation_20,'.b'); hold on
        [v,l1,l2]=plotyy(z_simulation_20,dose_simulation_20,z_simulation_20,G_20);
        c=[v,l1,l2];
        set(l1,'LineStyle','.','Color','b');
        set(l2,'Color','r');
        legend([l,l1,l2],'Dose measure','Dose simulation',strcat('Gamma Index',' %pass= ',num2str(passage_20),'%'))
        title('OCR 20 AFTER FIT');
        xlabel('Depth');    ylabel('Dose');
        ylabel(c(2),'Gamma index');
        grid on
   end
end


%display(G);
 
% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes1,'reset'); 
set(handles.measure,'string',''); set(handles.simul,'string','');
set(handles.dta,'string',''); set(handles.dose_diff,'string',''); 
set(handles.fit_value_first,'string','');
set(handles.step,'string',''); set(handles.fit_value_last,'string',''); 
set(handles.text22,'visible','off'); set(handles.text23,'visible','off');
set(handles.pop,'visible','off');    set(handles.text27,'visible','off');

clear all
function dta_Callback(hObject, eventdata, handles)
% hObject    handle to dta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dta as text
%        str2double(get(hObject,'String')) returns contents of dta as a double


% --- Executes during object creation, after setting all properties.
global dta;
dta=get(hObject,'string');  dta=str2num(dta);

function dta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dose_diff_Callback(hObject, eventdata, handles)
% hObject    handle to dose_diff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dose_diff as text
%        str2double(get(hObject,'String')) returns contents of dose_diff as a double


% --- Executes during object creation, after setting all properties.
global dose_diff;
dose_diff=get(hObject,'string');  dose_diff=str2num(dose_diff);

function dose_diff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dose_diff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in input_measure_ocr.

function ocr_Callback(hObject, eventdata, handles)
% hObject    handle to ocr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pdd.
%clear all
global flag_pdd;    global flag_ocr;
flag_pdd=0; flag_ocr=1;

set(handles.dta,'string',''); set(handles.dose_diff,'string','');
set(handles.text22,'visible','off'); set(handles.text23,'visible','on');
set(handles.fit_value_first,'string',''); set(handles.fit_value_last,'string','');
set(handles.step,'string','');
%set(handles.pop,'visible','on');    set(handles.text27,'visible','on');

function pdd_Callback(hObject, eventdata, handles)
% hObject    handle to pdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%clear all
global flag_pdd;    global flag_ocr;

flag_pdd=1; flag_ocr=0;

set(handles.dta,'string',''); set(handles.dose_diff,'string','');
set(handles.text22,'visible','on'); set(handles.text23,'visible','off');
set(handles.fit_value_first,'string',''); set(handles.fit_value_last,'string','');
set(handles.step,'string','');
set(handles.pop,'visible','off');   set(handles.text27,'visible','off');

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close();
clear all


% --- Executes on selection change in pop.
function pop_Callback(hObject, eventdata, handles)
% hObject    handle to pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop


% --- Executes during object creation, after setting all properties.

global flag;    global flag_ocr;

global z_meas_5; global z_meas_10;    global z_meas_20;  
global dose_meas_5;  global dose_meas_10; global dose_meas_20;
global z_simul_5; global z_simul_10; global z_simul_20;
global dose_simul_5;  global dose_simul_10;  global dose_simul_20;

flag=0;

if flag_ocr==1
switch get(handles.pop,'value')
    case 1
        cla(handles.axes1,'reset'); 
        axes(handles.axes1); title('OCR 5');
        plot(z_meas_5,dose_meas_5,'.k');   hold on
        plot(z_simul_5,dose_simul_5,'.b');    %title('OCR');
        legend('Dose measure','Dose simulation');   grid on
        xlabel('Depth'); ylabel('Dose');
        flag=1;
        
    case 2
        cla(handles.axes1,'reset'); 
        axes(handles.axes1); title('OCR 10');
        plot(z_meas_10,dose_meas_10,'.k');   hold on
        plot(z_simul_10,dose_simul_10,'.b');    %title('OCR');
        legend('Dose measure','Dose simulation');   grid on
        xlabel('Depth'); ylabel('Dose');
        flag=2;
        
    case 3
        cla(handles.axes1,'reset'); 
        axes(handles.axes1); title('OCR 20');
        plot(z_meas_20,dose_meas_20,'.k');   hold on
        plot(z_simul_20,dose_simul_20,'.b');    title('OCR');
        legend('Dose measure','Dose simulation');   grid on
        xlabel('Depth'); ylabel('Dose');
        flag=3;
end
end

function pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
