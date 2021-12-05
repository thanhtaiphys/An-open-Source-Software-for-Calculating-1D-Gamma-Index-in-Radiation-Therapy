function [z_fit,dose_fit]=fit_data(depth,dose,first_fit,step,last_fit)
%depth=0:1:30;  
%depth=transpose(depth);
%dose=depth.^2; 

chay=1; chay1=1;
z_fit=zeros(2,1);   dose_fit=zeros(2,1);

%%%%%%%%%%%%%%%%%% READ TO FIT %%%%%%%%%%%%%%%%%%%%%%%
fits=cell(size(depth,1)-1,1);
  for i=2:length(depth)
      for j=i-1:i
        z_fit(chay,1)=depth(j,1);
        dose_fit(chay,1)=dose(j,1);
        chay=chay+1;
      end
      chay=1;
      fits{chay1,1} = fit(z_fit,dose_fit,'linear');
      chay1=chay1+1;
  end
  %display(size(fits));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_new=first_fit:step:last_fit;  z_new=transpose(z_new);   pos=0;
dose_new=zeros(size(z_new));
for i=1:length(z_new)
    for j=1:length(depth)-1
        if z_new(i,1)>=depth(j,1) && z_new(i,1)<depth(j+1,1)
            pos=j;
        end
        if pos~=0
            if z_new(i,1)>depth(pos,1) && z_new(i,1)<depth(pos+1,1)
                fit_func=fits{pos,1};
                dose_new(i,1)=fit_func(z_new(i,1));
            end
            if z_new(i,1)==depth(pos,1)
                dose_new(i,1)=dose(j,1);
            end
            break;
       end
    end
    if z_new(i,1)==depth(length(depth),1)
        dose_new(i,1)=dose(length(depth),1);
        %pos=length(x);
    end
    %display(pos);
    %if x_new(i,1)==x(pos,1)
    %    y_new(i,1)=y(pos,1);
    %end
    pos=0;
end

z_fit=z_new;
dose_fit=dose_new;
%depth=transpose(depth); dose=transpose(dose); x_new=transpose(x_new); y_new=transpose(y_new);
%figure
%hold on
%plot(depth,dose,'r');    hold on
%plot(x_new,y_new,'b');

%clear all