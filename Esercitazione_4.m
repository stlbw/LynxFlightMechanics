%Helicopter Flight Mechanics
%Esercitazione 4 - Lynx
%----------- Created by Matheus Padilha -----------%
clear
close all
clc
%----------- Input -----------%
V = [0 40 80]'; %velocity in kn
V = V*0.5144; %velocity in m/s


disp('Helicopter Flight Mechanics')
disp('Esercitazione 4 - Lynx')
disp(' ')
disp('----Helicopter Dynamic Modes----')
disp(' ')

%% Longitudinal and latero-directional planes not coupled

[A_Long, A_Lat, A_coupled] = stability_derivatives_Lynx(V);

[eigvec_long, eigval_long, eigvec_lat, eigval_lat, eigvec_long_norm, eigvec_lat_norm, sys_long, sys_lat]=long_lat_ncoupled(A_Long,A_Lat);

%LONGITUDINAL PLANE
%Hover
eigvec_hover_long = eigvec_long(:,:,1);
eigval_hover_long = diag(eigval_long(:,:,1));

%40 kn
eigvec_40_long = eigvec_long(:,:,2);
eigval_40_long = diag(eigval_long(:,:,2));

%80 kn
eigvec_80_long = eigvec_long(:,:,3);
eigval_80_long = diag(eigval_long(:,:,3));

for i=1:length(V)
    eigval_long_sum(:,i) = diag(eigval_long(:,:,i)); 
end
figure()
pzplot(sys_long(:,:,1,1),sys_long(:,:,2,1),sys_long(:,:,3,1))
legend('0 kt', '40 kt' ,'80 kt','Location','best')
title('Longitudinal plane - uncoupled')
grid on

% eigval_long = [eigval_hover_long eigval_40_long eigval_80_long] %eigenvalues for each velocity (columns)
%LONGITUDINAL PLANE
%Hover
eigvec_hover_lat = eigvec_lat(:,:,1);
eigval_hover_lat = diag(eigval_lat(:,:,1));

%40 kn
eigvec_40_lat = eigvec_lat(:,:,2);
eigval_40_lat = diag(eigval_lat(:,:,2));

%80 kn
eigvec_80_lat = eigvec_lat(:,:,3);
eigval_80_lat = diag(eigval_lat(:,:,3));

for i=1:length(V)
    eigval_lat_sum(:,i) = diag(eigval_lat(:,:,i)); 
end
figure()
pzplot(sys_lat(:,:,1,1),sys_lat(:,:,2,1),sys_lat(:,:,3,1))
legend('0 kt', '40 kt' ,'80 kt','Location','best')
title('Latero-directional plane - uncoupled')
grid on
%% Coupled Planes 

[eigvec_coup, eigval_coup, eigvec_norm_coup, sys_coup] = long_lat_coupled(A_coupled);

figure()
pzplot(sys_coup(:,:,1,1),sys_coup(:,:,2,1),sys_coup(:,:,3,1))
legend('0 kt', '40 kt' ,'80 kt','Location','best')
title('Lynx flight dynamics - coupled planes')
grid on