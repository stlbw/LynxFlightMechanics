function [eigvec_long, eigval_long, eigvec_lat, eigval_lat, eigvec_long_norm, eigvec_lat_norm, sys_long, sys_lat] = long_lat_ncoupled(A_long, A_lat)
%LONG_LAT_NCOUPLED Returns the eigenvalues and eigenvectors for the
%longitudinal and latero-directional planes NOT COUPLED

%loop for computing eigenvalues and eigenvectors of [A], not coupled 
    for i=1:size(A_long,3)
        [eigvec_long(:,:,i),eigval_long(:,:,i)]=eig(A_long(:,:,i));
        [eigvec_lat(:,:,i),eigval_lat(:,:,i)]=eig(A_lat(:,:,i));

    % Normalized eigenvectors -> normalize each eigenvector with the "biggest"
    % component
        for m=1:length(A_long)
            eigvec_long_norm(:,m,i) = abs(eigvec_long(:,m,i))./max(abs(eigvec_long(:,m,i)));
        end

        for m=1:length(A_lat)
            eigvec_lat_norm(:,m,i) = abs(eigvec_lat(:,m,i))./max(abs(eigvec_lat(:,m,i)));
        end
    end

    %Creation of a State-space model for both longitudinal and
    %latero-directional planes

    B_long = zeros(length(A_long),2);
    C_long=ones(length(A_long),length(A_long));
    D_long=zeros(length(A_long),2);
    
    B_lat=zeros(length(A_lat),2);
    C_lat=ones(length(A_lat),length(A_lat));
    D_lat=zeros(length(A_lat),2);
    
    sys_long=ss(A_long,B_long,C_long,D_long);
    sys_lat=ss(A_lat,B_lat,C_lat,D_lat);
end

