function [eigvec_coup, eigval_coup, eigvec_norm_coup, sys_coup] = long_lat_coupled(A_coupled)
%LONG_LAT_COUPLED Returns the eigenvalues and eigenvectors for the
%longitudinal and latero-directional planes COUPLED
    %loop for computing eigenvalues and eigenvectors of [A], not coupled 
    for i=1:size(A_coupled,3)
        [eigvec_coup(:,:,i),eigval_coup(:,:,i)]=eig(A_coupled(:,:,i));

    % Normalized eigenvectors -> normalize each eigenvector with the "biggest"
    % component
        for m=1:length(A_coupled)
            eigvec_norm_coup(:,m,i) = abs(eigvec_coup(:,m,i))./max(abs(eigvec_coup(:,m,i)));
        end
    end

    B_coup = zeros(length(A_coupled),4);
    C_coup=ones(length(A_coupled),length(A_coupled));
    D_coup=zeros(length(A_coupled),4);
    sys_coup = ss(A_coupled,B_coup,C_coup,D_coup);
end