function x = tls(A,b)

[~,n] = size (A);
C= [A b]; % C: augmented matrix
[~,~,V] = svd(C); % Matrix V (eq. 11)

 

VXY = V(1:n,1+n:end); % Submatrix: first n rows and the n+1 to last column

VYY = V(1+n:end,1+n:end); % Take the bottom-right block of V.

if VYY ==0
    disp ('ERROR: not TLS solution');
    x = zeros(n,1);
    return;
end

x = -VXY/VYY;