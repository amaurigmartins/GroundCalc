function out = enforce_symm(R)
TOL = 1e-6;
[n, m] = size(R);
if n ~= m
    error('Hand me a square matrix, not this irregular monstrosity.');
end

for i = 1:n
    for j = i+1:n
        if abs(R(i,j) - R(j,i)) > TOL
            avg_value = mean([R(i,j), R(j,i)]);
            R(i,j) = avg_value;
            R(j,i) = avg_value;
        end
    end
end

out=R;
end