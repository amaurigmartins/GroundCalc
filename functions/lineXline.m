function [X] = lineXline(A,B)

%Find intersection of N lines in D-dimensional space, in least squares sense.
% X = lineXline(A,B)     -line starts & ends as N*D
% X = lineXline({x y..})    -lines as D long cell of starts & ends as 2*N
% [X,P,R] = lineXline(..)      -extra outputs
% [X,P,R,x,t,l] = lineXline(..)   -plot outputs
%X: Intersection point, in least squares sense, as 1*D
%P: Nearest point to the intersection on each line, as N*D
%R: Distance from intersection to each line, as N*1
%x: Intersection point X as D-long cell of coordinates {x y..}
%p: Nearest points P as D-long cell of coordinates {x y..} as N*1
%l: Initial lines A-to-B as D-long cell of coordinates {x y..} as 2*N
%
%Remarks:
%-Lines are assumed to be infinite in both directions.
%-Finds point nearest to all lines using minimum sum of squared distances.
%-For parallel lines returns an arbitrary point and a warning.
%-For lines of length zero returns NaNs.
%-Comments/issues/corrections email Serge: s3rg3y@hotmail.com
%
%Example: (3 lines, 2 dimensions)
% [X,P,R,x,p,l] = lineXline([0 0;3 0;0 4],[5 5;0 5;5 2]);
% plot(x{:},'*k',p{:},'.k',l{:})
% X2 = lineXline(l) %alternate input form, same results
%
%Example: (2 lines, 3 dimensions)
% [X,P,R,x,p,l] = lineXline(rand(2,3),rand(2,3));
% plot3(x{:},'*k',p{:},'.k',l{:})
%
%Example: (4 lines, 5 dimensions)
% [X,P,R] = lineXline(rand(4,5),rand(4,5))
%
%See also: mldivide
 
%convert cell input {x y..} to A,B form
if iscell(A)
    A = permute(cat(3,A{:}),[2 3 1]); %2*N*D > N*D*2
    [A,B] = deal(A(:,:,1),A(:,:,2));
end
 
%find intersection
V = B - A; %vectors from A to B
V = bsxfun(@rdivide,V,sqrt(sum(V.*V,2))); %normalized vectors
[N,D] = size(A); %number of points & dimensions
T = bsxfun(@minus,bsxfun(@times,V',reshape(V,[1 N D])),reshape(eye(D),[D 1 D])); %V.*V-1 as D*N*D
S = reshape(sum(T,2),[D D]); %sum T along N, as D*D
C = reshape(T,[D N*D])*A(:); %T*A, as D*1
X = mldivide(S,C)'; %solve for X: S*X=C, in least squares sense
 
%checks
% if any(isnan(V(:))) %zero length lines
%     warning('lineXline:ZeroLengthLine','One or more lines with zero length.')
% elseif rcond(S)<eps*1000 %parallel lines, stackoverflow.com/questions/13145948/how-to-find-out-if-a-matrix-is-singular
%     warning('lineXline:ParallelLines','Lines are near parallel.')
% end
 
%extra outputs
if nargout>=2
    U = sum(bsxfun(@times,bsxfun(@minus,X,A),V),2); %dot(X-A,V) distance from A to nearest point on each line
    P = A + bsxfun(@times,U,V); %nearest point on each line
end
if nargout>=3
    R = sqrt(sum(bsxfun(@minus,X,P).^2,2)); %distance from intersection to each line
end
 
%plot outputs
if nargout>=4
    x = num2cell(X); %intersection point X
end
if nargout>=5
    p = num2cell(P,1); %tangent points P
end
if nargout>=6
    l = mat2cell([A(:) B(:)]',2,ones(1,D)*N); %initial lines A,B using cell format {x y..}
end

end