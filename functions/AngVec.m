function Ang = AngVec( A, B )

%This function returns the angle between two vectors A and B

Ang = acosd( ( dot(A,B) )/( norm(A)*norm(B) ) );

end