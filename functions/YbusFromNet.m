function [Ybus] = YbusFromNet(NetList,MutualNetList,NumBus)

%This function calculates Ybus by using a NetList and, when necessary a
%MutualNetList. It also needs the number os buses as an input. The form of
%the matrix are:

%NetList = 
%[Branch number  Bus number  Bus number  Impedance between the buses; ...];

%MutalNetList = 
%[Branch number  Branch number  Mutual Impedance; ...];


if MutualNetList == 0
    MutualNetList = [1 2 0];
end

%Now we can obtain the matrix Zpr

Zpr = zeros(size(NetList,1));

for i=1:size(NetList,1)     %Filling diagonal values
    
    Zpr(i,i) = NetList(i,4);
        
end

for k=1:size(MutualNetList,1) %Filling out of diagonal values
    
    i = MutualNetList(k,1);
    j = MutualNetList(k,2);
    
    Zpr(i,j) = MutualNetList(k,3);
    Zpr(j,i) = Zpr(i,j);
    
end

%Now we obtain Ypr from Zpr

Ypr = inv(Zpr);

%Now we need to calculate matrix A

A = zeros(size(NetList,1),NumBus);

for k=1:size(NetList,1)
    for j=1:NumBus
        
        if NetList(k,2)==j         %from
            
            A(k,j) = 1;
            
        elseif NetList(k,3)==j     %to
            
            A(k,j) = -1;
            
        end
        
    end
end  

%With Ypr and A we can calculate Ybus

Ybus = A'*Ypr*A;

end