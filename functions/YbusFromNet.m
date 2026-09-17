function [Ybus] = YbusFromNet(NetList,MutualNetList,NumBus)

%This function calculates Ybus by using a NetList and, when necessary a
%MutualNetList. It also needs the number os buses as an input. The form of
%the matrix are:

%NetList = 
%[Branch number  Bus number  Bus number  Impedance between the buses; ...];

%MutalNetList = 
%[Branch number  Branch number  Mutual Impedance; ...];


% Without mutual impedances, assemble branch admittances directly. This
% avoids a dense branch-impedance inverse and also supports one branch.
if isempty(MutualNetList) || isequal(MutualNetList, 0)
    Ybus = zeros(NumBus);
    for k=1:size(NetList,1)
        from = NetList(k,2);
        to = NetList(k,3);
        y = 1/NetList(k,4);
        if from ~= 0
            Ybus(from,from) = Ybus(from,from) + y;
        end
        if to ~= 0
            Ybus(to,to) = Ybus(to,to) + y;
        end
        if from ~= 0 && to ~= 0
            Ybus(from,to) = Ybus(from,to) - y;
            Ybus(to,from) = Ybus(to,from) - y;
        end
    end
    return;
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