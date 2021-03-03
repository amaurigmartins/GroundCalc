function [ subdivcond ] = f( vetor_entrada, MAX_LEN )
% subdivcond: Subdivide o vetor_entrada em partições de comprimento
% inferior a MAX_LEN e retorna as coordenadas das extremidades, do ponto 
% médio, do comprimento e o raio de cada subdivisão.
%   Função adota os seguintes parâmetros de entrada:
%   xi, yi, zi - coordenadas da extremidade inicial do vetor a ser subdividido.
%   xf, yf, zf - coordenadas da extremidade final do vetor a ser subdividido.
%   r - raio do condutor (não é subdividido)
%   MAX_LEN - comprimento máximo de cada subdivisão.

vetor_saida=[0 0 0 0 0 0];

if size(vetor_entrada,2)==7
    
    for i=1:size(vetor_entrada,1)
        
        xi=vetor_entrada(i,1);
        xf=vetor_entrada(i,4);
        yi=vetor_entrada(i,2);
        yf=vetor_entrada(i,5);
        zi=vetor_entrada(i,3);
        zf=vetor_entrada(i,6);
		r=vetor_entrada(i,7);
        
        
        comprimento = sqrt((xf-xi)^2+(yf-yi)^2+(zf-zi)^2);
        
        if(comprimento <= MAX_LEN)
            if(vetor_saida(1,:)==[0 0 0 0 0 0])
                vetor_saida=[xi yi zi xf yf zf];
                L=[comprimento];
				R=[r];
            else
                vetor_saida=[vetor_saida; xi yi zi xf yf zf];
                L=[L; comprimento];
				R=[R; r];
            end
        else
            n_subdiv=ceil(comprimento/MAX_LEN);
            for j=0:n_subdiv-1
                
                xi_div=xi*(1-(j/n_subdiv))+xf*(j/n_subdiv);
                yi_div=yi*(1-(j/n_subdiv))+yf*(j/n_subdiv);
                zi_div=zi*(1-(j/n_subdiv))+zf*(j/n_subdiv);
                xf_div=xi*(1-((j+1)/n_subdiv))+xf*((j+1)/n_subdiv);
                yf_div=yi*(1-((j+1)/n_subdiv))+yf*((j+1)/n_subdiv);
                zf_div=zi*(1-((j+1)/n_subdiv))+zf*((j+1)/n_subdiv);
                comprimento_div = sqrt((xf_div-xi_div)^2+(yf_div-yi_div)^2+(zf_div-zi_div)^2);
                
                if(vetor_saida(1,:)==[0 0 0 0 0 0])
                    vetor_saida=[xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[comprimento_div];
					R=[r];
                else
                    vetor_saida=[vetor_saida; xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[L; comprimento_div];
					R=[R; r];
                end
            end
            
        end
        
    end
    
    vetor_mid = [(vetor_saida(:,1)+vetor_saida(:,4))./2 (vetor_saida(:,2)+vetor_saida(:,5))./2 (vetor_saida(:,3)+vetor_saida(:,6))./2];
    subdivcond = horzcat(vetor_saida,vetor_mid, L, R);
    
elseif size(vetor_entrada,2)==4
    
    
    for i=1:size(vetor_entrada,1)-1
        
        xi=vetor_entrada(i,1);
        xf=vetor_entrada(i+1,1);
        yi=vetor_entrada(i,2);
        yf=vetor_entrada(i+1,2);
        zi=vetor_entrada(i,3);
        zf=vetor_entrada(i+1,3);
		r=vetor_entrada(i,4);
        
        comprimento = sqrt((xf-xi)^2+(yf-yi)^2+(zf-zi)^2);
        
        if(comprimento <= MAX_LEN)
            if(vetor_saida(1,:)==[0 0 0 0 0 0])
                vetor_saida=[xi yi zi xf yf zf];
                L=[comprimento];
				R=[r];
            else
                vetor_saida=[vetor_saida; xi yi zi xf yf zf];
                L=[L; comprimento];
				R=[R; r];
            end
        else
            n_subdiv=ceil(comprimento/MAX_LEN);
            for j=0:n_subdiv-1
                
                xi_div=xi*(1-(j/n_subdiv))+xf*(j/n_subdiv);
                yi_div=yi*(1-(j/n_subdiv))+yf*(j/n_subdiv);
                zi_div=zi*(1-(j/n_subdiv))+zf*(j/n_subdiv);
                xf_div=xi*(1-((j+1)/n_subdiv))+xf*((j+1)/n_subdiv);
                yf_div=yi*(1-((j+1)/n_subdiv))+yf*((j+1)/n_subdiv);
                zf_div=zi*(1-((j+1)/n_subdiv))+zf*((j+1)/n_subdiv);
                comprimento_div = sqrt((xf_div-xi_div)^2+(yf_div-yi_div)^2+(zf_div-zi_div)^2);
                
                if(vetor_saida(1,:)==[0 0 0 0 0 0])
                    vetor_saida=[xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[comprimento_div];
					R=[r];
                else
                    vetor_saida=[vetor_saida; xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[L; comprimento_div];
					R=[R; r];
                end
            end
            
        end
        
    end
    
    
    vetor_mid = [(vetor_saida(:,1)+vetor_saida(:,4))./2 (vetor_saida(:,2)+vetor_saida(:,5))./2 (vetor_saida(:,3)+vetor_saida(:,6))./2];
    subdivcond = horzcat(vetor_saida,vetor_mid,L,R);
    
end
end

