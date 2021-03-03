function [ subdivcond ] = f( vec_in, MAX_LEN )


vec_out=[0 0 0 0 0 0];

if size(vec_in,2)==7
    
    for i=1:size(vec_in,1)
        
        xi=vec_in(i,1);
        xf=vec_in(i,4);
        yi=vec_in(i,2);
        yf=vec_in(i,5);
        zi=vec_in(i,3);
        zf=vec_in(i,6);
		r=vec_in(i,7);
        
        
        llen = sqrt((xf-xi)^2+(yf-yi)^2+(zf-zi)^2);
        
        if(llen <= MAX_LEN)
            if(vec_out(1,:)==[0 0 0 0 0 0])
                vec_out=[xi yi zi xf yf zf];
                L=[llen];
				R=[r];
            else
                vec_out=[vec_out; xi yi zi xf yf zf];
                L=[L; llen];
				R=[R; r];
            end
        else
            n_subdiv=ceil(llen/MAX_LEN);
            for j=0:n_subdiv-1
                
                xi_div=xi*(1-(j/n_subdiv))+xf*(j/n_subdiv);
                yi_div=yi*(1-(j/n_subdiv))+yf*(j/n_subdiv);
                zi_div=zi*(1-(j/n_subdiv))+zf*(j/n_subdiv);
                xf_div=xi*(1-((j+1)/n_subdiv))+xf*((j+1)/n_subdiv);
                yf_div=yi*(1-((j+1)/n_subdiv))+yf*((j+1)/n_subdiv);
                zf_div=zi*(1-((j+1)/n_subdiv))+zf*((j+1)/n_subdiv);
                llen_div = sqrt((xf_div-xi_div)^2+(yf_div-yi_div)^2+(zf_div-zi_div)^2);
                
                if(vec_out(1,:)==[0 0 0 0 0 0])
                    vec_out=[xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[llen_div];
					R=[r];
                else
                    vec_out=[vec_out; xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[L; llen_div];
					R=[R; r];
                end
            end
            
        end
        
    end
    
    vec_mid = [(vec_out(:,1)+vec_out(:,4))./2 (vec_out(:,2)+vec_out(:,5))./2 (vec_out(:,3)+vec_out(:,6))./2];
    subdivcond = horzcat(vec_out,vec_mid, L, R);
    
elseif size(vec_in,2)==4
    
    
    for i=1:size(vec_in,1)-1
        
        xi=vec_in(i,1);
        xf=vec_in(i+1,1);
        yi=vec_in(i,2);
        yf=vec_in(i+1,2);
        zi=vec_in(i,3);
        zf=vec_in(i+1,3);
		r=vec_in(i,4);
        
        llen = sqrt((xf-xi)^2+(yf-yi)^2+(zf-zi)^2);
        
        if(llen <= MAX_LEN)
            if(vec_out(1,:)==[0 0 0 0 0 0])
                vec_out=[xi yi zi xf yf zf];
                L=[llen];
				R=[r];
            else
                vec_out=[vec_out; xi yi zi xf yf zf];
                L=[L; llen];
				R=[R; r];
            end
        else
            n_subdiv=ceil(llen/MAX_LEN);
            for j=0:n_subdiv-1
                
                xi_div=xi*(1-(j/n_subdiv))+xf*(j/n_subdiv);
                yi_div=yi*(1-(j/n_subdiv))+yf*(j/n_subdiv);
                zi_div=zi*(1-(j/n_subdiv))+zf*(j/n_subdiv);
                xf_div=xi*(1-((j+1)/n_subdiv))+xf*((j+1)/n_subdiv);
                yf_div=yi*(1-((j+1)/n_subdiv))+yf*((j+1)/n_subdiv);
                zf_div=zi*(1-((j+1)/n_subdiv))+zf*((j+1)/n_subdiv);
                llen_div = sqrt((xf_div-xi_div)^2+(yf_div-yi_div)^2+(zf_div-zi_div)^2);
                
                if(vec_out(1,:)==[0 0 0 0 0 0])
                    vec_out=[xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[llen_div];
					R=[r];
                else
                    vec_out=[vec_out; xi_div yi_div zi_div xf_div yf_div zf_div];
                    L=[L; llen_div];
					R=[R; r];
                end
            end
            
        end
        
    end
    
    
    vec_mid = [(vec_out(:,1)+vec_out(:,4))./2 (vec_out(:,2)+vec_out(:,5))./2 (vec_out(:,3)+vec_out(:,6))./2];
    subdivcond = horzcat(vec_out,vec_mid,L,R);
    
end
end

