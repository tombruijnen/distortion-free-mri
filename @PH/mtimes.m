function output = mtimes(ph,data) 

if ph.adjoint==1 % A' - adjoint
    % Preallocate
    output=zeros(ph.Id);
    for avg=1:Id(12)
    for ex2=1:Id(11)
    for ex1=1:Id(10)
    for mix=1:Id(9)
    for loc=1:Id(8)
    for ech=1:Id(7)
    for ph=1:Id(6)
    for dyn=1:Id(5)
        output(:,:,:,:,dyn,ph,ech,loc,mix,ex1,ex2,avg)=exp(-1j*dyn*ph.A).*data(:,:,:,:,1);
    end
    end
    end
    end
    end
    end
    end
    end
else % A - forward
    % Preallocate
    output=zeros([ph.Id(1:4) 1 ph.Id(6:end)]);
    for avg=1:Id(12)
    for ex2=1:Id(11)
    for ex1=1:Id(10)
    for mix=1:Id(9)
    for loc=1:Id(8)
    for ech=1:Id(7)
    for ph=1:Id(6)
    for dyn=1:Id(5)
        output(:,:,:,:,1,ph,ech,loc,mix,ex1,ex2,avg)=output(:,:,:,:,1,ph,ech,loc,mix,ex1,ex2,avg)+...
            exp(1j*dyn*ph.A).*data(:,:,:,:,dyn,ph,ech,loc,mix,ex1,ex2,avg);        
    end
    end
    end
    end
    end
    end
    end
    end

end


% END
end  
