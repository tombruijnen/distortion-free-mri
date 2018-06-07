function krt = krt_matrix(N)
krt=[];
for n=1:N
    krt(:,n)=(1:N)+(n-1);
end
for n=1:numel(krt)
    if krt(n)>N
        krt(n)=krt(n)-N;
    end
end

% END
end