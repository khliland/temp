function X = GOuter(x)
X = x{1};
if length(x) > 1
    for i=2:length(x)
        if i==2
            X = GMPmem(X, x{i}',1);
        else
            X = GMP(X, x{i},0); %OBS - hva skjer her i det generelle tilfellet?
        end
    end
end