function fun = build_polynomial(x, y, n)

    p = polyfit(x,y,n);
    
    fun = @(x) 0;
    
    for iter = 1:n+1
        fun = @(x) fun(x) + p(iter) * x.^(n-iter+1);
    end
    
    return fun;
    
end