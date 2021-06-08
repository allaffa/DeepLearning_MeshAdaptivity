function [sfun] = smooth_fun(fun,ncycles,upd)

  ndata = size(fun,1);
  sfun = fun;

  if(upd==1) %Gauss-Seidel type of update

    for ic = 1:ncycles
      sfun(1) = 0.5*(sfun(1)+sfun(2));
      for ii = 2:ndata-1
        sfun(ii) = 0.25*sfun(ii-1)+0.5*sfun(ii)+0.25*sfun(ii+1);
      end
      sfun(ndata) = 0.5*(sfun(ndata-1)+sfun(ndata));
    end

  elseif(upd==2) %Jacobi type of update

    for ic = 1:ncycles
      sfun(1) = 0.5*(fun(1)+fun(2));
      for ii = 2:ndata-1
        sfun(ii) = 0.25*fun(ii-1)+0.5*fun(ii)+0.25*fun(ii+1);
      end
      sfun(ndata) = 0.5*(fun(ndata-1)+fun(ndata));
      fun = sfun;
    end

  else %Err

    error('Error: smooth_fun (upd)');

  end



