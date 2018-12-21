function [x,gradx,iter] = UNINNPBBNLS(A,B,x,tol,maxiter)
% Projected Barzilai-Borwein method for the Nonnegative 
% Least Squares Problem
% using GLL nonmonotone line search.
%  Coded on 10/21/2005, L.H.
%  
%  Modified on 11/19/2006
%  Modified on 11/26,27,28/2006
%
mm=10;  %mm nonmonotone line search parameter
lamax=10^20; lamin=10^-20;
gamma=10^-4; 
%initt = cputime;
ABt=A*B'; BBt=B*B';
gradx = x*BBt - ABt;
delta0 = -sum(sum(x*ABt'));
dQd0 = sum(sum(x*(BBt*x')));
%btb=sum(sum(b.*b));
%Remark: Exact function value is: f0=0.5*btb+delta0+0.5*dQd0;
%but using f0=delta0+0.5*dQd0 won't affect the result and it
% saves time if b has a large size.
 f0=delta0+0.5*dQd0;

%f0=0.5*(norm(A-UVWG))^2;
%initgrad = norm(gradx,'fro');
%fprintf(1,'Init gradient norm %f\n', initgrad); 
for iter=1:maxiter,
% stopping condition
  projnorm = norm(gradx(gradx<0 | x>0),'fro');
% Check the norm of the projected gradient 
    if projnorm<tol,
    break;
  end
  if iter==1
      func(iter)=f0;
  else
      func(iter)=fn;
  end
  jj=min(iter-1,mm-1); 
  fmax=max(func(iter-jj:iter));
% Update x
  if iter==1
  lambda=1/max(max(abs(gradx)));
  end
  xn= max(x - lambda*gradx, 0); 
  dx=xn-x;
  delta = sum(sum(dx.*gradx));
  dQd = sum(sum(dx*(BBt*dx')));
  alpha=1;
  fn=func(iter)+alpha*delta+0.5*alpha^2*dQd; 
%  alphat=-delta/dQd;
  while (fn > fmax +alpha*gamma*delta)
  % Use Backtracking Line Search
     alpha=alpha/4;
     fn=func(iter)+alpha*delta+0.5*alpha^2*dQd;
  end
  xn=x+alpha*dx;
  
% Compute the BB steplength 
  gradxn=xn*BBt-ABt;
  sx = xn-x;  yx=gradxn-gradx;
  x=xn; 
  gradx=gradxn;
  sts=sum(sum(sx.*sx));
  sty=sum(sum(sx.*yx));
  if sty <= 0
      lambda=lamax;
  else
      lambda=min(lamax,max(lamin,sts/sty));
  end
end

%usedtime=cputime-initt;
%fprintf(1,'\nIter = %d || Final proj-grad norm %f || UsedCPUtime %f\n', iter, projnorm, usedtime);