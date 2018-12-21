function [factors,G,niter] = nnTuckern(A,G0,tol,timelimit,maxiter,varargin)

% Coded on 3/21/2008
% Alternative Projected Barzilai-Borwein Method 
% It calls PBBLSG.m, UNIPBBLS.m

factors=varargin;N=size(varargin,2);
niter=0;
initgrad=0;
G=G0;
Au=tenmat(A,1,'bc'); Au=double(Au);

for i=1:N
    mode{i}=[1:(i-1) (i+1):N];
    Gu{i}=tenmat(G0,i,'bc'); Gu{i}=double(Gu{i});
    
%     ACt=[];
%     C=dnc(factors{mode{i}(1:end-1)});
%     k=1; s=size(factors{mode{i}(end)},2);
%     for j=1:size(C,2), ACt(:,k:(s*j))=Au{i}*dnc(C(:,j),factors{mode{i}(end)}); k=j*s + 1; end
%     ABt=ACt*Gu{i}';


    ABt=ttt(A,ttm(G0,varargin,mode{i}),mode{i});ABt = double(ABt);
    grad{i} = factors{i}*(Gu{i}*dnc(factors{mode{i}},'ttn')*Gu{i}') - ABt;
%     grad{i} = ttt(allprod,ttm(G0,varargin,mode{i}),mode{i}) - ttt(A,ttm(G0,varargin,mode{i}),mode{i});grad{i} = double(grad{i});
    initgrad=initgrad + norm(grad{i},'fro')^2;
end


U=factors{1};
C=dnc(factors{2:end-1});
ACt=U'*Au;
k=1;  s=size(factors{end},2);
for i=1:size(C,2), BtACt(:,k:(s*i))=ACt*dnc(C(:,i),factors{end}); k=i*s +1; end
BtB=U'*U;
for i=2:N, fac{i}=factors{i}'*factors{i}; end
CCt=dnc(fac{2:end});
gradG = -BtACt + BtB*Gu{1}*CCt;

initgrad=initgrad + norm(gradG(gradG<0 | Gu{1}>0),'fro')^2;
initgrad = sqrt(initgrad);




fprintf(1,'Init proj gradient norm %f\n', initgrad); 

tols{1} = max(0.001,tol)*initgrad; tolG=tols{1};
for i=2:N
    tols{i} = tols{1};
end


%Main Loop
initt = cputime;
for iter=1:maxiter,

  
  for i=1:N
      B=ttm(G,factors,mode{i}); B=tenmat(B,i,'bc'); B=double(B);
      Au=tenmat(A,i,'bc'); Au=double(Au);
      [factors{i},grad{i},iters{i}] = UNIPBBNLS(Au,B,factors{i},tols{i},1000);
      if iters{i}==1
          tols{i} = 0.1 * tols{i};
      end
      if iters{i}==1000
          fprintf(1,'Max iter in UNIPBBNLS in %d\n',i);
      end
  end
  
  Au=tenmat(A,1,'bc'); Au=double(Au);
  [G,gradG,iterG] = PBBNLSGn(Au,Gu{1},factors,N,tolG,1000);
  if iterG==1,
      tolG = 0.1 * tolG;
  end
  Gu{1}=G;
  zx=size(G0); Nreshape=[N:-1:2]; reshapesize=[zx(1) zx(Nreshape)];
  G=ipermute(reshape(G,reshapesize),[1 Nreshape]); G=tensor(G);
%   for i=1:N, Gu{i}=tenmat(G,i,'bc'); Gu{i}=double(Gu{i}); end
  
  niter=niter + iterG;
  for i=1:N
      niter=niter + iters{i};
  end

  % stopping condition
  projnorm = norm(gradG(gradG<0 | Gu{1}>0),'fro')^2;
  for i=1:N
       projnorm=projnorm + norm(grad{i}(grad{i}<0 | factors{i}>0),'fro')^2;
  end
  projnorm = sqrt(projnorm);
  if projnorm < tol*initgrad | cputime-initt > timelimit,
    break;
  end
  
  
  
%%  if (iterW==1 & iterH==1 & tolH + tolW < tol*initgrad),
%%    fprintf(1,'Failed to move\n'); break;
%%  end
  if rem(iter,10)==0, fprintf(1,'.'); end
end
usedtime=cputime-initt;
fprintf(1,'\nIter = %d || Final proj-grad norm %f || UsedCPUtime %f\n', iter, projnorm, usedtime);