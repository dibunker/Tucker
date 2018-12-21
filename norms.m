function [fnorm,pnorm] = norms(A,varargin)
%NORMS find final norm and projected gradient norm
if size(varargin,2)==2
    
    factors=varargin{1}; N=size(factors,2);
    G=varargin{2};
    fnorm=norm(A-ttm(G,factors));
    
    for i=1:N, mat{i}=factors{i}'*factors{i}; end
    gradG = ttm(G,mat) - ttm(A,factors,'t'); gradG=double(gradG);
    G=double(G);
    pnorm = norm(gradG(gradG<0 | G>0),'fro')^2;
    G=tensor(G); allprod = ttm(G,factors);
    for i=1:N
        mode=[1:(i-1) (i+1):N];
        grad{i} = ttt(allprod,ttm(G,factors,mode),mode) - ttt(A,ttm(G,factors,mode),mode);grad{i} = double(grad{i});
        pnorm=pnorm + norm(grad{i},'fro')^2;
    end
     pnorm = sqrt(pnorm);
     
else
    
    factors=varargin{1};
    r=size(factors{1},2); N=size(factors,2); pnorm=0;


%     Ap=ones([1 r])*khatrirao(factors,'r')';
%     Ap=reshape(Ap,size(A));
%     fnorm=norm(A-Ap);

    G=zeros([r r r]);
    for i=1:r, G(i,i,i)=1; end
    G=tensor(G);
    Ap=ttm(G,factors);
    fnorm=norm(A-Ap);

    for i=1:N
        mode=[(i+1):N 1:(i-1)];
        Au=tenmat(A,i,'bc'); Au=double(Au);
        B=khatrirao(factors{mode});
        grad{i} = factors{i}*B'*B - Au*B;
        pnorm=pnorm + norm(grad{i}(grad{i}<0 | factors{i}>0),'fro')^2;
    end
    pnorm = sqrt(pnorm);
end