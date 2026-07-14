function [CI, Chi2, alpha2] = PCR_conditional_independence_pvalue(X, Y, S, Data, test, alpha, ns)
if nargin < 5, test = 'LRT'; end
if nargin < 6, alpha = 0.05; end
if nargin < 7, ns = max(Data'); end
Data=Data';
N = size(Data,1);
qi=ns(S);
tmp=[1 cumprod(qi(1:end-1))];
qs=1+(qi-1)*tmp';
if isempty(qs),
    nij=zeros(ns(X),ns(Y));
    df=prod(ns([X Y])-1)*prod(ns(S));
else
    nijk=zeros(ns(X),ns(Y),1);
    tijk=zeros(ns(X),ns(Y),1);
    df=prod(ns([X Y])-1)*qs;
end
if (N<1*df)
    Chi2=-1;
    CI=0;
    alpha2 = 0;
elseif isempty(S)
    for i=1:ns(X),
        for j=1:ns(Y),
            nij(i,j)=length(find((Data(:,X)==i)&(Data(:,Y)==j))) ;
        end
    end
    restr=find(sum(nij,1)==0);
    if ~isempty(restr)
        nij=nij(:,find(sum(nij,1)));
    end
    tij=sum(nij,2)*sum(nij,1)/N ;
 switch test
    case 'pearson',
        tmpij=nij-tij;
        [xi yj]=find(tij<10);
        for i=1:length(xi),
           tmpij(xi(i),yj(i))=abs(tmpij(xi(i),yj(i)))-0.5;
        end
        warning off;
        tmp=(tmpij.^2)./tij;
        warning on;
        tmp(find(tmp==Inf))=0;
    case 'LRT',
        warning off;
        tmp=nij./tij;
        warning on;
        tmp(find(tmp==Inf | tmp==0))=1;
        tmp(find(tmp~=tmp))=1;
        tmp=2*nij.*log(tmp);
    otherwise,
        error(['unrecognized test ' test]);
    end
    Chi2=sum(sum(tmp));
    alpha2=1-chisquared_prob(Chi2,df);
    CI=(alpha2>=alpha) ;
else
    SizeofSSi=1;
    for exemple=1:N,
        i=Data(exemple,X);
        j=Data(exemple,Y);
        Si=Data(exemple,S)-1;
        if exemple==1
            SSi(SizeofSSi,:)=Si;
            nijk(i,j,SizeofSSi)=1;
        else
            flag=0;
            for iii=1:SizeofSSi
                if isequal(SSi(iii,:),Si)
                    nijk(i,j,iii)=nijk(i,j,iii)+1;
                    flag=1;
                end
            end
            if flag==0
                SizeofSSi=SizeofSSi+1;
                SSi(SizeofSSi,:)=Si;
                nijk(i,j,SizeofSSi)=1;
            end
        end
    end
    nik=sum(nijk,2);
    njk=sum(nijk,1);
    N2=sum(njk);
    for k=1:SizeofSSi
        if N2(:,:,k)==0
            tijk(:,:,k)=0;
        else
            tijk(:,:,k)=nik(:,:,k)*njk(:,:,k)/N2(:,:,k);
        end
    end
    switch test
    case 'pearson',
        tmpijk=nijk-tijk;
        [xi yj]=find(tijk<10);
        for i=1:length(xi),
            tmpijk(xi(i),yj(i))=abs(tmpijk(xi(i),yj(i)))-0.5;
        end
        warning off;
        tmp=(tmpijk.^2)./tijk;
        warning on;
        tmp(find(tmp==Inf))=0;
    case 'LRT',
        warning off;
        tmp=nijk./tijk;
        warning on;
        tmp(find(tmp==Inf | tmp==0))=1;
        tmp(find(tmp~=tmp))=1;
        tmp=2*nijk.*log(tmp);
    otherwise,
        error(['unrecognized test ' test]);
    end
    Chi2=sum(sum(sum(tmp)));
    alpha2=1-chisquared_prob(Chi2,df);
    CI=(alpha2>=alpha) ;
end
clear tijk
clear nijk
clear nij
clear tij
clear tmpijk
