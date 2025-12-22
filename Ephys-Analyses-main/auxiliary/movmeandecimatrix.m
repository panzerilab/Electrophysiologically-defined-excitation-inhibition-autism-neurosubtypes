function [ XmatResDec ] = movmeandecimatrix ( Xmat, decimfactor)

szs=size(Xmat);

% TIME HAS TO LAY ALONG 1st DIM
% FACTOR HAS TO DIVIDE 1st SIZE

if(rem(szs(1),decimfactor)~=0)
    error('Signal samples does not divide Decim Factor, check decimation does not overlap.');
end

Xmatmvm=movmean(Xmat,decimfactor,1);
XmatResDec=eval(['Xmatmvm(ceil(decimfactor/2):decimfactor:end' repmat(',:',1,length(szs)) ')']);

%plot(Xmat(:,1,1),'b*-'); hold on;
%plot(Xmatmvm(:,1,1),'ro-');
%stem(ceil(decimfactor/2):decimfactor:size(Xmat,1),XmatResDec(:,1,1),'gx-');

end