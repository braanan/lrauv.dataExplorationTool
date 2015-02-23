function cleanx = kickout( x )

x(abs(x)>180)=NaN;

% mid = nanmedian(x);
% sig = nanstd(x);
cleanx = x;

% cleanx(cleanx > (mid+10*(sig)) | cleanx < (mid-10*(sig))) = NaN;

end