function y=SoftThresh(x,p)
y=(abs(x)-p).*sign(x).*(abs(x)>p);
y(isnan(y))=0;
end 