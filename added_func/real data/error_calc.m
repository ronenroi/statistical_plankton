a=abs(recon_proj-unscaled_projections);
b=squeeze(max(unscaled_projections));
b=squeeze(max(b));
c=squeeze(sum(a));
c=squeeze(sum(c));
d=c./b;
e=sum(d);
e=e/157/132
f=d/157/132;
std(f)

