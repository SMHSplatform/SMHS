function [GLRLMout]=GLRLM_Features(IMG)
%--------------------------------------------------------------------------
% this program select a roi, qunatize to lower bit level and computing 
% gray level run length matrix and seven texture parameters viz., 
%    1. short run emphasis (sre) 
%    2. long run emphasis(lre)
%    3. gray level non-uniformity (gln)
%    4. gray level non-uniformity normalized(glnn)
%    5. run percentage (rp)
%    6. run length non-uniformity (rln)
%    7. run length non-uniformity normalized (rlnn)
%    8. run variance (rv)
%    9. gray level variance (glv)
%    10. low gray level run emphasis (lgre)
%    11. high gray level run emphasis (hgre)
%    12. Short Run Low Gray Level Emphasis (SRLGLE)
%    13. Short Run High Gray Level Emphasis (SRHGLE)
%    14. Long Run Low Gray Level Emphasis (LRLGLE)
%    15.  Long Run High Gray Level Emphasis (LRHGLE)
%--------------------------------------------------------------------------
% author: courage
%--------------------------------------------------------------------------
GLRLMout.sre=zeros(1,1); % short run emphasis
GLRLMout.lre=zeros(1,1); % long run emphasis
GLRLMout.gln=zeros(1,1); % gray level non-uniformity
GLRLMout.glnn=zeros(1,1); % Gray Level Non-Uniformity Normalized (GLNN)
GLRLMout.rp=zeros(1,1); % run percentage
GLRLMout.rln=zeros(1,1); % run length non-uniformity
GLRLMout.rlnn=zeros(1,1); % run length non-uniformity normalized
GLRLMout.rv=zeros(1,1); % run variance
GLRLMout.glv=zeros(1,1); % run variance
GLRLMout.lgre=zeros(1,1); % low gray level run emphasis
GLRLMout.hgre=zeros(1,1); % high gray level run emphasis
GLRLMout.srlgle=zeros(1,1); % run variance
GLRLMout.srhgle=zeros(1,1); % run variance
GLRLMout.lrlgle=zeros(1,1); % low gray level run emphasis
GLRLMout.lrhgle=zeros(1,1); % high gray level run emphasis

% --------- image input------------------
im=IMG;
%figure
%imshow(im)
%im1=imcrop(im);
%im2=im1(1:128,1:128);
%im2=double(im2);
%[m,n]=size(im2);
im2=im;
[m,n]=size(im2);
% --------- image quantization to 4 bits (16 gray levels)------------------
imax=max(max(im2));
imin=min(min(im2));
newim=im2-imin;
nmax=max(max(newim));
nmin=min(min(newim));
q=round(nmax/16);
[m,n]=size(newim);
quant=0;
for i=1:m
    for j=1:n
        k = newim(i,j);
        for b = 1:16
            if (i>quant)&(i<=quant+q)
                newim(i,j)=b/16;
                quant=quant+q;
            end            
        end
    end
end
newmax=max(max(newim));
newim1=newim/newmax;
newim2=round(newim1*16)+1;
dir=0; 
dist1=1;
if (dir == 1)
    newim2=newim2';
end
mx = max(max(newim2));
mn = min(min(newim2));
gl = (mx-mn)+1;
[p,q] = size(newim2);
n=p*q;
count=1;
c=1;
col=1;
grl(mx,p)=0;
maxcount(p*q)=0;
mc=0;
%---------------------computing gray level run length matrix---------------

for j=1:p
    for k=1:q-dist1
        mc=mc+1;
        g=newim2(j,k);
        f=newim2(j,k+dist1);
        if (g==f)&(g~=0)
            count=count+1;
            c=count;
            col=count;
            maxcount(mc)=count;
        else grl(g,c)=grl(g,c)+1;
            col=1;
            count=1;
            c=1;
        end
        grl(f,col)=grl(f,col)+1;
        count=1;
        c=1;
    end   
end
i=(mx:mn);
m=grl(mn:mx,:);
m1=m';
maxrun=max(max(maxcount));
s=0;
g(gl)=0;
r(p)=0;
for u=1:gl
    for v=1:p
        g(u)=g(u)+m(u,v);
        s=s+m(u,v);
    end
end
for u1=1:p
    for v1=1:gl
        r(u1)=r(u1)+m1(u1,v1);
    end
end
[dim,dim1]=size(g);
sre=0; lre=0; gln=0; rln=0; rp=0; lgre=0; hgre=0; u1=0; rv=0; u2=0;glv=0;
for h1=1:maxrun
    sre=sre+(r(h1)/(h1*h1));
    lre=lre+(r(h1)*(h1*h1));
    rln=rln+(r(h1)*r(h1));
    rp=rp+r(h1);
    u1=u1+(r(h1)*h1);
    rv=rv+(r(h1)*(h1-u1)*(h1-u1));
end
sre1=sre/s;
lre1=lre/s;
rln1=rln/s;
rlnn1=rln1/s; %归一化
rp1=rp/n;
for h2=1:gl
    gln=(gln+g(h2)^2);
    lgre=lgre+(g(h2)/(h2*h2));
    hgre=hgre+(h2*h2)*g(h2);
    u2=u2+(g(h2)*h2);
    glv=glv+(g(h2)*(h2-u2)*(h2-u2));
end
gln1=gln/s;
glnn1=gln1/s; %归一化
lgre1=lgre/s;
hgre1=hgre/s;

srlgle=0; srhgle=0;lrlgle=0;lrhgle=0;
for h1=1:maxrun
    for h2=1:gl
        srlgle=srlgle+m(h2,h1)/((h2*h2)*(h1*h1));
        srhgle=srhgle+m(h2,h1)*(h2*h2)/(h1*h1);
        lrlgle=lrlgle+m(h2,h1)*(h1*h1)/(h2*h2);
        lrhgle=lrhgle+m(h2,h1)*((h2*h2)*(h1*h1));
    end
end
srlgle1=srlgle/s; %归一化
srhgle1=srhgle/s; 
lrlgle1=lrlgle/s;
lrhgle1=lrhgle/s;

% ---------------------------display the parameters------------------------
%disp(sprintf('%6.4f',sre1))
%disp(sprintf('%6.4f',lre1))
%disp(sprintf('%6.4f',gln1))
%disp(sprintf('%6.4f',rp1))
%disp(sprintf('%6.4f',rln1))
%disp(sprintf('%6.4f',lgre1))
%disp(sprintf('%6.4f',hgre1))

GLRLMout.sre=sre1; % short run emphasis
GLRLMout.lre=lre1; % long run emphasis
GLRLMout.gln=gln1; % gray level non-uniformity
GLRLMout.glnn=glnn1; % Gray Level Non-Uniformity Normalized (GLNN)
GLRLMout.rp=rp1; % run percentage
GLRLMout.rln=rln1; % run length non-uniformity
GLRLMout.rlnn=rlnn1; % run length non-uniformity normalized
GLRLMout.rv=rv; % run variance
GLRLMout.glv=glv; % run variance
GLRLMout.lgre=lgre1; % low gray level run emphasis
GLRLMout.hgre=hgre1; % high gray level run emphasis
GLRLMout.srlgle=srlgle1; % run variance
GLRLMout.srhgle=srhgle1; % run variance
GLRLMout.lrlgle=lrlgle1; % low gray level run emphasis
GLRLMout.lrhgle=lrhgle1; % high gray level run emphasis
