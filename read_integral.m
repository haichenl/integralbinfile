clear;
fileID = fopen('integral.bin');
natom = fread(fileID,1,'int');

A1 = fread(fileID,3,'int');
nh1 = A1(1);
nrow1 = A1(2);
ncol1 = A1(3);
s = fread(fileID,[nrow1,ncol1],'double');

A2 = fread(fileID,3,'int');
nh2 = A2(1);
nrow2 = A2(2);
ncol2 = A2(3);
ke = fread(fileID,[nrow2,ncol2],'double');

nbasis = nrow1;
en = zeros(nbasis, nbasis, natom);
for i=1:natom
    A3 = fread(fileID,3,'int');
    nh3 = A3(1);
    nrow3 = A3(2);
    ncol3 = A3(3);
    en(:,:,i) = fread(fileID,[nrow3,ncol3],'double');
end

nuniq = (nbasis.*(nbasis+1).*(nbasis.^2+nbasis+2))./8;
h2 = zeros(nbasis, nbasis, nbasis, nbasis);
for i=1:nuniq
    ind = fread(fileID,4,'int');
    p = ind(1) + 1;
    q = ind(2) + 1;
    r = ind(3) + 1;
    s = ind(4) + 1;
    tmp = fread(fileID,1,'double');
    h2(p, q, r, s) = tmp;
    h2(r, s, p, q) = tmp;
    h2(p, q, s, r) = tmp;
    h2(s, r, p, q) = tmp;
    h2(q ,p, r, s) = tmp;
    h2(r, s, q, p) = tmp;
    h2(q, p, s, r) = tmp;
    h2(s, r, q, p) = tmp;
end

fclose(fileID);
