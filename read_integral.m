clear;
fileID = fopen('NOTINTEGRAL.BIN');
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
for iatom=1:natom
    A3 = fread(fileID,3,'int');
    nh3 = A3(1);
    nrow3 = A3(2);
    ncol3 = A3(3);
    en(:,:,iatom) = fread(fileID,[nrow3,ncol3],'double');
end

nuniq = (nbasis.*(nbasis+1).*(nbasis.^2+nbasis+2))./8;
h2 = zeros(nbasis, nbasis, nbasis, nbasis);
for iuni=1:nuniq
    ind = fread(fileID,4,'int');
    i = ind(1) + 1;
    j = ind(2) + 1;
    k = ind(3) + 1;
    l = ind(4) + 1;
    tmp = fread(fileID,1,'double');
    h2(i, j, k, l) = tmp;
    h2(k, l, i, j) = tmp;
    h2(i, j, l, k) = tmp;
    h2(l, k, i, j) = tmp;
    h2(j ,i, k, l) = tmp;
    h2(k, l, j, i) = tmp;
    h2(j, i, l, k) = tmp;
    h2(l, k, j, i) = tmp;
end

e_nuc = fread(fileID,1,'double');
do_env = fread(fileID,1,'int');
if(do_env)
    num_ptq = fread(fileID,1,'int');
    env = zeros(nbasis, nbasis, num_ptq);
    for i_ptq=1:num_ptq
        dimtmp = fread(fileID,3,'int');
        env(:,:,i_ptq) = fread(fileID,[nrow3,ncol3],'double');
    end
end

fclose(fileID);
