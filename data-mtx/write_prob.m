function write_prob(name,A,B,L)
[ai,aj,av] = find(A);
arrays2mtxfile([name '-A.mtx'],ai,aj,av);

[ai,aj,av] = find(B);
arrays2mtxfile([name '-B.mtx'],ai,aj,av);

[ai,aj,av] = find(L);
arrays2mtxfile([name '-L.mtx'],ai,aj,av);
