# GeneIndexer
This is a solution to Problem B of the Summer Camp Math Modeling, 2015 Shenzhen Cup. Solve it for my hobby.

#Problem Description

This problem is known as K-mer index problem (the problem description saids).

Given gene as several strings, each 100 chars in length, 1,000,000 strings in total. For given K, search a gene pattern length K in all gene.

eg.: gene = CTGTACTGTAT, K = 5, possible pattern can be: {CTGTA，TGTAC，GTACT，TACTG，ACTGT，TGTAT}

The problem is asking for:

* create index for quick searching such pattern, K is constant for each index.
* onces the index is created, query should be fast, and memory resage should be as lower as possible.
* analyse the complexity of time and memory of building index and query
* for 8GB memory, what's the maxium K supported
* a program is judeged by the following KPI (ordery by importance): query speed, index memory, maxium K supported for 8GB memroy, time for building index

And KPI of my solution (when K = 30 on 4 cores machine, you may check the code in main.cpp):
* query speed: 2us (worsst case), 1us in adverage.
* total memory: around 400MB.
* maxium K: 400MB is engouh for all possible K.
* time for building index: < 8s (including reading data files, single thread for reading file, and 8 threads for building index).
