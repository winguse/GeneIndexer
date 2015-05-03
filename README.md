# GeneIndexer
This is a solution to Problem B of the Summer Camp Math Modeling, 2015 Shenzhen Cup. Solve it for my hobby.

##Problem Description

This problem is known as K-mer index problem (the problem description said).

Given gene as several strings, each 100 chars in length, 1,000,000 strings in total. For given K, search a gene pattern length K in all gene.

eg.: gene = CTGTACTGTAT, K = 5, possible pattern can be: {CTGTA，TGTAC，GTACT，TACTG，ACTGT，TGTAT}

The problem is asking for:

* create index for quick searching such pattern, K is constant for each index.
* once the index is created, query should be as faster as possible, and memory usage should be as lower as possible.
* analyse the complexity of time and memory of building index and query.
* for 8GB memory, what's the maximum K supported.
* a program is judged by the following KPI (order by importance): query speed, index memory, maximum K supported for 8GB memory, time for building index

And KPI of my solution (when K = 30 on 4 cores machine, you may check the code in main.cpp):
* query speed: 2us (worst case), 1us in average.
* total memory: around 400MB.
* maximum K: 400MB is enough for all possible K.
* time for building index: < 8s (including reading data files, single thread for reading file, and 8 threads for building index).

Command to build it using g++:
```
g++ -std=c++0x -m64 -pthread -O2 main.cpp gene_indexer.cpp prime_generator.cpp murmur3.cpp -o main
```

## Basis of the solution

This problem is asking for mapping of pattern to position, and pattern is one of the segment in a continous string. So its easy to save memory in hash table by storage only position as value. That is, hash(pattern) => position => original gene string.

Position is a limited set: (pattern_length - K) * gene_count, let 4 byte to storage it is good enough.

For hash algorithm, if K is small, just convert it as integer, eg.: K = 10, can be convert to 0 to 4^10; if K is large, use some hash algorithm, I would choose [murmur3](https://github.com/PeterScott/murmur3).

To save more, the original gene string can be storage as binary, but resulting in a little bit slow.

To speed up indexing, we can use multi-threading. I did not implement lock-free hash table here, instead, I just simply devide the data in to several partitions, one per thread. As using one thread for querying, it should be a little slower, but very tinny impact.

# Original Chinese Description

## 2015 年“深圳杯”数学建模夏令营
### B 题：DNA 序列的 k-mer index 问题
这个问题来自 DNA 序列的 k-mer index 问题。

给定一个 DNA 序列，这个系列只含有4个字母ATCG，如 S =“CTGTACTGTAT”。给定一个整数值 k ，从 S 的第一个位置开始，取一连续k个字母的短串，称之为 k-mer （如 k= 5，则此短串为 CTGTA ）， 然后从 S 的第二个位置， 取另一k-mer（如 k= 5，则此短串为 TGTAC ），这样直至 S 的末端，就得一个集合，包含全部k-mer 。 如对序列 S 来说，所有 5-mer 为

> ｛CTGTA，TGTAC，GTACT，TACTG，ACTGT，TGTAT｝

通常这些 k-mer 需一种数据索引方法，可被后面的操作快速访问。例如，对 5-mer 来说，当查询 CTGTA ，通过这种数据索引方法，可返回其在 DNA 序列 S 中的位置为｛1，6｝。

####问题
现在以文件形式给定 100万个 DNA 序列，序列编号为 1 - 1000000，每个基因序列长度为100 。

1. 要求对给定k， 给出并实现一种数据索引方法，可返回任意一个k-mer所在的DNA序列编号和相应序列中出现的位置。每次建立索引，只需支持一个k值即可，不需要支持全部k值。

2. 要求索引一旦建立，查询速度尽量快，所用内存尽量小。

3. 给出建立索引所用的计算复杂度，和空间复杂度分析。

4. 给出使用索引查询的计算复杂度，和空间复杂度分析。

5. 假设内存限制为 8G ，分析所设计索引方法所能支持的最大k值和相应数据查询效率。

6. 按重要性由高到低排列，将依据以下几点，来评价索引方法性能 

  * 索引查询速度

  * 索引内存使用

  * 8G内存下，所能支持的k值范围

  * 建立索引时间
