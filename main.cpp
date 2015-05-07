#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <thread>
#if defined (__linux__)
#include <sys/time.h>
#else
#include <ctime>
#endif

#include "gene_indexer.h"

using namespace std;

#if defined (__linux__)
struct timeval start_time, last_time;
#else
clock_t start_time, last_time;
#endif

void log(const char *format, ...)
{
#if defined (__linux__)
	struct timeval current_time;
	gettimeofday(&current_time, NULL);
#else
	clock_t current_time = clock();
#endif
	printf(
		"%8.3f %8.3f | ",
#if defined (__linux__)
		current_time.tv_sec - start_time.tv_sec + (current_time.tv_usec - start_time.tv_usec) / 1e6 ,
		current_time.tv_sec - last_time.tv_sec + (current_time.tv_usec - last_time.tv_usec) / 1e6
#else
		(current_time - start_time) / (double) CLOCKS_PER_SEC, (current_time - last_time) / (double) CLOCKS_PER_SEC
#endif
	);
	va_list args;
	va_start (args, format);
	vprintf (format, args);
	va_end (args);

#if defined (__linux__)
	gettimeofday(&last_time, NULL);
#else
	last_time = clock();
#endif
}

void print_gene(const void *gene, const int gene_length)
{
	const char* gene_chars = "ATCG";
	char *_data = (char*)gene;
	for(int i = 0; i < gene_length; i++)
	{
		int bit_idx = ((i & 3) << 1);
		putchar(gene_chars[(_data[i >> 2] & (3 << bit_idx)) >> bit_idx]);
	}
}

int main()
{
#if defined (__linux__)
	gettimeofday(&start_time, NULL);
	gettimeofday(&last_time, NULL);
#else
	start_time = last_time = clock();
#endif

	log("start\n");

	int         k               = 30;
	int         gene_count      = 1000000;
	int         gene_length     = 100;
	int         thread_count    = 8;
	bool        compress_mode   = false;
	int         hash_table_size = 13312517; // this value is base on the config above

	const char  *file_name      = "100w";

	//scanf("%d", &thread_count);

	int gene_per_indexer = gene_count / thread_count;

	gene_indexer **indexers = new gene_indexer*[thread_count];


	for(int i = 0; i < thread_count - 1; i++)
	{
		indexers[i] = new gene_indexer(k, gene_per_indexer, gene_length, compress_mode, hash_table_size);
	}
	indexers[thread_count - 1] = new gene_indexer(k, gene_count - gene_per_indexer * (thread_count - 1), gene_length, compress_mode, hash_table_size);

	log("m_hash_table_size = %d\n", indexers[thread_count - 1]->m_hash_table_size);

	log("allocate indexer memory\n");

	FILE *file = fopen(file_name, "r");
	char buffer[1024];
	int current_idx = 0;
	////while(~fscanf(file, "%s", buffer)) // slow
	while(fgets(buffer, 1024, file))
	{
		if(buffer[0] == '>') continue;
		if(buffer[0] == '\0') break;
		if(current_idx >= thread_count)
		{
			log("data issue!");
			return 1;
		}
		indexers[current_idx]->add(buffer);
		if(indexers[current_idx]->is_full())
		{
			log("indexers[%d] is full.\n", current_idx);
			current_idx++;
		}
	}
	fclose(file);

	log("%d gene was read.\n", gene_count);

	if(thread_count == 1)
	{
		indexers[0]->build_index();
	}
	else
	{
		thread **threads = new thread*[thread_count];
		for(int i = 0; i < thread_count; i++)
		{
			threads[i] = new thread(&gene_indexer::build_index, indexers[i]);
		}
		for(int i = 0; i < thread_count; i++)
		{
			threads[i]->join();
			delete threads[i];
		}
		delete threads;
	}

	for(int i = 0; i < thread_count; i++)
	{
		log("indexers[%d] build is done, input to search.\n", i);
		log("m_duplicate_gene_pattern_count  = %d\n", indexers[i]->m_duplicate_gene_pattern_count);
		log("m_re_address_count              = %d\n", indexers[i]->m_re_address_count);
		log("m_max_re_address_times          = %d\n", indexers[i]->m_max_re_address_times);
		log("m_once_hit_count                = %d\n", indexers[i]->m_once_hit_count);
		log("m_direct_hash_mode              = %d\n", indexers[i]->m_direct_hash_mode);
	}

	while(true)
	{
		log("please input: ");
		if(!gets(buffer) || buffer[0] == '\0') break;
		log("got    input: %s\n", buffer);
		uint32_t pos = 0xFFFFFFFF;
		for(int i = 0; i < thread_count; i++)
		{
			pos = indexers[i]->search(buffer);
			if(pos != 0xFFFFFFFF)
			{
				pos += i * gene_per_indexer * gene_length;
				log("[%02d] it is at %d, %d\n", i, pos / gene_length, pos % gene_length);
				break;
			}
		}
		if(pos == 0xFFFFFFFF)
		{
			log("it was not found.\n");
		}
	}

	for(int i = 0; i < thread_count; i++)
	{
		delete indexers[i];
	}

	delete indexers;
	log("program terminated.\n");
	return 0;
}
