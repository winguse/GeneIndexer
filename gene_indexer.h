#ifndef GENE_INDEXER_H
#define GENE_INDEXER_H
#include "prime_generator.h"
#include "murmur3.h"
#include <cstring>
#include <cstdint>


const uint32_t DEFAULT_BUFFER_SIZE  = 1024; // should be greater than m_pattern_buffer_size
const uint32_t UINT_32_MAX          = 0xFFFFFFFF;

class gene_indexer
{
public:
	gene_indexer(const uint32_t pattern_length, const uint32_t gene_count, const uint32_t gene_length, const bool compress_data = true, const uint32_t hash_table_size = 0);
	void add(const char* gene_string);
	void export_gene_pattern(const uint32_t position, void* gene_buff) const;
	void __export_gene_pattern(const uint32_t position, const uint32_t length, void* gene_buff) const;
	void build_index();
	uint32_t search(const char *gene_pattern_string);
	bool is_full() const;
	virtual ~gene_indexer();

	//diagnostic
	uint32_t    m_duplicate_gene_pattern_count;
	uint32_t    m_re_address_count;
	uint32_t    m_max_re_address_times;
	uint32_t    m_once_hit_count;

	bool        m_direct_hash_mode;
	bool        m_compress_data;

	uint32_t    m_hash_table_size;
protected:
	uint32_t get_hash_table_slot(const void* gene_buff, const uint32_t logical_address = UINT_32_MAX);
	bool are_equal(const void *gene_buff_a, const void *gene_buff_b, uint32_t gene_length = UINT_32_MAX) const;
	void to_binary_gene(const char *gene_pattern_string, char* gene_buff) const;
private:
	uint32_t    m_pattern_length;
	uint32_t    m_gene_count;
	uint32_t    m_gene_length;

	uint8_t*    m_data;
	uint32_t    m_data_length;
	uint32_t*   m_hash_table;

	uint32_t    m_pattern_buffer_size;
	uint32_t    m_pattern_actual_byte_count;

	uint8_t     m_gene_2_byte[128];
};

#endif // GENE_INDEXER_H
