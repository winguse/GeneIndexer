#include "gene_indexer.h"


gene_indexer::gene_indexer(const uint32_t pattern_length, const uint32_t gene_count, const uint32_t gene_length, const bool compress_data, const uint32_t hash_table_size)
{
	m_pattern_length    = pattern_length;
	m_gene_count        = gene_count;
	m_gene_length       = gene_length;
	m_compress_data     = compress_data;
	m_hash_table_size   = hash_table_size;

	uint32_t max_address_count = (m_gene_length - m_pattern_length + 1) * m_gene_count;

	if(m_hash_table_size == 0)
	{
		prime_generator *pg = new prime_generator(max_address_count << 1);
		m_hash_table_size = max_address_count * 1.5;
		while(!pg->is_prime(m_hash_table_size))m_hash_table_size++;
		delete pg;
		pg = NULL;
	}

	uint32_t tmp = m_hash_table_size;
	for(uint32_t i = 0; i < m_pattern_length; i++)
	{
		tmp >>= 2; // equal to: tmp /= 4
	}
	m_direct_hash_mode = tmp > 0; // m_hash_table_size > power(4, m_pattern_length);

	if(m_direct_hash_mode)
	{
		m_hash_table_size = 1 << (m_pattern_length << 1);
#ifdef DEBUG
		printf("direct hash mode is enabled, m_hash_table_size = %d\n");
#endif // DEBUG
	}

	m_hash_table = new uint32_t[m_hash_table_size]; // this can save more memory if the address is narrow, eg, short. or if address is more than 2^32, this should update
	memset(m_hash_table, -1, m_hash_table_size << 2);

	int max_data_length = gene_count * gene_length;
	if(m_compress_data)
	{
		max_data_length >>= 2;
		max_data_length++;
	}
	m_data = new uint8_t[max_data_length];
	memset(m_data, 0, max_data_length);
	m_data_length = 0;

	if(m_compress_data)
	{
		tmp = m_pattern_length >> 2; // pattern length / 4 == bytes length of a pattern
		m_pattern_actual_byte_count = tmp;
		if(tmp & 3) m_pattern_actual_byte_count++;
		m_pattern_buffer_size = tmp;
	}
	else
	{
		m_pattern_actual_byte_count = m_pattern_length;
		m_pattern_buffer_size       = m_pattern_length;
	}

	m_pattern_buffer_size >>= 3;
	if(tmp & 7) m_pattern_buffer_size++;
	m_pattern_buffer_size <<= 3; // at least 8 bytes as we use uint64_t for optimization

	memset(m_gene_2_byte, -1, sizeof(m_gene_2_byte));
	m_gene_2_byte['A'] = 0b00;
	m_gene_2_byte['T'] = 0b01;
	m_gene_2_byte['C'] = 0b10;
	m_gene_2_byte['G'] = 0b11;

	m_duplicate_gene_pattern_count  = 0;
	m_re_address_count              = 0;
	m_max_re_address_times          = 0;
	m_once_hit_count                = 0;
}

gene_indexer::~gene_indexer()
{
	delete m_hash_table;
	delete m_data;
}


void gene_indexer::add(const char* gene_string)
{
	for(int i = 0; gene_string[i]; i++)
	{
		uint8_t x = m_gene_2_byte[(uint8_t)gene_string[i]];
		if(x == 0xFF)continue;
		if(m_compress_data)
		{
			m_data[m_data_length >> 2] |= (x << ((m_data_length & 3) << 1));
		}
		else
		{
			m_data[m_data_length] = gene_string[i];
		}
		m_data_length++;
	}
}
bool gene_indexer::is_full() const
{
	return m_data_length == m_gene_count * m_gene_length;
}

uint32_t gene_indexer::get_hash_table_slot(const void* gene_buff, const uint32_t logical_address)
{
	uint32_t hashed_address;
	char gene_buff_2[DEFAULT_BUFFER_SIZE];
	if(!m_direct_hash_mode)
	{
		uint64_t hashed_address_128[2];
		MurmurHash3_x64_128(gene_buff, m_pattern_actual_byte_count, 0, hashed_address_128);
		hashed_address = hashed_address_128[0] % m_hash_table_size;

		uint32_t current_re_address_times = 0;
		while(true)
		{
			if(m_hash_table[hashed_address] == UINT_32_MAX)
			{
				break;
			}

			export_gene_pattern(m_hash_table[hashed_address], gene_buff_2);
			if(are_equal(gene_buff, gene_buff_2))
			{
				if(logical_address != UINT_32_MAX) m_duplicate_gene_pattern_count++;
				break;
			}
			if(logical_address != UINT_32_MAX)
			{
				current_re_address_times++;
				m_re_address_count++;
			}
			hashed_address++;
			if(hashed_address == m_hash_table_size) hashed_address = 0;
		}
		if(logical_address != UINT_32_MAX)
		{
			if(current_re_address_times == 0) m_once_hit_count++;
			if(current_re_address_times > m_max_re_address_times) m_max_re_address_times = current_re_address_times;
		}

	}
	else
	{
		if(m_compress_data)
		{
			hashed_address = *((uint32_t*)gene_buff);
		}
		else
		{
			to_binary_gene((char*)gene_buff, gene_buff_2);
			hashed_address = *((uint32_t*)gene_buff_2);
		}
	}

	if(logical_address != UINT_32_MAX && m_hash_table[hashed_address] == UINT_32_MAX)
	{
		m_hash_table[hashed_address] = logical_address;
	}
	return m_hash_table[hashed_address];
}

void gene_indexer::export_gene_pattern(const uint32_t position, void* gene_buff) const
{
	memset(gene_buff, 0, m_pattern_buffer_size);
	__export_gene_pattern(position, m_pattern_length, gene_buff);
}

void gene_indexer::__export_gene_pattern(const uint32_t position, const uint32_t length, void* gene_buff) const
{
	char *out = (char*) gene_buff;
	if(m_compress_data)
	{
		for(uint32_t i = 0; i < length; i++)
		{
			uint32_t byte_idx = (position + i) >> 2;
			uint32_t slot_idx = (position + i) & 3;
			uint32_t bit_idx = slot_idx << 1;
			out[i >> 2] |= ((m_data[byte_idx] & (3 << bit_idx)) >> bit_idx) << ((i & 3) << 1);
		}
	}
	else
	{
		for(uint32_t i = 0; i < length; i++)
		{
			out[i] = m_data[position + i];
		}
	}
}

bool gene_indexer::are_equal(const void *gene_buff_a, const void *gene_buff_b, uint32_t gene_length) const
{
	if(gene_length == UINT_32_MAX) gene_length = m_pattern_length;
	uint32_t max_idx = gene_length >> 5; // using uint64_t, 8 bytes, 32 gene
	if(gene_length & 31) max_idx++;
	uint64_t* a =  (uint64_t*)gene_buff_a;
	uint64_t* b =  (uint64_t*)gene_buff_b;
	for(uint32_t idx = 0; idx < max_idx; idx++)
	{
		if(*a != *b)
		{
			return false;
		}
		a++;
		b++;
	}
	return true;
}

void gene_indexer::build_index()
{
	uint32_t position, logical_position;
	char gene_buff[DEFAULT_BUFFER_SIZE];
	for(uint32_t i = 0; i < m_gene_count; i++)
	{
		for(uint32_t j = 0; j < m_gene_length - m_pattern_length + 1; j++)
		{
			position = i * m_gene_length + j;
			// logic address can be the same as position, redusing its size can resulting in lower memory
			// eg.: logical_position = i * (m_gene_length - m_pattern_length + 1) + j;
			logical_position = position;
			export_gene_pattern(position, gene_buff);
			get_hash_table_slot(gene_buff, logical_position);
		}
	}
}

uint32_t gene_indexer::search(const char *gene_pattern_string)
{
	if(m_compress_data)
	{
		char gene_buff[DEFAULT_BUFFER_SIZE];
		to_binary_gene(gene_pattern_string, gene_buff);
		return get_hash_table_slot(gene_buff);
	}
	return get_hash_table_slot(gene_pattern_string);
}


void gene_indexer::to_binary_gene(const char *gene_pattern_string, char *gene_buff) const
{
	memset(gene_buff, 0, m_pattern_buffer_size);
	for(int i = 0; gene_pattern_string[i]; i++)
	{
		uint8_t x = m_gene_2_byte[(uint8_t)gene_pattern_string[i]];
		if(x == 0xFF)continue;
		gene_buff[i >> 2] |= (x << ((i & 3) << 1));
	}
}
