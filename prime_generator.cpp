#include "prime_generator.h"

prime_generator::prime_generator(int max_number)
{
	m_max_number = max_number;
	m_is_prime = new bool[m_max_number];
	for(int i = 0; i < m_max_number; i++)
	{
		m_is_prime[i] = true;
	}

	m_is_prime[0] = m_is_prime[1] = false;

	int max_start = sqrt(m_max_number);
	for(int i = 2; i <= max_start; i++)
	{
		if(!m_is_prime[i]) continue;
		for(int j = i + i; j < m_max_number; j += i)
			m_is_prime[j] = false;
	}

}

bool prime_generator::is_prime(int x)
{
	return m_is_prime[x];
}

prime_generator::~prime_generator()
{
	delete m_is_prime;
}
