#ifndef PRIME_GENERATOR_H
#define PRIME_GENERATOR_H
#include <cmath>

class prime_generator
{
public:
	prime_generator(int max_number);
	bool is_prime(int x);
	virtual ~prime_generator();
protected:
private:
	int     m_max_number;
	bool*   m_is_prime;
};

#endif // PRIME_GENERATOR_H
