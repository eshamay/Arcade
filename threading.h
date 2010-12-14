#include <boost/thread/thread.hpp>

namespace threads {

	int block_low (const int& id, const int& p, const int& n) {
		return id*n/p;
	}

	int block_high (const int& id, const int& p, const int& n) {
		return block_low (id+1,p,n) - 1;
	}

	int block_size (const int& id, const int& p, const int& n) {
		return block_low (id+1,p,n) - block_low(id,p,n);
	}

	int block_owner (int& id, int& p, int& n) {
		return (p*(id+1)-1)/n;
	}

}	// namespace threads
