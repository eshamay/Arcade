//#include <boost/thread/thread.hpp>
#include "pthread.h"

namespace threads {

	/* id = process id
	 * p  = number of processes
	 * n  = size of the problem (e.g. number of indices of an array, etc.)
	 *
	 * numbering starts at 0 for index
	 */
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
