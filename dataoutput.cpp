#include "dataoutput.h"

namespace md_analysis {

	void StarStatusBarUpdater::_updateStatus () {
		if (!fmod(_count, this->_frequency * 10))
			std::cout << std::endl << _count << "/" << this->_maxcount << " ) ";
		if (!fmod(_count, this->_frequency)) {
			std::cout << "*";
		}

		fflush (stdout);
	}


	void PercentProgressBar::_updateStatus () {
		int percent = (int)(_count / _maxcount * 100.0);
		this->progressbar (percent);
	}

	void PercentProgressBar::progressbar (int percent) {
		std::cout << "\r"; // carriage return back to beginning of line
		std::cout << bars.str() << " " << slash[x] << " " << _count << " / " << _maxcount << " == " << percent << "%"; // print the bars and percentage
		x++; // increment to make the slash appear to rotate
		if(x == 4)
			x = 0; // reset slash animation
		//if (!fmod(_count, this->_frequency * 10))
			//bars << slash[x];
	}

}// namespace md analysis
