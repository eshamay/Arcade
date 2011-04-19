#ifndef DATAOUTPUT_H_
#define DATAOUTPUT_H_

#include "patterns.h"
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <cmath>

namespace md_analysis {

	class StatusUpdater : public patterns::observer::observer {
		protected:
			virtual void _updateStatus () = 0;
			double _frequency;
			double _count;
			double _maxcount;

		public:

			StatusUpdater () : _count(0) { }
			StatusUpdater (const double frequency, const double maxcount, const double startingcount = 0) 
				: _frequency(frequency), _count(startingcount), _maxcount(maxcount) { }

			virtual void Set (const double frequency, const double maxcount, const double startingcount = 0) {
				_frequency = frequency;
				_maxcount = maxcount;
				_count = startingcount;
			}
			
			// every time the updater is called the count is updated, and then output is performed based on the specific frequency supplied
			virtual void notify () {
				_count++;
				this->_updateStatus ();
			}
	};


	class PercentProgressBar : public StatusUpdater {

		protected:
			virtual void _updateStatus ();
			void progressbar (int percent);
			std::stringstream bars;
			int x;
			std::string slash[4];

		public:
			PercentProgressBar () : StatusUpdater(), x(0) { 
				slash[0] = "\\";
				slash[1] = "-";
				slash[2] = "/";
				slash[3] = "|";
				bars << "|";
			}

			PercentProgressBar (const double frequency, const double maxcount, const double startingcount = 0) 
				: StatusUpdater (frequency, maxcount, startingcount), x(0) { 
					slash[0] = "\\";
					slash[1] = "-";
					slash[2] = "/";
					slash[3] = "|";
					bars << "|";
				}

			virtual ~PercentProgressBar () { }

	};


	class StarStatusBarUpdater : public StatusUpdater {
		protected:
			virtual void _updateStatus ();

		public:
			StarStatusBarUpdater () : StatusUpdater() { }
			virtual ~StarStatusBarUpdater () { }

			StarStatusBarUpdater (const double frequency, const double maxcount, const double startingcount = 0) 
				: StatusUpdater (frequency, maxcount, startingcount) { }

	};	// Star status bar updater

}	// namespace md_analysis


#endif
