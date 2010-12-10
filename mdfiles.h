#ifndef MD_FILES_H_
#define MD_FILES_H_

#include <cstdio>
#include <string>
#include <vector>
#include "vecr.h"

namespace md_files {

	class CoordinateFile {

		public:

			CoordinateFile (const std::string path, int const c_size) 
				:
					_file ((FILE *)NULL),
					_path(path),
					_size(c_size),
					_coords (_size*3, 0.0),
					_frame(0),
					_eof(false)	{

						_file = fopen64 (path.c_str(), "r");
						if (_file == (FILE *)NULL)
						{
							printf ("Couldn't load the Coordinate file %s\n", path.c_str());
							exit(1);
						}

						// and then map those elements to vectors
						_vectors.clear();
						for (int i = 0; i < _size; i++)
							_vectors.push_back (Eigen::Map<VecR> (&_coords[3*i]));
					}

			virtual ~CoordinateFile () {
				fclose (_file);
			}

			char * ReadLine () { 
				fgets (_line, 1000, _file);
				return _line;
			}

			char * Line () {
				return _line;
			}

			virtual void LoadNext () = 0;

			// retrieves coordinates as VecR (3-element vectors)
			const coord_t& Coordinate (const int index) const { return _vectors[index]; }
			const coord_t& operator() (const int index) const { return _vectors[index]; }

			coord_it begin () const { return _vectors.begin(); }
			coord_it end () const { return _vectors.end(); }

			// retrieves the coordinate array
			const std::vector<double>& GetArray () const { return _coords; }

			int size () 	const { return _size; }

			bool eof () 	const { return _eof; }
			int Frame () 	const { return _frame; }

		protected:
			FILE *_file;				// the XYZ file listing all the atom coordinates
			std::string _path;

			int				_size;		// number of coordinates to parse in each frame

			std::vector<double>								_coords;				// array of atomic coordinates
			coord_set_t												_vectors;				// set of vectors representing positions

			char _line[1000];
			int 			_frame;		// The current frame (number of timesteps processed)
			bool			_eof;		// end of file marker for the coord file

	};	// Coordinate file


}	// namespace md_files

#endif
