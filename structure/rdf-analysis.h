#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "analysis.h"
#include "molecule-analysis.h"

namespace md_analysis {

		class RDFAnalyzer : public AnalysisSet {
			public:
			typedef Analyzer system_t;

				RDFAnalyzer (system_t * t) :
					AnalysisSet (t,
						std::string("RDF Analysis"),
						std::string("rdf.so2-O.wat-H.dat")),
					histo(0.5, 6.0, 0.05) { }

				void Analysis ();
				void DataOutput ();

			protected:

				histogram_utilities::Histogram1D<double> histo;
		};

		class RDFAgent {
			protected:
				Histogram1DAgent histo;
				std::string filename;

			public:

				RDFAgent (std::string fn, const double minimum, const double maximum, const double resolution) :
					histo (fn, minimum, maximum, resolution),
					filename (fn) { }

				void operator() (const double d) { histo(d); }

				std::vector<double> CalcRDF ();
				//double Max () const { return histo.max; }
				//double Min () const { return histo.min; }
				//double Res () const { return histo.res; }
				void SetOutputFilename (std::string fn) { filename = fn; }
				void OutputData ();
		};


		class RDFByDistanceAnalyzer : public molecule_analysis::SuccinicAcidAnalysis {
			protected:
				std::vector<RDFAgent> rdfs;
				double posmin, posmax, posres;

			public:
				RDFByDistanceAnalyzer (Analyzer * t) :
					SuccinicAcidAnalysis (t,
							std::string ("RDF v. Distance analysis"),
							std::string ("")),
					posmin (-12.0), posmax(4.0), posres(2.0) 

			{	// positions of slices of the surface/slab
				rdfs.clear();
				rdfs.resize (8, RDFAgent (std::string (""), 0.5, 6.0, 0.05));	// RDF parameters

				double pos;
				for (int i = 0; i < rdfs.size(); i++) {
					pos = posres * i + posmin;
					std::stringstream sstr;
					sstr.clear();
					std::string filenum;
					filenum.clear();
					sstr << pos;
					filenum = sstr.str();
					std::string filepath (std::string("./alcohol-oxygen-water-oxygen.distance-rdfs/rdf.") + filenum + ".dat");
					rdfs[i].SetOutputFilename (filepath);
				}
			}

				void SuccinicAcidCalculation (alkane::SuccinicAcid *);
				void DataOutput () {
					for (int i = 0; i < rdfs.size(); i++) {
						rdfs[i].OutputData();
					}
				}

				RDFAgent * FindRDFAgent (const double pos);

		};

} // namespace md_analysis
#endif
