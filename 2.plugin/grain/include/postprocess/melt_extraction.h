/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_postprocess_melt_extraction_h
#define _aspect_postprocess_melt_extraction_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  using namespace dealii;
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the compositional
     * fields, if any.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class MeltExtraction : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some temperature statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
		
		double
		read_from_txt(const std::string path,
	                  const unsigned int mpi_communicator)const;

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;
		
		private:
		
		unsigned int melt_fraction_refinement;
		double max_x;
		unsigned int repetitions_x;
		std::string melt_info_path;
		bool vertical_trans;
		bool include_mantle_melting;
		bool include_crustal_melting;
		//20251104
		bool include_carbonate_melting;
		bool dike_mode;
		double dike_bottom_temperature;
		double max_y;
		unsigned int repetitions_y;
		double max_refinement_length_x;
		double max_refinement_length_y;
		//20251201
		bool constant_dike_length;
		bool constant_dike_x;
		double dike_top_temperature;
		double dike_width;
		double dike_length;
					   
    };
  }
}


#endif
