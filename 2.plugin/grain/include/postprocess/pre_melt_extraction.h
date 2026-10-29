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


#ifndef _aspect_postprocess_pre_melt_extraction_h
#define _aspect_postprocess_pre_melt_extraction_h

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
    class PreMeltExtraction : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;
		
		private:
		
		double depth_max;
		double y_add;
		unsigned int melt_fraction_refinement;
		double max_x;
		unsigned int repetitions_x;
		double max_y;
		unsigned int repetitions_y;
		unsigned int min_container_amount_y;
		unsigned int melt_container_amount_x;
		unsigned int melt_container_amount_y;
		unsigned int flag_amount;
		double melt_limit;
		double melt_factor;
		bool include_mantle_melting;
		bool include_crustal_melting;
		std::string convey_flag_mantle_melting_path;
		std::string container_mantle_melting_path;
		std::string convey_flag_crustal_melting_path;
		std::string container_crustal_melting_path;
		std::string melt_info_path;
		//20250428
		bool include_carbonate_melting;
		std::string container_carbon_melting_path;
		//20251104
		bool dike_mode;
		
		
		
					   
    };
  }
}


#endif
