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


#ifndef _aspect_postprocess_melt_convey_1_h
#define _aspect_postprocess_melt_convey_1_h

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
    class MeltConveyOne : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some temperature statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
		
		
		void
		write_to_container(const std::string path,
	                      const unsigned int crood,
	                      const unsigned int mpi_communicator,
				          const double write_information)const;
		
		double
		read_from_container(const std::string path,
	                       const unsigned int crood,
	                       const unsigned int mpi_communicator)const;
		
		double
		read_from_txt(const std::string path,
	                  const unsigned int mpi_communicator)const;
		
		void
		convey_flag_update(const std::string path,
	                       const unsigned int crood_x,
					       const unsigned int crood_y_in_flag_container,
					       const unsigned int mpi_communicator,
					       const unsigned int flag)const;
						   
        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;
		
		private:
		
		double depth_max;
		unsigned int melt_fraction_refinement;
		double max_x;
		unsigned int repetitions_x;
		double max_y;
		unsigned int repetitions_y;
		unsigned int min_container_amount_y;
		std::string convey_flag_mantle_melting_path;
		std::string container_mantle_melting_path;
		std::string convey_flag_crustal_melting_path;
		std::string container_crustal_melting_path;
		std::string melt_info_path;
		double container_limit;
		
		bool include_mantle_melting;
		bool include_crustal_melting;
		//20250428
		bool include_carbonate_melting;
		bool vertical_trans;
		std::string carbon_container_path;
		//202602023
		bool salt_mode;
					   
    };
  }
}


#endif
