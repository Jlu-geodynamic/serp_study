/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_diffusion_creep_grainsize_h
#define _aspect_material_model_rheology_diffusion_creep_grainsize_h


#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      /**
       * Data structure for diffusion creep parameters.
       */
      struct DiffusionCreepGrainSizeParameters
      {
        /**
         * The diffusion creep prefactor, activation energy, activation volume
         * and grain size exponent.
         */
        double prefactor_diff;
        double activation_energy_diff;
        double activation_volume_diff;
        double stress_exponent_diff;
        double grain_size_exponent_diff;
		double water_exponent_diff;
        double prefactor_dis;
		double activation_energy_dis;
		double activation_volume_dis;
		double stress_exponent_dis;
		double water_exponent_dis;
		double grain_growth_exponent;
        double grain_growth_rate_constant;
        double grain_growth_activation_energy;
        double grain_growth_activation_volume;
        double boundary_area_change_work_fraction;
		double geometric_constant;
		double grain_boundary_energy;
      };

      template <int dim>
      class DiffusionCreepGrainSize : public ::aspect::SimulatorAccess<dim>
      {
          public:
          /**
           * Constructor.
           */
          DiffusionCreepGrainSize();

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and will be checked against the parsed
           * parameters.
           */
          void
          parse_parameters (ParameterHandler &prm,
		                    const double input_max_x,
		                    const double input_max_y,
							const bool input_add_cratontic,
							const bool input_center_cratontic,
						    const bool input_exclude_mantle_lower,
							const std::vector<double> input_layer_thickness,
						    const std::vector<double> input_cratontic_thicknesses,
							const std::vector<double> input_cratontic_boundary,
                            const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition =
                              std::shared_ptr<std::vector<unsigned int>>());


          /**
           * Compute the creep parameters for the diffusion creep law.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          const DiffusionCreepGrainSizeParameters
          compute_creep_parameters (const unsigned int composition,
                                    const std::vector<double> &phase_function_values = std::vector<double>(),
                                    const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;
          //20240702
		  double
          compute_grainsize(const double second_strain_rate_invariant,
							const double pressure,
							const double temperature,
							const unsigned int composition,
							const double last_grain_size,
							const double posx,
							const double posy,
							const double water_content,
							const double initial_water_content,
							const double yield_stress,
							const double melt_fraction,
							const std::vector<double>& phase_function_values = std::vector<double>(),
							const std::vector<unsigned int>& n_phases_per_composition = std::vector<unsigned int>())const;
          /**
           * Compute the viscosity based on the diffusion creep law.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_viscosity (const double pressure,
                             const double temperature,
                             const unsigned int composition,
                             const double grain_size,
							 const double water_content,
							 const double initial_water_content,
                             const std::vector<double> &phase_function_values = std::vector<double>(),
                             const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;

          /**
            * Compute the strain rate and first stress derivative
            * as a function of stress based on the diffusion creep law.
            * If @p expected_n_phases_per_composition points to a vector of
            * unsigned integers this is considered the number of phase transitions
            * for each compositional field and viscosity will be first computed on
            * each phase and then averaged for each compositional field.
            */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double stress,
                                              const double pressure,
                                              const double temperature,
                                              const double grain_size,
											  const double water_content,
                                              const DiffusionCreepGrainSizeParameters creep_parameters) const;

		  /**
           * Depth of banning grain calculation.
           */
          double banned_depth_meter;
		  /**
           * Grain size of grain calculation banned place.
           */
		  std::vector<double> initial_grain_size;
		  double minimum_grain_size;
		  double max_grain_size;

					
        private:
          /**
           * List of diffusion creep prefactors A.
           */
          std::vector<double> prefactors_diffusion;
          /**
           * List of diffusion creep stress exponents n (usually = 1).
           */
          std::vector<double> stress_exponents_diffusion;
          /**
           * List of diffusion creep grain size exponents m.
           */
          std::vector<double> grain_size_exponents_diffusion;
          /**
           * List of diffusion creep activation energies E.
           */
          std::vector<double> activation_energies_diffusion;
          /**
           * List of diffusion creep activation volumes V.
           */
          std::vector<double> activation_volumes_diffusion;
		  std::vector<double> water_exponents_diffusion;
		  
		  std::vector<double> prefactors_dislocation;
		  std::vector<double> activation_energies_dislocation;
		  std::vector<double> activation_volumes_dislocation;
		  std::vector<double> stress_exponents_dislocation;
		  std::vector<double> water_exponents_dislocation;
		  
		  std::vector<double> grain_growth_exponent;
		  std::vector<double> grain_growth_rate_constant;
		  std::vector<double> grain_growth_activation_energy;
		  std::vector<double> grain_growth_activation_volume;
		  std::vector<double> boundary_area_change_work_fraction;
		  std::vector<double> geometric_constant;
		  std::vector<double> grain_boundary_energy;
		  std::vector<double> flag;
		  std::vector<double> thicknesses;
		  double edot_test;
		  double box_y;
		  //20240702
		  bool add_cratontic;
		  std::vector<double> cratontic_boundary;
		  std::vector<double> cratontic_thicknesses;
		  //20240716
		  bool exclude_mantle_lower;
		  bool center_cratontic;
		  //20241104
		  bool input_from_flank;
		  double box_x;
      };
    }
  }
}
#endif
