/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_visco_plastic_grainsize_h
#define _aspect_material_model_rheology_visco_plastic_grainsize_h
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/rheology/diffusion_creep_grainsize.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/rheology/dislocation_creep_grainsize.h>
#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/rheology/strain_dependent.h>
#include <aspect/material_model/rheology/frank_kamenetskii.h>
#include <aspect/material_model/rheology/peierls_creep.h>
#include <aspect/material_model/rheology/constant_viscosity_prefactors.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/simulator_access.h>

#include<deal.II/fe/component_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for the plastic parameters weakened (or hardened)
     * by strain to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class PlasticAdditionalOutputsGrainSize : public NamedAdditionalMaterialOutputs<dim>
    {
    public:
        PlasticAdditionalOutputsGrainSize(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Cohesions at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> cohesions;

        /**
         * Internal angles of friction at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> friction_angles;

        /**
         * The area where the viscous stress exceeds the plastic yield stress,
         * and viscosity is rescaled back to the yield envelope.
         */
        std::vector<double> yielding;

        //两种机制计算得到的粘度比例
        std::vector<double> creep_visco_ratio;

        //软流圈
        std::vector<double> volume_0;
		
		std::vector<double> eta_dislocation;
		
		std::vector<double> eta_diffusion;
		
		std::vector<double> eta_for_viscous;
		
		std::vector<double> plastic_stress;
		
		std::vector<double> viscous_stress;
		
		std::vector<double> yield_strength;
		
    };

    /**
       * A data structure with the output of calculate_isostrain_viscosities.
       */
    struct IsostrainViscosities
    {
      /**
       * The composition viscosity.
       */
      std::vector<double> composition_viscosities;

      /**
       * The composition yielding.
       */
      std::vector<bool> composition_yielding;
    };

    namespace Rheology
    {
	  //20251110
	  //新增岩墙相关：岩墙模式下累积应变类型
	  enum DikemodeWeakeningMechanism
      {
        dike_none,
        dike_plastic_weakening_with_plastic_strain_only,
        dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain,
		dike_plastic_weakening_with_plastic_strain_and_viscous_weakening_with_grain_size
      };

      template <int dim>
      class ViscoPlasticGrainSize : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ViscoPlasticGrainSize();
		  
		  /* std::vector<double> 
		  melt_convey_function(const double posx,
		                       const double points,
							   const std::vector<double> &composition,
							   const std::vector<double> &raw_volume_fractions)const; */

          std::vector<double>
	      fractions_modify(const std::vector<double> &temp_volume_fractions) const;
		  
		  //20251104
		  double
		  read_from_txt(const std::string path,
	                    const unsigned int mpi_communicator)const;
	      bool
	      dike_depend(const double posx,
	                  const double posy,
	                  const unsigned int mpi_communicator)const;
		  //20251110
          std::array<double, 3>
          dikemode_compute_strain_weakening_factors(const unsigned int j,
                                                    const std::vector<double> &composition) const;
		
          void
          dikemode_fill_reaction_outputs(const MaterialModel::MaterialModelInputs<dim> &in,
                                         const int i,
                                         const double min_strain_rate,
                                         const bool plastic_yielding,
										 const bool in_dike,
									     const SymmetricTensor<2,dim> strain_rate_dev,
                                         MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * This function calculates viscosities assuming that all the compositional fields
           * experience the same strain rate (isostrain).
           */
          IsostrainViscosities
          calculate_isostrain_viscosities ( const MaterialModel::MaterialModelInputs<dim> &in,
                                            const unsigned int i,
                                            const std::vector<double> &volume_fractions,
                                            const std::vector<double> &phase_function_values = std::vector<double>(),
                                            const std::vector<unsigned int> &n_phases_per_composition =
                                              std::vector<unsigned int>()) const;

          /**
           * A function that fills the viscosity derivatives in the
           * MaterialModelOutputs object that is handed over, if they exist.
           * Does nothing otherwise.
           */
          void compute_viscosity_derivatives(const unsigned int point_index,
                                             const std::vector<double> &volume_fractions,
                                             const std::vector<double> &composition_viscosities,
                                             const MaterialModel::MaterialModelInputs<dim> &in,
                                                   MaterialModel::MaterialModelOutputs<dim>& out,
                                             const std::vector<double> &phase_function_values = std::vector<double>(),
                                             const std::vector<unsigned int> &n_phases_per_composition =
                                               std::vector<unsigned int>()) const;

          /**
           * A function that returns a ComponentMask that represents all compositional
           * fields that should be considered 'volumetric', that is representing a
           * physical proportion of the material, e.g. volume fraction of peridotite
           * (as opposed to non-volumetric quantities like the amount of finite-strain).
           */
          ComponentMask get_volumetric_composition_mask() const;


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
		   
		  //20251104
		  //20251201 新增 constant_dike_length
          void
          parse_parameters (ParameterHandler &prm,
		                    const double max_x,
		                    const double max_y,
							const bool add_cratontic,
							const bool center_cratontic,
							const bool exclude_mantle_lower,
							const bool input_dike_mode,
							const bool input_constant_M,
							const bool input_constant_dike_length,
							const double input_M_for_dike,
							const double input_dike_width,
							const double input_dike_length,
							const double spreading_rate,
							const std::vector<double> layer_thickness,
							const std::vector<double> cratontic_thicknesses,
							const std::vector<double> cratontic_boundary,
                            const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition =
                              std::shared_ptr<std::vector<unsigned int>>());

		  
          void fill_grainsize_outputs(const MaterialModel::MaterialModelInputs<dim>& in,
                                     MaterialModel::MaterialModelOutputs<dim>& out,
                                     const unsigned int i,
									 const bool plastic_yielding,
                                     const std::vector<double>& volume_fractions,
                                     const std::vector<double>& phase_function_values = std::vector<double>(),
                                     const std::vector<unsigned int>& n_phases_per_composition =
                                     std::vector<unsigned int>()) const;
          /**
           * Create the additional material model outputs object that contains the
           * plastic outputs.
           */
          void
          create_plastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * A function that fills the plastic additional output in the
           * MaterialModelOutputs object that is handed over, if it exists.
           * Does nothing otherwise.
           */
          void fill_plastic_outputs (const unsigned int point_index,
                                     const std::vector<double> &volume_fractions,
                                     const bool plastic_yielding,
                                     const MaterialModel::MaterialModelInputs<dim> &in,
                                     MaterialModel::MaterialModelOutputs<dim> &out,
									 const std::vector<double>& phase_function_values = std::vector<double>(),
                                     const std::vector<unsigned int>& n_phases_per_composition =
                                     std::vector<unsigned int>()) const;

          /**
           * Reference viscosity used by material models using this rheology.
           */
          double ref_visc;

          /**
           * Minimum strain rate used to stabilize the strain rate dependent rheology.
           */
          double min_strain_rate;
		   
		  /**
           * Ratio of fractions modification.
           */
		  double fraction_ratio;

          /**
           * Enumeration for selecting which viscosity averaging scheme to use.
           */
          MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

          /**
           * Object for computing the strain dependence of the rheology model.
           */
          Rheology::StrainDependent<dim> strain_rheology;

          /**
           * Object for computing viscoelastic viscosities and stresses.
           */
          Rheology::Elasticity<dim> elastic_rheology;

          /**
           * Whether to include viscoelasticity in the constitutive formulation.
           */
          bool use_elasticity;
		  
		  //20241104
		  double box_x;
		  bool input_from_flank;


        private:

          /**
           * Reference strain rate for the first non-linear iteration
           * in the first time step.
           */
          double ref_strain_rate;

          /**
           * Minimum and maximum viscosities used to improve the
           * stability of the rheology model.
           */
          double min_visc;
          double max_visc;
		  
		  std::vector<double> init_water_content;

          /**
           * Enumeration for selecting which type of viscous flow law to use.
           * Select between diffusion, dislocation, frank_kamenetskii or composite.
           */
          enum ViscosityScheme
          {
            diffusion,
            dislocation,
            frank_kamenetskii,
            composite
          } viscous_flow_law;

          /**
           * Enumeration for selecting which type of yield mechanism to use.
           * Select between Drucker Prager and stress limiter.
           */
          enum YieldScheme
          {
            stress_limiter,
            drucker_prager
          } yield_mechanism;

          /**
           * Whether to allow negative pressures to be used in the computation
           * of plastic yield stresses and viscosities. If false, the minimum
           * pressure in the plasticity formulation will be set to zero.
           */
          bool allow_negative_pressures_in_plasticity;

          /**
           * List of exponents controlling the behaviour of the stress limiter
           * yielding mechanism.
           */
          std::vector<double> exponents_stress_limiter;

          /**
           * Temperature gradient added to temperature used in the flow law.
           */
          double adiabatic_temperature_gradient_for_viscosity;

          /**
           * Objects for computing viscous creep viscosities.
           */
          Rheology::DiffusionCreepGrainSize<dim> diffusion_creep_grainsize;
          Rheology::DislocationCreepGrainSize<dim> dislocation_creep_grainsize;
          std::unique_ptr<Rheology::FrankKamenetskii<dim> > frank_kamenetskii_rheology;

          /**
           * Whether to include Peierls creep in the constitutive formulation.
           */
          bool use_peierls_creep;

          /**
           * Object for computing Peierls creep viscosities.
           */
          std::unique_ptr<Rheology::PeierlsCreep<dim> > peierls_creep;

          /**
           * Object for computing the viscosity multiplied by a constant prefactor.
           * This multiplication step is done just prior to calculating the effective
           * viscoelastic viscosity or plastic viscosity.
           */
          Rheology::ConstantViscosityPrefactors<dim> constant_viscosity_prefactors;

          /*
           * Object for computing plastic stresses, viscosities, and additional outputs
           */
          Rheology::DruckerPrager<dim> drucker_prager_plasticity;

          /*
           * Input parameters for the drucker prager plasticity.
           */
          Rheology::DruckerPragerParameters drucker_prager_parameters;
		  
		  //20250122
		  double switch_min_visc;
		  double switch_time_in_year;
		  
		  //20250916
		  double input_total_force;
		  
		  //20251104
		  std::string melt_info_path;
		  double box_y;
		  unsigned int repetitions_x;
		  unsigned int repetitions_y;
		  unsigned int melt_fraction_refinement;
		  bool dike_mode;
		  bool constant_M;
		  double M_for_dike;
		  double dike_width;
		  double dike_length;
		  double spreading_rate;
		  //20251110
		  DikemodeWeakeningMechanism dike_weakening_mechanism;
		  //20251201
		  bool constant_dike_length;
		  //20260320
		  double serp_weakening;
		  

      };
    }
  }
}
#endif
