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


#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/postprocess/melt_convey_2.h>
#include </fs2/home/liuzhonglan/wy/lib_extra/melt20260814/grain/include/postprocess/pre_melt_extraction.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <iostream>
#include <dirent.h>
#include <sys/stat.h>


namespace aspect
{
  namespace Postprocess
  {

	  
	template <int dim>
	void
	MeltConveyTwo<dim>::
	write_to_container(const std::string path,
	                  const unsigned int crood,
	                  const unsigned int mpi_communicator,
				      const double write_information)const
	{
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::stringstream get_name;
		    get_name << crood;
		    std::string temp_txt_name = get_name.str();
		    std::ofstream write((path + temp_txt_name + ".txt").c_str());
			while(!write)
			{
				write.close();
				sleep(1.);
				std::ofstream write((path + temp_txt_name + ".txt").c_str());
			}
		    write << write_information;
		    write.close();
		}
	}
	
	template <int dim>
	double
	MeltConveyTwo<dim>::
	read_from_container(const std::string path,
	                   const unsigned int crood,
	                   const unsigned int mpi_communicator)const
	{
		double temp_input = 0;
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::stringstream get_name;
		    get_name << crood;
		    std::string temp_txt_name = get_name.str();
			std::ifstream data_from_txt((path + temp_txt_name + ".txt").c_str());
			std::string line;
			getline(data_from_txt, line);
			while(line == "")
			{
				data_from_txt.close();
				sleep(1.);
				std::ifstream data_from_txt((path + temp_txt_name + ".txt").c_str());
				getline(data_from_txt, line);
			}
			std::stringstream ss(line);
			ss >> temp_input;
			data_from_txt.close();
			return temp_input;
		} 
		return temp_input;
	}	
    	
		
    template <int dim>
	double
	MeltConveyTwo<dim>::
	read_from_txt(const std::string path,
	              const unsigned int mpi_communicator)const
	{
		double temp_input = 0;
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			std::ifstream data_from_txt(path.c_str());
			std::string line;
			getline(data_from_txt, line);
			while(line == "")
			{
				data_from_txt.close();
				sleep(1.);
				std::ifstream data_from_txt(path.c_str());
				getline(data_from_txt, line);
			}
			std::stringstream ss(line);
			ss >> temp_input;
			data_from_txt.close();
			return temp_input;
		} 
		return temp_input;
	}	
	
	template <int dim>
	void
	MeltConveyTwo<dim>::
    convey_flag_update(const std::string path,
	                   const unsigned int crood_x,
					   const unsigned int crood_y_in_flag_container,
					   const unsigned int mpi_communicator,
					   const unsigned int flag)const
	{
		if(this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == mpi_communicator)
		{
			const int flag_place = 2 * crood_y_in_flag_container;
			const double write_information = flag;
			std::stringstream get_name;
			get_name << crood_x;
			std::string temp_txt_name = get_name.str();
			std::fstream flag_container((path + temp_txt_name + ".txt").c_str(), std::ios::in | std::ios::out);
			while(!flag_container)
			{
				flag_container.close();
				sleep(1.);
				std::fstream flag_container((path + temp_txt_name + ".txt").c_str(), std::ios::in | std::ios::out);
			}
			std::stringstream ss;
			ss << write_information;
			flag_container.seekp(flag_place, std::ios::beg);
			flag_container << ss.str();
			flag_container.close(); 
		}
	}
    	
    template <int dim>
    std::pair<std::string,std::string>
    MeltConveyTwo<dim>::execute (TableHandler &statistics)
    {
      if (this->n_compositional_fields() == 0)
        return std::pair<std::string,std::string>();

      // create a quadrature formula based on the compositional element alone.
      // be defensive about determining that a compositional field actually exists
      AssertThrow (this->introspection().base_elements.compositional_fields
                   != numbers::invalid_unsigned_int,
                   ExcMessage("This postprocessor cannot be used without compositional fields."));
	  
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

	  MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
	  

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
			//获取mpi通道编号
			unsigned int mpi_comm_id = 0;
			unsigned int out_comm_id = 0;
		    while (out_comm_id < 5000)
		    {
			    if (Utilities::MPI::this_mpi_process (this->get_mpi_communicator()) == out_comm_id)
			    {
				    mpi_comm_id = out_comm_id;
				    break;
			    } 
			    out_comm_id ++;
		    } 
			//填充输入
            fe_values.reinit (cell);
			in.reinit(fe_values, cell, this->introspection(), this->get_solution());
			const unsigned int mid_point = n_q_points / 2;
			const double posx = in.position[mid_point][0];
			const double posy = in.position[mid_point][1];
			const double depth = this->get_geometry_model().depth(in.position[mid_point]);
	        const unsigned int crood_x = std::floor( posx / ((max_x / repetitions_x) / std::pow(2, melt_fraction_refinement)));
			const unsigned int crood_y = std::floor( posy / ((max_y / repetitions_y) / std::pow(2, melt_fraction_refinement)));
			const int crood_y_in_flag_container = crood_y - min_container_amount_y;
			
			//2.熔体运输
			//将容器内的熔体运输到指定位置
			//熔体在提取时已经限定了窗口，因此在计算岩浆侵位时不需要考虑窗口的问题。
			//20250701
		    //调整成分场调用逻辑，可在prm文件中删除不需要的成分场
			std::vector<double> sd1_values(n_q_points, 0.);
			std::vector<double> sd2_values(n_q_points, 0.);
			std::vector<double> upper_values(n_q_points, 0.);
			std::vector<double> lower_values(n_q_points, 0.);
			std::vector<double> mantle_upper_values(n_q_points, 0.);
			std::vector<double> mantle_middle_values(n_q_points, 0.);
			std::vector<double> mantle_lower_values(n_q_points, 0.);
			std::vector<double> target_mantle_melting_values(n_q_points, 0.);
			std::vector<double> target_crustal_melting_values(n_q_points, 0.);
			std::vector<double> melt_fraction_values(n_q_points, 0.);
			
			//20260127 添加膏盐模式相关的相，膏岩模式与其他模式共用sd1
			std::vector<double> salt_sd_init_values(n_q_points, 0.);
			std::vector<double> salt_values(n_q_points, 0.);
			std::vector<double> salt_crust_values(n_q_points, 0.);
			std::vector<double> salt_mantle_l_values(n_q_points, 0.);
			std::vector<double> salt_mantle_d_values(n_q_points, 0.);
			
			double sd1_average = 0;
			double sd2_average = 0;
			double upper_average = 0;
			double lower_average = 0;
			//20240701 在传输判定时，所有岩石圈地幔绑定到一起
			double mantle_average = 0;
			
			//20260127 添加膏盐模式相关的相，膏岩模式与其他模式共用sd1
			//在传输判定时，所有岩石圈地幔绑定到一起，所有沉积层和膏盐绑定到一起
			double salt_sd_average = 0.;
			double salt_crust_average = 0.;
			double salt_mantle_average = 0.;
			
			double target_mantle_melting_average = 0;
			double target_crustal_melting_average = 0;
			double fraction_max = 0;
			double area_in_km = 0;
			//20260312 添加温度限制，按等温线拆分到melt convey 2-4里
			//double temp_T = in.temperature[mid_point];
			//20260203 先用膏盐模式判定，减少if次数
			//sediment_1为两模式共用
			if(this->introspection().compositional_name_exists("sediment_1"))
			{
				const unsigned int sd1_id = this->introspection().compositional_index_for_name("sediment_1");
			    fe_values[this->introspection().extractors.compositional_fields[sd1_id]].get_function_values (this->get_solution(), sd1_values);
			}
			//20260318 修改判定逻辑 默认洋壳和熔融程度存在
			const unsigned int melt_fraction_id = this->introspection().compositional_index_for_name("melt_fraction");
		    fe_values[this->introspection().extractors.compositional_fields[melt_fraction_id]].get_function_values (this->get_solution(), melt_fraction_values);
			const unsigned int target_mantle_melting_id = this->introspection().compositional_index_for_name("new_crust");
			fe_values[this->introspection().extractors.compositional_fields[target_mantle_melting_id]].get_function_values (this->get_solution(), target_mantle_melting_values);
			/* if (include_mantle_melting && this->introspection().compositional_name_exists("new_crust"))
			{
				const unsigned int target_mantle_melting_id = this->introspection().compositional_index_for_name("new_crust");
				fe_values[this->introspection().extractors.compositional_fields[target_mantle_melting_id]].get_function_values (this->get_solution(), target_mantle_melting_values);
			}
			if (include_crustal_melting && this->introspection().compositional_name_exists("core_complex"))
			{
				const unsigned int target_crustal_melting_id = this->introspection().compositional_index_for_name("core_complex");
				fe_values[this->introspection().extractors.compositional_fields[target_crustal_melting_id]].get_function_values (this->get_solution(), target_crustal_melting_values);
			}
			if(this->introspection().compositional_name_exists("melt_fraction"))
			{
				unsigned int melt_fraction_id = this->introspection().compositional_index_for_name("melt_fraction");
				fe_values[this->introspection().extractors.compositional_fields[melt_fraction_id]].get_function_values (this->get_solution(), melt_fraction_values);
			}		*/
			
			if(!salt_mode) 
			{
				//20260320 默认上地壳和上地幔存在
				const unsigned int upper_id = this->introspection().compositional_index_for_name("upper");
				fe_values[this->introspection().extractors.compositional_fields[upper_id]].get_function_values (this->get_solution(), upper_values);
				const unsigned int mantle_upper_id = this->introspection().compositional_index_for_name("mantle_upper");
				fe_values[this->introspection().extractors.compositional_fields[mantle_upper_id]].get_function_values (this->get_solution(), mantle_upper_values); 
					
				//20260320 暂时屏蔽
				/* if(this->introspection().compositional_name_exists("sediment_2"))
				{
					const unsigned int sd2_id = this->introspection().compositional_index_for_name("sediment_2");
					fe_values[this->introspection().extractors.compositional_fields[sd2_id]].get_function_values (this->get_solution(), sd2_values); 
				} */
				/* if(this->introspection().compositional_name_exists("upper"))
				{
					
				} */
				if(this->introspection().compositional_name_exists("lower"))
				{
					const unsigned int lower_id = this->introspection().compositional_index_for_name("lower");
					fe_values[this->introspection().extractors.compositional_fields[lower_id]].get_function_values (this->get_solution(), lower_values);
				}
				/* if(this->introspection().compositional_name_exists("mantle_upper"))
				{
					
				} */
				if(this->introspection().compositional_name_exists("mantle_middle"))
				{
					const unsigned int mantle_middle_id = this->introspection().compositional_index_for_name("mantle_middle");
					fe_values[this->introspection().extractors.compositional_fields[mantle_middle_id]].get_function_values (this->get_solution(), mantle_middle_values); 
				}
				if(this->introspection().compositional_name_exists("mantle_lower"))
				{
					const unsigned int mantle_lower_id = this->introspection().compositional_index_for_name("mantle_lower");
					fe_values[this->introspection().extractors.compositional_fields[mantle_lower_id]].get_function_values (this->get_solution(), mantle_lower_values); 
				}
				//20260320 暂时屏蔽
				/* if (include_crustal_melting && this->introspection().compositional_name_exists("core_complex"))
			    {
				const unsigned int target_crustal_melting_id = this->introspection().compositional_index_for_name("core_complex");
				fe_values[this->introspection().extractors.compositional_fields[target_crustal_melting_id]].get_function_values (this->get_solution(), target_crustal_melting_values);
			    } */
			}
			else
			{
				//20260318 默认crust和mantle_L存在
				const unsigned int crust_id = this->introspection().compositional_index_for_name("crust");
				fe_values[this->introspection().extractors.compositional_fields[crust_id]].get_function_values (this->get_solution(), salt_crust_values); 
				const unsigned int mantle_l_id = this->introspection().compositional_index_for_name("mantle_L");
				fe_values[this->introspection().extractors.compositional_fields[mantle_l_id]].get_function_values (this->get_solution(), salt_mantle_l_values); 
				//20260127 添加膏盐模式相关的相，膏岩模式与其他模式共用sd1
				if(this->introspection().compositional_name_exists("sediemnt_init"))
				{
					const unsigned int sd_init_id = this->introspection().compositional_index_for_name("sediemnt_init");
					fe_values[this->introspection().extractors.compositional_fields[sd_init_id]].get_function_values (this->get_solution(), salt_sd_init_values); 
				}
				if(this->introspection().compositional_name_exists("salt"))
				{
					const unsigned int salt_id = this->introspection().compositional_index_for_name("salt");
					fe_values[this->introspection().extractors.compositional_fields[salt_id]].get_function_values (this->get_solution(), salt_values); 
				}
				/* if(this->introspection().compositional_name_exists("crust"))
				{
					
				}
				if(this->introspection().compositional_name_exists("mantle_L"))
				{
					
				} */
				if(this->introspection().compositional_name_exists("mantle_D"))
				{
					const unsigned int mantle_d_id = this->introspection().compositional_index_for_name("mantle_D");
					fe_values[this->introspection().extractors.compositional_fields[mantle_d_id]].get_function_values (this->get_solution(), salt_mantle_d_values); 
				}
			}

            			
			 for (unsigned int q = 0; q < n_q_points; q++)
			{
				sd1_average += std::max(0., sd1_values[q] * fe_values.JxW(q) / 1e6);
				sd2_average += std::max(0., sd2_values[q] * fe_values.JxW(q) / 1e6);
				upper_average += std::max(0., upper_values[q] * fe_values.JxW(q) / 1e6);
				lower_average += std::max(0., lower_values[q] * fe_values.JxW(q) / 1e6);
				mantle_average += std::max(0., (mantle_upper_values[q] + mantle_middle_values[q] + mantle_lower_values[q]) * fe_values.JxW(q) / 1e6);
				
				//20260127 添加膏盐模式相关的相，膏岩模式与其他模式共用sd1
				//20260310 膏盐与地壳相共同计算
				//salt_sd_average += std::max(0., sd1_values[q] * fe_values.JxW(q) / 1e6) + std::max(0., salt_sd_init_values[q] * fe_values.JxW(q) / 1e6) + std::max(0., salt_values[q] * fe_values.JxW(q) / 1e6);
				//salt_crust_average += std::max(0., salt_crust_values[q] * fe_values.JxW(q) / 1e6);
				salt_sd_average += std::max(0., sd1_values[q] * fe_values.JxW(q) / 1e6) + std::max(0., salt_sd_init_values[q] * fe_values.JxW(q) / 1e6);
				salt_crust_average += std::max(0., salt_values[q] * fe_values.JxW(q) / 1e6) + std::max(0., salt_crust_values[q] * fe_values.JxW(q) / 1e6);
				salt_mantle_average += std::max(0., salt_mantle_l_values[q] * fe_values.JxW(q) / 1e6) + std::max(0., salt_mantle_d_values[q] * fe_values.JxW(q) / 1e6);
				
				target_mantle_melting_average += std::max(0., target_mantle_melting_values[q] * fe_values.JxW(q) / 1e6);
				target_crustal_melting_average += std::max(0., target_crustal_melting_values[q] * fe_values.JxW(q) / 1e6);
				area_in_km += fe_values.JxW(q) / 1e6;
				fraction_max = std::max(fraction_max, std::abs(melt_fraction_values[q]));
			} 
								
			//地幔熔融传输条件2：新生地壳的边缘区域（主要成分非上地壳）
			//暂时规定两种熔融产生的物质不能混在一起
			const bool mantle_convey_start_2 = upper_average / area_in_km < 0.3
											 && target_mantle_melting_average / area_in_km > 0.02
											 //20250222
											 //&& target_mantle_melting_average / area_in_km < 0.9
											 && target_mantle_melting_average / area_in_km < 0.5
											 
											 && mantle_average / area_in_km > 0.02
											 //20250220
											 //20250226
											 && sd1_average / area_in_km < 0.02
											 && sd2_average / area_in_km < 0.02
											 && target_crustal_melting_average / area_in_km < 0.02;
											 
			//膏岩模式下地幔熔融传输条件2：新生地壳的边缘区域（主要成分非地壳，即在莫霍面附近及其下方产生）
			const bool salt_convey_start_2 = 
			                                  //20260310
			                                  //salt_crust_average / area_in_km < 0.3
											  //20260317
											  salt_crust_average / area_in_km < 0.8
			                               && target_mantle_melting_average / area_in_km > 0.02
										   //20260314
										   //&& target_mantle_melting_average / area_in_km < 0.5
										   //&& target_mantle_melting_average / area_in_km < 0.3
										   //20260317
										   && target_mantle_melting_average / area_in_km < 0.8
										   && salt_mantle_average / area_in_km > 0.02
										   && salt_sd_average / area_in_km < 0.02
										   && target_crustal_melting_average / area_in_km < 0.02;
											 
			//地壳熔融传输条件2：核杂岩边缘区域（主要成分非下地壳）
			//暂时规定两种熔融产生的物质不能混在一起
			const bool crustal_convey_start_2 = lower_average / area_in_km < 0.3
			                                 && target_crustal_melting_average / area_in_km > 0.02
											 && target_crustal_melting_average / area_in_km < 0.5
											 && mantle_average / area_in_km < 0.02
											 && sd1_average / area_in_km < 0.02
											 && sd2_average / area_in_km < 0.02
											 && target_mantle_melting_average / area_in_km < 0.02;
											 
											
			//考虑到aspect规定通道0为master读写通道，为防止输出结果时出现读写冲突，当mpi通道为0时不占用通道向容器内读写
			//即模型的最左侧（根据精度情况，约为最左侧10~15km）不启用熔体运输功能
			if (include_mantle_melting 
			     && ((!this->introspection().compositional_name_exists("salt") && mantle_convey_start_2) 
				 || (this->introspection().compositional_name_exists("salt") && salt_convey_start_2)) 
				 && fraction_max < 1e-2
			     && this->get_timestep_number() > 0
				 && depth <= depth_max
				 && mpi_comm_id > 0
				 && crood_y_in_flag_container >= 0)
			{
				
			    //20250428
				//含碳熔融传输
				//仅在垂直提取模式中使用
				const double temp_carbon_melt_fraction = vertical_trans && include_carbonate_melting ? read_from_container(carbon_container_path, crood_x, mpi_comm_id) : 0.;
				//规定含碳熔融提取量达到网格总面积2%时即可开始提取
				const double carbon_melt_fraction = 0.02;
			    const double temp_mantle_melt_fraction = read_from_container(container_mantle_melting_path, crood_x, mpi_comm_id);
			    if (temp_mantle_melt_fraction >= container_limit + area_in_km - target_mantle_melting_average || temp_carbon_melt_fraction >= container_limit + (area_in_km - target_mantle_melting_average) * carbon_melt_fraction)
				{
					//写入搬运标记
			        convey_flag_update(convey_flag_mantle_melting_path, crood_x, crood_y_in_flag_container, mpi_comm_id, 1);
				    //更新熔体值
					const double out_mantle_melt_fraction = temp_mantle_melt_fraction - (area_in_km - target_mantle_melting_average);
				    write_to_container(container_mantle_melting_path, crood_x, mpi_comm_id, out_mantle_melt_fraction);
				}
			}
			
			if (include_crustal_melting 
			     && crustal_convey_start_2 
				 && fraction_max < 1e-2
			     && this->get_timestep_number() > 0
				 && depth <= depth_max
				 && mpi_comm_id > 0
				 && crood_y_in_flag_container >= 0)
		    {
				const double temp_crustal_melt_fraction = read_from_container(container_crustal_melting_path, crood_x, mpi_comm_id);
			    if (temp_crustal_melt_fraction >= container_limit + area_in_km - target_crustal_melting_average)
				{
					//写入搬运标记
			        convey_flag_update(convey_flag_crustal_melting_path, crood_x, crood_y_in_flag_container, mpi_comm_id, 1);
				    //更新熔体值
					const double out_crustal_melt_fraction = temp_crustal_melt_fraction - (area_in_km - target_crustal_melting_average);
				    write_to_container(container_crustal_melting_path, crood_x, mpi_comm_id, out_crustal_melt_fraction);
			    }
			}
          }
	  //预留位置，如有需要可视化的参数再进行添加
      statistics.add_value ("Melt Convey 2", 0);
	  {
        const char *columns[] = {"Melt Convey 2"};
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }
      std::ostringstream output;
      output.precision(4);
	  output << "step 2";
      return std::pair<std::string, std::string> ("Processing of melt convey:",
                                                  output.str());
    }
	
	template <int dim>
    void
    MeltConveyTwo<dim>::parse_parameters (ParameterHandler &prm)
	{
		prm.enter_subsection("Postprocess");
		{
			prm.enter_subsection("Melt extraction");
			{
				depth_max = prm.get_double ("Max depth of melt convey function");
				//20240614
				//melt_fraction_refinement = prm.get_integer ("Mesh refinement of melt fraction");
				//20240712
				container_limit = prm.get_double ("Container limit of new crust generation");
				//20250428
				vertical_trans = prm.get_bool("Vertical transportation mode");
			}
			prm.leave_subsection ();
		}
		prm.leave_subsection ();
		
		//20240614
		prm.enter_subsection("Mesh refinement");
		{
			melt_fraction_refinement = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
		}
		prm.leave_subsection ();
		
		prm.enter_subsection("Geometry model");
		{
			prm.enter_subsection("Box");
			{
				max_x = prm.get_double ("X extent");
				repetitions_x = prm.get_integer ("X repetitions");
				max_y = prm.get_double ("Y extent");
				repetitions_y = prm.get_integer ("Y repetitions");
				
				min_container_amount_y = std::ceil((max_y - depth_max - 5e3) / ((max_y / repetitions_y) / std::pow(2, melt_fraction_refinement)));
				convey_flag_mantle_melting_path = this->get_output_directory() + "convey_flag_mantle_melting/";
				container_mantle_melting_path = this->get_output_directory() + "container_mantle_melting/";
				convey_flag_crustal_melting_path = this->get_output_directory() + "convey_flag_crustal_melting/";
				container_crustal_melting_path = this->get_output_directory() + "container_crustal_melting/";
				melt_info_path = this->get_output_directory() + "melt_info/";
				//20250428
				carbon_container_path = this->get_output_directory() + "container_carbonate_melting/";
			}
			prm.leave_subsection ();
		}
		prm.leave_subsection ();
		
		prm.enter_subsection("Material model");
		{
			prm.enter_subsection ("Visco Plastic Grain Size");
			{
				include_mantle_melting = prm.get_bool ("Include mantle melting");
				include_crustal_melting = prm.get_bool("Include crustal melting");
				//20250428
				include_carbonate_melting = prm.get_bool("Include carbonate melting");
			}
			prm.leave_subsection();
		}
		prm.leave_subsection();
		
		//20260203
		prm.enter_subsection ("Initial composition model");
		{
			prm.enter_subsection("Lithosphere with rift");
			{
				//20260203
			    salt_mode = prm.get_bool ("Salt mode");
			}
			prm.leave_subsection();
		}
		prm.leave_subsection();
	}
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MeltConveyTwo,
                                  "melt convey 2",
                                  "A postprocessor that computes some statistics about "
                                  "the compositional fields, if present in this simulation. "
                                  "In particular, it computes maximal and minimal values of "
                                  "each field, as well as the total mass contained in this "
                                  "field as defined by the integral "
                                  "$m_i(t) = \\int_\\Omega c_i(\\mathbf x,t) \\; \\text{d}x$.")
  }
}
