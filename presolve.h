// Copyright (C) 2024 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#ifndef PRESOLVE_H
#define PRESOLVE_H

#include <string>

#include "Eigen/Dense"

namespace simplex
{
	class LinearProgram;

	class Presolve
	{
	public:
		Presolve(LinearProgram & lp, bool reductions_enabled = true)
			: m_lp(lp), m_reductions_enabled(reductions_enabled) {}

		~Presolve() = default;

		void run();

	private:
		LinearProgram & m_lp;
		bool m_reductions_enabled;

		void eliminate_lbs();
		void apply_reductions() const;
		void fix_variable(Eigen::Index var_id, double values);

		void set_shift(const std::string & var_name, double value);
		double get_shift(const std::string & var_name) const;
	};
}

#endif
