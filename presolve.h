// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#ifndef _PRESOLVE_H_
#define _PRESOLVE_H_

#include <string>

namespace simplex
{
	class LinearProgram;

	class Presolve
	{
	public:
		Presolve(LinearProgram & lp, bool reductions_enabled = true)
			: m_lp(lp), b_reductions_enabled(reductions_enabled) {}

		~Presolve() = default;

		void run();

	private:
		LinearProgram & m_lp;
		bool b_reductions_enabled;

		void eliminate_lbs();
		void apply_reductions();
		void fix_variable(const std::string & var_name, double values);

		void set_shift(const std::string & var_name, double value);
		double get_shift(const std::string & var_name) const;
	};
}

#endif
