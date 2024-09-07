
#include "InitialBasis.h"
#include "Simplex.h"
#include "Logger.h"

namespace simplex
{
    Basis InitialBasis::get_basis() const
    {
        Basis init_basis;

        m_lp.add_slack_variables_for_inequality_constraints();

        // Construct artificial variables for Phase I.
        const auto art_var_id = m_lp.add_artificial_variables_for_first_phase(init_basis);

        if (art_var_id > 0)
        {
            // note: number of variables N has changed after adding artificial variables

            // Solve Phase I LP.
            LOG(debug) << "*** PHASE I ***";
            m_lp.initialize_tableau();

            // Objective function: sum of artificial variables.
            m_lp.vector_c.setZero();
            for (int i = 0; i < art_var_id; i++)
            {
                const auto var_name = LinearProgram::get_artificial_variable(i);
                m_lp.vector_c[m_lp.variable_name_to_id[var_name]] = 1.0;
            }

            // TODO: filling the non-basic columns (the passed initial basis may contain only the basic indices)
            for (auto i = 0u; i < m_lp.get_num_vars(); i++)
            {
                if (!init_basis.contains(i))
                {
                    init_basis.insert_to_non_basic(i);
                }
            }

            // m_lp.print_tableau();
            m_simplex.simplex(init_basis);

            // Remove artificial variables.
            for (auto & [label, constraint] : m_lp.constraints)
            {
                const auto it = m_lp.m_artificials.find(label);
                if (it == m_lp.m_artificials.end())
                {
                    continue;
                }

                const auto var_id = it->second;
                const auto & var_name = m_lp.variable_id_to_name[var_id];

                // Check if some artificial variable remained in the basis.
                if (init_basis.contains(var_id))
                {
                    // TODO: handle this case

                    // init_basis.erase(var_name);
                    // for (map<string, int>::iterator jt = variable_name_to_id.begin(); jt !=
                    // variable_name_to_id.end(); ++jt) {
                    //     if (init_basis.count(jt->second) == 0 && (strncmp(jt->first.c_str(), ARTIFICIAL, 12) != 0)) {
                    //         init_basis.insert(jt->second);
                    //         printf("replaced variable %d (%s) by %d (%s) in the basis\n",  var_id, var_name.c_str(),
                    //         jt->second, jt->first.c_str()); break;
                    //     }
                    // }

                    std::stringstream msg;
                    msg << "Artificial variable " << var_name << " in the basis.";
                    throw msg.str().c_str();
                }

                constraint.remove_term(var_name);
                m_lp.remove_variable(var_id);
                m_lp.variable_name_to_id.erase(var_name);
                init_basis.remove(var_id);
                LOG(debug) << "removed artificial variable: " << var_name;
            }

            // note: number of variables N has changed after removing artificial variables

            m_lp.vector_b.setZero();
            m_lp.vector_c.setZero();
            m_lp.matrix_A.setZero();
        }
        else
        {
            LOG(debug) << "Using slacks for initial basis, skipping Phase I.";

            // TODO: filling the non-basic columns (the passed initial basis may contain only the basic indices)
            for (auto i = 0u; i < m_lp.get_num_vars(); i++)
            {
                if (!init_basis.contains(i))
                {
                    init_basis.insert_to_non_basic(i);
                }
            }
        }

        return init_basis;
    }
}