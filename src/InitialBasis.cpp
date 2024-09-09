
#include "InitialBasis.h"
#include "Logger.h"

namespace
{
    const std::string ARTIFICIAL = "_ARTIFICIAL_";

    void fill_non_basic_variables(simplex::Basis & basis, const simplex::LinearProgram & lp)
    {
        for (auto i = 0u; i < lp.get_num_vars(); i++)
        {
            if (!basis.contains(i))
            {
                basis.insert_to_non_basic(i);
            }
        }
    }
}

namespace simplex
{
    size_t InitialBasis::add_artificial_variables_for_initial_basis(Basis & basis)
    {
        int art_var_id = 0;
        for (auto & [label, constraint] : m_lp.constraints)
        {
            // Note: at this point there should be no inequality constraints
            if (constraint.type == '<' || constraint.type == '>')
            {
                throw "All constraints must be equality at this point.";
            }

            const auto rhs = constraint.rhs;

            const auto it = m_lp.m_slacks.find(label);
            if (it != m_lp.m_slacks.end() && rhs >= 0.0)
            {
                // Slack is available for this constraint, but can be used only if the right-hand side is nonnegative.
                basis.insert(it->second);
            }
            else
            {
                // Artificial variable is introduced with negative sign if the right-hand side is negative.
                const std::string var_name(ARTIFICIAL + std::to_string(art_var_id++));
                constraint.add_term(var_name, rhs >= 0.0 ? 1.0 : -1.0);
                const auto var_id = m_lp.add_variable(var_name);
                basis.insert(var_id);
                m_artificials[label] = var_id;
                LOG(debug) << "added artificial variable: " << var_name;
            }
        }
        LOG(debug) << "added " << art_var_id << " artificial variables";
        return art_var_id;
    }

    Basis InitialBasis::get_basis()
    {
        Basis init_basis;

        m_lp.add_slack_variables_for_inequality_constraints();

        // Construct artificial variables for Phase I.
        const auto art_var_id = add_artificial_variables_for_initial_basis(init_basis);

        if (art_var_id > 0)
        {
            // Solve Phase I LP.
            LOG(debug) << "*** PHASE I ***";
            m_lp.initialize_tableau();

            // Objective function: sum of artificial variables.
            m_lp.vector_c.setZero();
            for (size_t i = 0u; i < art_var_id; i++)
            {
                const std::string var_name(ARTIFICIAL + std::to_string(i));
                m_lp.vector_c[m_lp.variable_name_to_id[var_name]] = 1.0;
            }

            fill_non_basic_variables(init_basis, m_lp);

            m_simplex.simplex(init_basis);

            // Remove artificial variables.
            for (auto & [label, constraint] : m_lp.constraints)
            {
                const auto it = m_artificials.find(label);
                if (it == m_artificials.end())
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
                    msg << "artificial variable " << var_name << " in the basis.";
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
            LOG(debug) << "Using all slacks for initial basis, skipping Phase I.";

            fill_non_basic_variables(init_basis, m_lp);
        }

        return init_basis;
    }
}