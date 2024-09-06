// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "Basis.h"
#include "Logger.h"

namespace simplex
{
    //Basis::Basis(size_t N)
    //{
    //    m_non_basic_cols.reserve(N);
    //    for (size_t i = 0; i < N; ++i)
    //    {
    //        m_non_basic_cols.push_back(i);
    //    }
    //}

    void Basis::insert(Eigen::Index var_id)
    {
        m_basic_cols.push_back(var_id);
        const auto it = std::find(m_non_basic_cols.begin(), m_non_basic_cols.end(), var_id);
        if (it != m_non_basic_cols.end()) m_non_basic_cols.erase(it);
    }

    void Basis::insert_to_non_basic(Eigen::Index var_id)
    {
        m_non_basic_cols.push_back(var_id);
    }

    void Basis::remove(Eigen::Index var_id)
    {
        const auto it1 = std::find(m_basic_cols.begin(), m_basic_cols.end(), var_id);
        if (it1 != m_basic_cols.end())
        {
            m_basic_cols.erase(it1);
            return;
        }
        const auto it2 = std::find(m_non_basic_cols.begin(), m_non_basic_cols.end(), var_id);
        if (it2 != m_non_basic_cols.end()) m_non_basic_cols.erase(it2);
    }

    void Basis::clear()
    {
        m_basic_cols.clear();
        m_non_basic_cols.clear();
    }

    bool Basis::contains(Eigen::Index var_id) const
    {
        return std::find(m_basic_cols.begin(), m_basic_cols.end(), var_id) != m_basic_cols.end();
    }

    void Basis::update(Eigen::Index entering_index, Eigen::Index leaving_index)
    {
        const auto leaving_col = m_basic_cols[leaving_index];
        const auto entering_col = m_non_basic_cols[entering_index];
        m_basic_cols[leaving_index] = entering_col;
        m_non_basic_cols[entering_index] = leaving_col;

        LOG(debug) << "leaving basis: " << leaving_col << " (" << leaving_index
                       << "), entering basis: " << entering_col << " (" << entering_index << ")";
    }
}