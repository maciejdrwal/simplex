// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#ifndef BASIS_H
#define BASIS_H

#include "Eigen/Dense"

#include <vector>
#include <string>

namespace simplex
{
    struct Basis
    {
        [[nodiscard]] const std::vector<Eigen::Index> & get_basic_columns() const { return m_basic_cols; }
        [[nodiscard]] const std::vector<Eigen::Index> & get_non_basic_columns() const { return m_non_basic_cols; }

        void insert(Eigen::Index var_id);
        void insert_to_non_basic(Eigen::Index var_id);
        void remove(Eigen::Index var_id);
        void clear();
        bool contains(Eigen::Index var_id) const;
        void update(Eigen::Index entering_index, Eigen::Index leaving_index);
        std::string show() const;

    private:
        std::vector<Eigen::Index> m_basic_cols;
        std::vector<Eigen::Index> m_non_basic_cols;
    };

}  // namespace simplex

#endif
