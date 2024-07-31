//// Copyright (C) 2024 Maciej Drwal
////
//// Permission is granted to copy and distribute verbatim copies and modified
//// versions of this file, provided that the copyright notice and this permission
//// notice are preserved on all copies and modified versions of this file.
////

#include <algorithm>
#include <locale>
#include <string>

#include "Parser.h"

namespace
{
    std::string to_lower(std::string str)
    {
        std::transform(str.begin(), str.end(), str.begin(), [](char c)
                       { return std::tolower(c); });
        return str;
    }

    std::string parse_line(std::string::const_iterator &iter, std::string::const_iterator end)
    {
        std::string line;
        while (iter != end)
        {
            if (*iter == '\n')
            {
                ++iter;
                break;
            }
            line += *iter;
            ++iter;
        }
        return line;
    }

    void set_sense(std::string line, simplex::LinearProgram &lp)
    {
        line = to_lower(line);

        if (line.size() < 3)
        {
            throw "Expected keyword: Minimize or Maximize";
        }

        if (line.substr(0, 3) == "min")
        {
            lp.set_sense('m');
        }
        if (line.substr(0, 3) == "max")
        {
            lp.set_sense('M');
        }
    }

    void expect_keyword(std::string line, std::string_view keyword)
    {
        if (to_lower(line) != keyword)
        {
            throw "Expected keyword: " + std::string(keyword);
        }
    }

    std::string strip_spaces(std::string &line)
    {
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
        return line;
    }

    std::pair<std::string, std::string> extract_label_and_expr(const std::string &line)
    {
        const auto pos = line.find(':');
        if (pos == std::string::npos)
        {
            return {"", line};
        }
        return {line.substr(0, pos), line.substr(pos + 1)};
    }

    std::vector<std::string> tokenize_expr(const std::string &expr)
    {
        std::vector<std::string> tokens;
        std::string token;
        for (const char c : expr)
        {
            if ((c == '+' || c == '-') && !token.empty())
            {
                tokens.push_back(token);
                token.clear();
            }
            token += c;
        }
        tokens.push_back(token);
        return tokens;
    }

    void set_objective(const std::string &expr, simplex::LinearProgram &lp)
    {
        const auto tokens = tokenize_expr(expr);
        for (const auto &token : tokens)
        {
            const auto pos = std::find_if(token.begin(), token.end(), [](char c)
                                          { return std::isalpha(c); });
            if (pos == token.end())
            {
                throw "Invalid token in objective function: " + token;
            }
            const auto empty = pos == token.begin();
            const auto single = pos - token.begin() == 1 && (token[0] == '+' || token[0] == '-');
            const auto coeff = (!empty && !single) ? std::stod(token.substr(0, pos - token.begin())) : (token[0] == '-') ? -1.0
                                                                                                                         : 1.0;
            const auto var_name = token.substr(pos - token.begin());

            if (lp.has_variable(var_name))
            {
                throw "Variable already defined in objective function: " + var_name;
            }

            lp.add_variable(var_name, coeff);
        }
    }

    void add_constraint(const std::string &label, const std::string &expr, simplex::LinearProgram &lp)
    {
        auto tokens = tokenize_expr(expr);
        if (tokens.empty())
        {
            throw "Empty constraint expression";
        }

        const auto pos = tokens.back().find_first_of("<=>");
        if (pos == std::string::npos)
        {
            throw "Invalid constraint expression: " + tokens.back();
        }
        const auto rhs = std::stod(tokens.back().substr(pos + 1));
        simplex::Constraint constraint(tokens.back().at(pos), rhs);

        tokens.back() = tokens.back().substr(0, pos); // note: this invalidates pos

        for (const auto &token : tokens)
        {
            const auto pos1 = std::find_if(token.begin(), token.end(), [](char c)
                                           { return std::isalpha(c); });
            const auto empty = pos1 == token.begin();
            const auto single = pos1 - token.begin() == 1 && (token[0] == '+' || token[0] == '-');
            const auto coeff = (!empty && !single) ? std::stod(token.substr(0, pos1 - token.begin())) : (token[0] == '-' ? -1.0 : 1.0);
            const auto var_name = token.substr(pos1 - token.begin());
            constraint.add_term(var_name, coeff);
        }

        lp.constraints.emplace(label, constraint);
    }

    bool parse_impl(const std::string &input, simplex::LinearProgram &lp)
    {
        enum ParserState
        {
            SetSense,
            Objective,
            SubjectTo,
            Constraints,
            Bounds
        } state = SetSense;

        auto it = input.begin();
        const auto end = input.end();

        while (it != end)
        {
            std::pair<std::string, std::string> label_and_expr;
            std::string line = parse_line(it, end);
            strip_spaces(line);
            if (line.empty())
            {
                continue;
            }

            switch (state)
            {
            case SetSense:
                set_sense(line, lp);
                state = Objective;
                break;

            case Objective:
                label_and_expr = extract_label_and_expr(line);
                set_objective(label_and_expr.second, lp);
                state = SubjectTo;
                break;

            case SubjectTo:
                expect_keyword(line, "subjectto");
                state = Constraints;
                break;
            case Constraints:
                label_and_expr = extract_label_and_expr(line);
                strip_spaces(label_and_expr.second);
                if (to_lower(label_and_expr.second) == "bounds")
                {
                    state = Bounds;
                }
                if (to_lower(label_and_expr.second) == "end")
                {
                    return true;
                }
                add_constraint(label_and_expr.first, label_and_expr.second, lp);
                break;
            case Bounds:
                break;
            }
        }
        return false;
    }

    std::string remove_comments(std::string_view input)
    {
        std::string result;
        size_t pos_start = 0, pos_cur = 0;

        while (true)
        {
            pos_cur = input.find_first_of("/\\#", pos_start);
            if (pos_cur == std::string::npos)
            {
                pos_cur = input.length();
                break;
            }
            result += input.substr(pos_start, pos_cur - pos_start);
            pos_start = input.find_first_of('\n', pos_cur) + 1;
        }
        result += input.substr(pos_start, pos_cur - pos_start);

        return result;
    }
}

namespace simplex
{
    bool run_parser_lp(std::string_view input, LinearProgram &lp)
    {
        const std::string stripped_input = remove_comments(input);
        return parse_impl(stripped_input, lp);
    }
}
