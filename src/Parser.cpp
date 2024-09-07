//// Copyright (C) 2024 Maciej Drwal
////
//// Permission is granted to copy and distribute verbatim copies and modified
//// versions of this file, provided that the copyright notice and this permission
//// notice are preserved on all copies and modified versions of this file.
////

#include "Parser.h"
#include "Utils.h"
#include "Logger.h"

#include <algorithm>
#include <locale>
#include <string>

namespace
{
    std::string to_lower(std::string str)
    {
        std::transform(str.begin(), str.end(), str.begin(), [](char c) { return std::tolower(c); });
        return str;
    }

    std::string strip_chars(std::string line, const std::function<bool(char)> & pred)
    {
        line.erase(std::remove_if(line.begin(), line.end(), pred), line.end());
        return line;
    }

    std::string strip_spaces(std::string & line)
    {
        return strip_chars(line, [](char c) { return isspace(c); });
    }

    std::string parse_line(std::string::const_iterator & iter, std::string::const_iterator end)
    {
        std::string line;
        while (line.empty() && iter != end)
        {
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
            line = strip_spaces(line);
        }
        return line;
    }

    std::pair<std::string, std::string> parse_line_until(std::string::const_iterator & iter, std::string::const_iterator end, const std::string & tag)
    {
        std::string text;
        while (iter != end)
        {
            const auto line = parse_line(iter, end);
            const auto line_tlsc = to_lower(strip_chars(line, [](char c) { return c == ':'; })); 
            if (line_tlsc == tag || line_tlsc == "end")
            {
                return { text, line_tlsc };
            }
            text += line;
        }
        return { text, "" };
    }

    std::pair<std::string, std::string> parse_constraint_line_until(std::string::const_iterator & iter, std::string::const_iterator end, const std::string & tag)
    {
        std::string text;
        while (iter != end)
        {
            const auto line = parse_line(iter, end);
            const auto line_tlsc = to_lower(strip_chars(line, [](char c) { return c == ':'; }));
            if (line_tlsc == tag || line_tlsc == "end")
            {
                return { text, line_tlsc };
            }
            text += line;
            if (line.find_first_of("<=>") != std::string::npos)
            {
				break;
			}
        }
        return { text, "" };
    }

    void set_sense(std::string line, simplex::LinearProgram & lp)
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

    std::string as_keyword(const std::string & line)
    {
        std::string tmp_line(line);
        return to_lower(strip_chars(tmp_line, [](char c) { return c == ':'; }));
    }

    void expect_keyword(std::string line, std::string_view keyword)
    {
        line = to_lower(strip_chars(line, [](char c) { return c == ':'; }));
        if (line != keyword)
        {
            throw "Expected keyword: " + std::string(keyword);
        }
    }

    std::pair<std::string, std::string> extract_label_and_expr(const std::string & line)
    {
        const auto pos = line.find(':');
        if (pos == std::string::npos)
        {
            return { "", line };
        }
        return { line.substr(0, pos), line.substr(pos + 1) };
    }

    std::vector<std::string> tokenize_expr(const std::string & expr)
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

    void set_objective(const std::string & expr, simplex::LinearProgram & lp)
    {
        const auto tokens = tokenize_expr(expr);
        for (const auto & token : tokens)
        {
            const auto pos = std::find_if(token.begin(), token.end(), [](char c) { return std::isalpha(c); });
            if (pos == token.end())
            {
                throw "Invalid token in objective function: " + token;
            }
            const auto empty = pos == token.begin();
            const auto single = pos - token.begin() == 1 && (token[0] == '+' || token[0] == '-');
            const auto coeff = (!empty && !single) ? std::stod(token.substr(0, pos - token.begin()))
                               : (token[0] == '-') ? -1.0
                                                   : 1.0;
            const auto var_name = token.substr(pos - token.begin());

            if (lp.has_variable(var_name))
            {
                throw "Variable already defined in objective function: " + var_name;
            }

            lp.add_variable(var_name, coeff);
        }
    }

    void add_constraint(const std::string & label, const std::string & expr, simplex::LinearProgram & lp)
    {
        auto pos = expr.find_first_of("<=>");
        const auto sign = expr.at(pos);
        if (pos == std::string::npos)
        {
            throw "Invalid constraint expression: " + expr.back();
        }
        const auto lhs = expr.substr(0, pos);
        while (!isdigit(expr.at(pos + 1)) && expr.at(pos + 1) != '-')
        {
            ++pos;
        }
        const auto rhs = expr.substr(pos + 1);

        auto tokens = tokenize_expr(lhs);
        if (tokens.empty())
        {
            throw "Empty constraint expression";
        }

        simplex::Constraint constraint(sign, stod(rhs));

        //tokens.back() = tokens.back().substr(0, tokens.back().find(sign));  // note: this invalidates pos

        for (const auto & token : tokens)
        {
            const auto pos1 = std::find_if(token.begin(), token.end(), [](char c) { return std::isalpha(c); });
            const auto empty = pos1 == token.begin();
            const auto single = pos1 - token.begin() == 1 && (token[0] == '+' || token[0] == '-');
            const auto coeff =
                (!empty && !single) ? std::stod(token.substr(0, pos1 - token.begin())) : (token[0] == '-' ? -1.0 : 1.0);
            const auto var_name = token.substr(pos1 - token.begin());
            if (!lp.has_variable(var_name))
            {
                lp.add_variable(var_name);
            }
            constraint.add_term(var_name, coeff);
        }

        lp.add_constraint(label, std::move(constraint));
    }

    void parse_bounds(const std::string & line, simplex::LinearProgram & lp)
    {
        std::vector<std::string> tokens;
        const auto pos1 = line.find_first_of("<=");
        const auto pos2 = line.find_last_of("<=");
        const auto single_bound = pos2 - pos1 <= 1;
        const auto term_1st = line.substr(0, pos1);
        const auto term_3rd = line.substr(pos2 + 1);

        double low_value = 0.0, high_value = std::numeric_limits<double>::infinity();
        std::string var_name;
        if (single_bound)
        {
            // check if the first term is a number
            if (std::find_if(term_1st.begin(), term_1st.end(), ::isalpha) == term_1st.end())
            {
                low_value = std::stod(term_1st);
                var_name = term_3rd;
            }
            else
            {
                high_value = std::stod(term_3rd);
                var_name = term_1st;
            }
        }
        else
        {
            const auto term_2nd = line.substr(pos1 + 2, pos2 - pos1 - 3);
            low_value = std::stod(term_1st);
            high_value = std::stod(term_3rd);
            var_name = term_2nd;
        }

        const auto var_id = lp.get_variable_id(var_name);

        if (!utils::is_float_zero(low_value))
        {
            lp.set_lower_bound(var_id, low_value);
        }

        if (high_value < std::numeric_limits<double>::infinity())
        {
            lp.set_upper_bound(var_id, high_value);
        }
    }

    bool parse_impl(const std::string & input, simplex::LinearProgram & lp)
    {
        enum ParserState
        {
            SetSense,
            Objective,
            Constraints,
            Bounds
        } state = SetSense;

        auto it = input.begin();
        const auto end = input.end();

        while (it != end)
        {
            std::pair<std::string, std::string> label_and_expr;
            std::string line, tag;

            switch (state)
            {
                case SetSense:
                    line = parse_line(it, end);
                    set_sense(line, lp);
                    state = Objective;
                    break;

                case Objective:
                    std::tie(line, tag) = parse_line_until(it, end, "subjectto");
                    label_and_expr = extract_label_and_expr(line);
                    set_objective(label_and_expr.second, lp);
                    state = Constraints;
                    break;
                case Constraints:
                    std::tie(line, tag) = parse_constraint_line_until(it, end, "bounds");
                    if (as_keyword(tag) == "bounds")
                    {
                        state = Bounds;
                        break;
                    }
                    if (as_keyword(tag) == "end")
                    {
                        return true;
                    }
                    label_and_expr = extract_label_and_expr(line);
                    add_constraint(label_and_expr.first, label_and_expr.second, lp);
                    break;
                case Bounds:
                    line = parse_line(it, end);
                    if (as_keyword(line) == "end")
                    {
                        return true;
                    }
                    parse_bounds(line, lp);
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
}  // namespace

namespace simplex
{
    bool run_parser_lp(std::string_view input, LinearProgram & lp)
    {
        MEASURE_TIME_START(Parser);
        const std::string stripped_input = remove_comments(input);
        const auto result = parse_impl(stripped_input, lp);
        LOG(info) << PRINT_ELAPSED_TIME(Parser);
        return result;
    }
}  // namespace simplex
