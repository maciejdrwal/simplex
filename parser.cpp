// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <boost/spirit/home/x3.hpp>

#include "parser.h"
#include "simplex.h"
#include "utils.h"

namespace simplex
{
    struct ParserState
    {
        ParserState(LinearProgram &_lp) : lp(_lp), var_coeff_cache(1.0), sign_cache(1.0), label_cache("") {}

        LinearProgram &lp;
        double var_coeff_cache;
        double sign_cache;
        std::string label_cache;
    };

    namespace x3 = boost::spirit::x3;
    namespace ascii = boost::spirit::x3::ascii;

    namespace parser
    {
        using x3::double_;
        using x3::eoi;
        using x3::lexeme;
        using x3::lit;
        using x3::no_case;
        using x3::skip;

        using ascii::alnum;
        using ascii::alpha;
        using ascii::char_;
        using ascii::space;

        // Semantic actions definitions
        auto set_minimize = [](auto &ctx)
        { _val(ctx).lp.set_sense('m'); };
        auto set_maximize = [](auto &ctx)
        { _val(ctx).lp.set_sense('M'); };
        auto set_obj_label = [](auto &ctx)
        { _val(ctx).lp.set_objective_label(_attr(ctx)); };
        auto add_obj_fun_coeff = [](auto &ctx)
        { _val(ctx).var_coeff_cache = _attr(ctx); };
        auto add_obj_fun_var = [](auto &ctx)
        {
            // create new variable and store associated obj.fun. coefficient
            auto var_name = _attr(ctx);
            if (_val(ctx).lp.has_variable(var_name))
            {
                throw "Parser: variable already defined in objective function.";
            }
            _val(ctx).lp.add_variable(var_name);
            _val(ctx).lp.objective_name_coeff[var_name] =
                _val(ctx).sign_cache * _val(ctx).var_coeff_cache;
        };
        auto neg_var_coeff = [](auto &ctx)
        { _val(ctx).sign_cache = -1; };
        auto term_parsed = [](auto &ctx)
        {
            _val(ctx).sign_cache = 1;
            _val(ctx).var_coeff_cache = 1;
        };
        auto add_obj_fun_const = [](auto &ctx)
        { _val(ctx).lp.set_obj_value_shift(_attr(ctx)); };

        auto set_constr_label = [](auto &ctx)
        {
            // create labeled constraint
            Constraint constraint;
            _val(ctx).label_cache = _attr(ctx);
            _val(ctx).lp.constraints.emplace(_attr(ctx), constraint);
        };
        auto add_constr_coeff = [](auto &ctx)
        { _val(ctx).var_coeff_cache = _attr(ctx); };
        auto add_constr_var = [](auto &ctx)
        {
            auto &constr_label = _val(ctx).label_cache;
            if (constr_label.empty())
            {
                // no label was defined by user; create new constraint with default label
                Constraint constraint;
                constr_label = std::string("CONSTR") + utils::tostr<size_t>(_val(ctx).lp.constraints.size());
                _val(ctx).label_cache = constr_label;
                _val(ctx).lp.constraints.emplace(constr_label, constraint);
            }

            auto var_name = _attr(ctx);

            // NOTE: it should be possible to accumulate the repeated variable's coefficients
            // if (_val(ctx).lp.constraints[constr_label].has_variable(var_name))
            // {
            //     throw(std::string("Parser: variable already defined in constraint") + constr_label).c_str();
            // }

            if (!_val(ctx).lp.has_variable(var_name))
            {
                // adding new variable via constraint
                _val(ctx).lp.add_variable(var_name);
            }

            auto coeff = _val(ctx).var_coeff_cache;
            auto sign = _val(ctx).sign_cache;
            _val(ctx).lp.constraints[constr_label].add_term(var_name, sign * coeff);
        };

        auto add_constr_type = [](auto &ctx)
        {
            _val(ctx).lp.constraints[_val(ctx).label_cache].type = _attr(ctx);
        };
        auto add_constr_rhs = [](auto &ctx)
        {
            double _rhs = _attr(ctx);
            _val(ctx).lp.constraints[_val(ctx).label_cache].rhs = _rhs;
            auto op = _val(ctx).lp.constraints[_val(ctx).label_cache].type;
            if ((op == '<' && _rhs < 0.0) ||
                (op == '>' && _rhs > 0.0) ||
                (op == '='))
                _val(ctx).lp.set_all_inequalities(false);
            _val(ctx).label_cache = "";
        };

        auto add_lb_constr = [](auto &ctx)
        {
            auto lb = _val(ctx).var_coeff_cache; // now it contains LB
            auto var_name = _attr(ctx);
            _val(ctx).label_cache = var_name; // store var_name, in case we need to add UB
            if (utils::isfloatzero(lb))
            {
                return;
            }
            _val(ctx).lp.var_lbnd[var_name] = lb;
        };
        auto add_ub_constr = [](auto &ctx)
        {
            auto var_name = _val(ctx).label_cache; // now it contains var_name
            auto ub = _attr(ctx);
            _val(ctx).lp.var_ubnd[var_name] = ub;
        };
        auto store_ub_var_name = [](auto &ctx)
        { _val(ctx).label_cache = _attr(ctx); };
        auto add_inf_lb = [](auto &ctx)
        {
            auto lb = std::numeric_limits<double>::min();
            _val(ctx).var_coeff_cache = lb;
        };
        auto add_inf_ub = [](auto &ctx)
        {
            auto var_name = _val(ctx).label_cache; // now it contains var_name
            auto ub = std::numeric_limits<double>::max();
            _val(ctx).lp.var_ubnd[var_name] = ub;
        };

        // Grammar rules definitions
        // TODO: handle properly the constant objective function, e.g., obj: 0

        struct keywords_t : x3::symbols<x3::unused_type>
        {
            keywords_t()
            {
                add("Subject")("Bounds")("End")("Inf")("Infinity");
            }
        } const keywords;

        const auto distinct_keywords = lexeme[no_case[keywords] >> !(alnum | '_')];

        x3::rule<struct identifier_tag, std::string> const identifier("identifier");
        const auto identifier_def = lexeme[+(alpha | '_') >> *(alnum | '_')] - distinct_keywords;

        BOOST_SPIRIT_DEFINE(identifier);

        const auto obj_fun_term =
            -(double_[add_obj_fun_coeff]) >> identifier[add_obj_fun_var];

        const auto obj_fun_expr =
            -(identifier >> ':')[set_obj_label] >> -(char_('+') | char_('-')[neg_var_coeff]) >> obj_fun_term[term_parsed] >> *((char_('+') | char_('-')[neg_var_coeff]) > obj_fun_term[term_parsed]);

        const auto constr_term =
            -(double_[add_constr_coeff]) >> identifier[add_constr_var];

        const auto constraint_expr =
            -(identifier >> ':')[set_constr_label] >>
            (-(char_('+') | char_('-')[neg_var_coeff]) >> constr_term[term_parsed] >> *((char_('+') | char_('-')[neg_var_coeff]) > constr_term[term_parsed]) >> (char_("<>="))[add_constr_type] >> -char_('=') >> double_[add_constr_rhs]);

        const auto bound_expr =
            ((double_[add_constr_coeff] | ('-' >> (no_case["Inf"] | no_case["Infinity"])[add_inf_lb])) >> lit("<=") >> identifier[add_lb_constr] >> -(lit("<=") >> double_[add_ub_constr])) |
            (identifier[store_ub_var_name] >> lit("<=") >> (double_[add_ub_constr] | (-char_('+') >> (no_case["Inf"] | no_case["Infinity"])[add_inf_ub])));

        x3::rule<class lp_rules, ParserState> const lp_rules("lp_rules");
        const auto lp_rules_def =
            skip(space)[(no_case["Minimize"][set_minimize] | no_case["Maximize"][set_maximize]) >> -char_(':') >> obj_fun_expr >> no_case["Subject To"] >> -char_(':') >> +constraint_expr >> -(no_case["Bounds"] >> -char_(':') >> +bound_expr) >> no_case["End"] >> -char_('.') >> eoi];

        BOOST_SPIRIT_DEFINE(lp_rules);
    }

    std::string remove_comments(const std::string &input)
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

    bool run_parser_lp(const std::string &input, LinearProgram &lp)
    {
        const std::string stripped_input = remove_comments(input);
        auto iter = stripped_input.begin();
        const auto end = stripped_input.end();

        ParserState parser_state{lp};

        bool result = parse(iter, end, parser::lp_rules, parser_state);

        if (result && iter == end)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}
