// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include "parser.h"
#include "utils.h"

#include <boost/spirit/home/x3.hpp>

#include <iostream>

namespace x3 = boost::spirit::x3;
namespace ascii = boost::spirit::x3::ascii;

namespace parser 
{
    using x3::eps;
    using x3::lit;
    using x3::lexeme;
    using x3::double_;
    using x3::eol;
    using x3::eoi;
    using x3::skip;
    using x3::no_case;

    using ascii::char_;
    using ascii::alpha;
    using ascii::alnum;
    using ascii::space;

    // Semantic actions definitions
    auto set_minimize       = [](auto & ctx) { _val(ctx).set_sense('m'); };
    auto set_maximize       = [](auto & ctx) { _val(ctx).set_sense('M'); };
    auto set_obj_label      = [](auto & ctx) { _val(ctx).set_objective_label(_attr(ctx)); };
    auto add_obj_fun_coeff  = [](auto & ctx) { _val(ctx).set_var_coeff_cache(_attr(ctx)); };
    auto add_obj_fun_var    = [](auto & ctx) 
    { 
       _val(ctx).add_variable(_attr(ctx)); 
       _val(ctx).objective_name_coeff[_attr(ctx)] = 
            _val(ctx).get_sign_cache() 
            ? -_val(ctx).get_var_coeff_cache()
            :  _val(ctx).get_var_coeff_cache();  
     };
    auto neg_var_coeff      = [](auto & ctx) { _val(ctx).set_sign_cache(true); };
    auto term_parsed        = [](auto & ctx) 
    { 
        _val(ctx).set_sign_cache(false);
        _val(ctx).set_var_coeff_cache(1); 
    };
    auto add_obj_fun_const  = [](auto & ctx) { _val(ctx).set_obj_value_shift(_attr(ctx)); };

    auto set_constr_label   = [](auto & ctx) 
    {
        // create labeled constraint
        auto constr_ptr = make_shared<Constraint>();
        _val(ctx).set_label_cache(_attr(ctx)); 
        _val(ctx).constraints[_attr(ctx)] = constr_ptr;
    };
    auto add_constr_coeff   = [](auto & ctx) { _val(ctx).set_var_coeff_cache(_attr(ctx)); };
    auto add_constr_var     = [](auto & ctx) 
    {
        auto constr_label = _val(ctx).get_label_cache();
        if (constr_label == "") {
            // no label was defined by user; create new constraint with default label
            auto constr_ptr = make_shared<Constraint>();
            constr_label = string("CONSTR") + tostr<size_t>(_val(ctx).constraints.size());
            _val(ctx).set_label_cache(constr_label);
            _val(ctx).constraints[constr_label] = constr_ptr;
        }

        if (!_val(ctx).has_variable(_attr(ctx))) {
            // adding new variable via constraint
            _val(ctx).add_variable(_attr(ctx));
        }

        auto coeff = _val(ctx).get_var_coeff_cache();
        auto sign  = _val(ctx).get_sign_cache() ? -1 : 1;
        _val(ctx).constraints[constr_label]->name_coeff[_attr(ctx)] = sign * coeff;
    };

    auto add_constr_type = [](auto & ctx) 
    { 
        _val(ctx).constraints[_val(ctx).get_label_cache()]->type = _attr(ctx); 
    };
    auto add_constr_rhs  = [](auto & ctx) 
    { 
        _val(ctx).constraints[_val(ctx).get_label_cache()]->rhs = _attr(ctx);
        _val(ctx).set_label_cache("");
    };

    auto add_lb_constr   = [](auto & ctx) 
    { 
        auto lb = _val(ctx).get_var_coeff_cache();  // now it contains LB
        auto var_name = _attr(ctx);
        _val(ctx).set_label_cache(var_name);  // store var_name, in case we need to add UB
        if (_isfloatzero(lb)) return;
        if (lb < 0.0) {
            // negative lower bound: substitute z = x + lb, z >= 0
            _val(ctx).add_shift(var_name, lb);
        }
        else {
            // add ordinary constraint x >= lb
            auto constr_label = string("LBND") + tostr<size_t>(_val(ctx).constraints.size());
            auto constr_ptr = make_shared<Constraint>('>', lb);
            _val(ctx).constraints[constr_label] = constr_ptr;
            _val(ctx).constraints[constr_label]->name_coeff[var_name] = 1.0;
        }
    };
    auto store_ub_var_name = [](auto & ctx) { _val(ctx).set_label_cache(_attr(ctx)); };
    auto add_ub_constr     = [](auto & ctx)
    {
        auto var_name = _val(ctx).get_label_cache();    // now it contains var_name
        auto ub = _attr(ctx);
        char _type = '<';
        if (ub < 0.0) {
            // negative upper bound: substitute z = -x
            _val(ctx).add_shift(var_name, 0.0);
            ub = -ub;
            _type = '>';
        }
        // add ordinary constraint x <= ub (or -x >= -ub)
        auto constr_ptr = make_shared<Constraint>(_type, ub);
        auto constr_label = string("UBND") + tostr<size_t>(_val(ctx).constraints.size());
        _val(ctx).constraints[constr_label] = constr_ptr;
        _val(ctx).constraints[constr_label]->name_coeff[var_name] = 1.0;
    };

    auto print_attr = [](auto & ctx) { std::cout << _attr(ctx) << std::endl; };

    // Grammar rules definitions
    // TODO: handle properly the constant objective function, e.g., obj: 0

    struct keywords_t : x3::symbols<x3::unused_type> {
        keywords_t() {
            add("Subject")
               ("Bounds")
               ("End");
        }
    } const keywords;

    auto const distinct_keywords = lexeme[no_case[keywords] >> !(alnum | '_')];

    x3::rule<struct identifier_tag, std::string> const identifier ("identifier");
    auto const identifier_def = lexeme[+(alpha | '_') >> *(alnum | '_')] - distinct_keywords;

    BOOST_SPIRIT_DEFINE(identifier);

    auto const obj_fun_term =
           -(double_ [add_obj_fun_coeff])
        >> identifier [add_obj_fun_var]
    ;

    auto const obj_fun_expr =
           -(identifier [set_obj_label] >> ':') 
        >> -(char_('+') | char_('-')[neg_var_coeff]) >> obj_fun_term [term_parsed]
        >> *((char_('+') | char_('-')[neg_var_coeff]) > obj_fun_term [term_parsed])
    ;

    auto const constr_term =
           -(double_ [add_constr_coeff])
        >> identifier [add_constr_var]
    ;

    auto const constraint_expr = 
        -(identifier >> ':') [set_constr_label]
        >> 
        (
               -(char_('+') | char_('-')[neg_var_coeff]) >> constr_term [term_parsed]
            >> *((char_('+') | char_('-')[neg_var_coeff]) > constr_term [term_parsed])
            >> (char_("<>=")) [add_constr_type] >> -char_('=')
            >> double_ [add_constr_rhs]
        )
    ;

    auto const bound_expr =
        (double_[add_constr_coeff] 
            >> lit("<=") >> identifier [add_lb_constr] 
            >> -(lit("<=") >> double_[add_ub_constr]))
        |
        (identifier [store_ub_var_name] 
            >> lit("<=") >> double_[add_ub_constr])
    ;

    x3::rule<class lp_rules, LinearProgram> const lp_rules ("lp_rules");
    auto const lp_rules_def =
        skip(space)[
               (no_case["Minimize"] [set_minimize] | no_case["Maximize"] [set_maximize])
            >> -char_(':')
            >> obj_fun_expr
            >> no_case["Subject To"] >> -char_(':')
            >> +constraint_expr
            >> -(no_case["Bounds"] >> -char_(':') >> +bound_expr)
            >> no_case["End"] >> -char_('.')
            >> eoi
        ]
    ;

    BOOST_SPIRIT_DEFINE(lp_rules);
}

std::string remove_comments(const std::string & input)
{
    std::string result("");
    size_t pos_start = 0, pos_cur = 0;

    while (true) 
    {
        pos_cur = input.find_first_of("/\\#", pos_start);
        if (pos_cur == string::npos) {
            pos_cur = input.length();
            break;
        }
        result += input.substr(pos_start, pos_cur - pos_start);
        pos_start = input.find_first_of('\n', pos_cur) + 1;
    }
    result += input.substr(pos_start, pos_cur - pos_start);

    return result;
}

bool run_parser_lp(const std::string & input, LinearProgram & lp)
{
    std::string stripped_input = remove_comments(input);
    typedef std::string::const_iterator iterator_type;
    iterator_type iter = stripped_input.begin();
    iterator_type const end = stripped_input.end();
    
    bool r = parse(iter, end, parser::lp_rules, lp);

    if (r && iter == end) {
        std::cout << "Parser success." << endl;
        return true;
    }
    else
    {
        std::cout << "Parser fail." << endl;
        return false;
    }
}

