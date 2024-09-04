#include <gtest/gtest.h>

#include "../src/LinearProgram.h"

class UnitTestsBase : public ::testing::Test
{
protected:
    void SetUp() override { }
    void TearDown() override { }
};

TEST_F(UnitTestsBase, ConstraintTests1)
{
    simplex::Constraint constraint('>', 9);
    constraint.add_term("x1");
    constraint.add_term("y23", -3.21);
    constraint.add_term("new_var", 0.0);
    constraint.add_term("new_var1", -1.0);

    EXPECT_EQ(constraint.rhs, 9);
    EXPECT_EQ(constraint.type, '>');
    ASSERT_TRUE(constraint.has_variable("x1"));
    ASSERT_TRUE(constraint.has_variable("y23"));
    ASSERT_FALSE(constraint.has_variable("new_var"));
    ASSERT_TRUE(constraint.has_variable("new_var1"));

    EXPECT_EQ(constraint.get_coefficient("x1"), 1.0);
    EXPECT_EQ(constraint.get_coefficient("new_var"), std::nullopt);
    EXPECT_EQ(constraint.get_coefficient("y23"), -3.21);

    constraint.remove_term("y23");
    ASSERT_FALSE(constraint.has_variable("y23"));
    EXPECT_EQ(constraint.get_coefficient("y23"), std::nullopt);
}

TEST_F(UnitTestsBase, LinearProgramTests1)
{
    simplex::LinearProgram lp;

	lp.add_variable("x1");
    lp.add_variable("x2", 1.0);
    lp.add_variable("x3", -2.0);

    auto c1 = simplex::Constraint('<', 5).add_term("x1").add_term("x2", 2.0);

    lp.add_constraint("C1", std::move(c1));

	lp.initialize_tableau();
    lp.print_tableau();

    EXPECT_THROW(std::ignore = lp.get_constraint("C2"), std::string);
    const auto & test_c1 = lp.get_constraint("C1");
    ASSERT_TRUE(test_c1.has_variable("x2"));
    EXPECT_EQ(test_c1.get_coefficient("x2").value(), 2.0);
}
