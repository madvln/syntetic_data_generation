#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
using namespace pde_solvers;

#include <iostream>
#include <fstream>
#include <filesystem>
#include <time.h>
#include <algorithm>

/// @brief ���������� �������� ������ � ������� TestBundle.TestName
inline std::string get_test_string() {
    auto test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    auto test_string = std::string(test_info->test_case_name()) + "." + std::string(test_info->name());
    return test_string;
}

/// @brief ������� ���� ��� ����� � ����� 
/// ���������� ���� � ��������� �����
inline std::string prepare_test_folder()
{
    std::string path = std::string("../testing_out/") + get_test_string() + "/";
    std::filesystem::create_directories(path);
    return path;
}

#include "../research/test_data_generation.h"



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef _WIN32
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
    std::locale::global(std::locale("en_US.UTF-8"));
    int res = RUN_ALL_TESTS();
    return res;
}

