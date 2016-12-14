#include <iostream>
#include <sstream>
#include <string>

#include "src/conf.h"

#define ASSERT(cond, err_string) if (!cond) { print_test_fail(err_string, __func__, __LINE__); }
#define ASSERT_EQ(a, b, err_string) assert_eq(a, b, err_string, __func__, __LINE__);
#define ASSERT_EQ_VEC(a, b, err_string) assert_eq_vec(a, b, err_string, __func__, __LINE__);

using namespace std;

static int num_failed = 0;
static bool test_failed = false;

void print_test_fail(string err_string, string file, int line)
{
    cerr << "  [FAIL] " << file << ':' << line << ' ' << err_string << '\n';
    test_failed = true;
}

template<typename T1, typename T2>
void assert_eq(T1 recv, T2 expt, string err_string, string file, int line)
{
    ostringstream oss;
    oss << err_string << " (expected: " << expt << ", received: " << recv << ')';
    if (recv != static_cast<T1>(expt)) {
        print_test_fail(oss.str(), file, line);
    }
}

template<typename T>
void assert_eq_vec(vector<T> recv, vector<T> expt, string err_string, string file, int line)
{
    if (recv != expt) {
        ostringstream oss;
        oss << err_string << " (expected: [";
        for (auto v : expt) {
            oss << v << ' ';
        }
        oss << "], received: [";
        for (auto v : recv) {
            oss << v << ' ';
        }
        oss << "])";

        print_test_fail(oss.str(), file, line);
    }
}

void finalize_test()
{
    if (test_failed) {
        ++num_failed;
        test_failed = false;
    }
}

void test_system_init()
{
    system_conf conf {3};
    ASSERT_EQ(conf.num_atoms, 0, "number of atoms is not zero-initialized");
    ASSERT_EQ(conf.xs.size(), 0, "number of position elements is not 0");
    ASSERT_EQ(conf.vs.size(), 0, "number of velocity elements is not 0");
    ASSERT_EQ(conf.fs.size(), 0, "number of force elements is not 0");
    ASSERT_EQ(conf.xs.capacity(), 9, "capacity of vector is not 3*3 elems as input");

    finalize_test();
}

void test_system_add_atom()
{
    system_conf conf {3};
    auto num = conf.add_atom(1.0, 2.0, 3.0);
    ASSERT_EQ(num, 1, "new atom number was not correctly returned after add");
    num = conf.add_atom(1.0, 2.0, 3.0);
    ASSERT_EQ(num, 2, "new atom number was not correctly returned after add");

    ASSERT_EQ(conf.num_atoms, 2, "atom number was not updated after add");
    ASSERT_EQ(conf.xs.size(), 6, "wrong number of position elements");
    ASSERT_EQ(conf.vs.size(), 6, "wrong number of velocity elements");
    ASSERT_EQ(conf.fs.size(), 6, "wrong number of force elements");

    vector<double> xs {1.0, 2.0, 3.0, 1.0, 2.0, 3.0};
    vector<double> zs {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    ASSERT_EQ_VEC(conf.xs, xs, "added atom positions are not those input");
    ASSERT_EQ_VEC(conf.vs, zs, "added atom velocities are not zero-initialized");
    ASSERT_EQ_VEC(conf.fs, zs, "added atom forces are not zero-initialized");

    finalize_test();
}

int main(int argc, char* argv[])
{
    cout << "Running tests in '" << __FILE__ << "'.\n";

    test_system_init();
    test_system_add_atom();

    if (num_failed == 0) {
        cout << "All tests passed.\n\n";
    }
    else {
        cout << num_failed << " test(s) failed.\n\n";
    }

    return 0;
}
