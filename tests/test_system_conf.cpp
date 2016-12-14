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

void print_test_fail(const string err_string, const string file, const int line)
{
    cerr << "  [FAIL] " << file << ':' << line << ' ' << err_string << '\n';
    test_failed = true;
}

// Compare to input objects and see that they are equivalent. They may
// be of different types, in which case a static_cast must be possible
// to perform.
template<typename T1, typename T2>
void assert_eq(const T1 recv, const T2 expt, const string err_string, const string file, const int line)
{
    ostringstream oss;
    oss << err_string << " (expected: " << expt << ", received: " << recv << ')';
    if (recv != static_cast<T1>(expt)) {
        print_test_fail(oss.str(), file, line);
    }
}

// Compare two input arrays and assert that they are identical.
// The objects must be of identical types implementing size(), cbegin()
// and cend() methods for iteration over the values. Their size and
// values are both controlled.
template<typename T>
void assert_eq_vec(const T recv, const T expt, const string err_string, const string file, const int line)
{
    bool not_equal = false;

    if (recv.size() != expt.size()) {
        not_equal = true;
    }
    else {
        auto recvit = recv.cbegin();

        for (auto e : expt) {
            if (e != *recvit++) {
                not_equal = true;
                break;
            }
        }
    }

    if (not_equal) {
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
    SystemConf conf {3};
    ASSERT_EQ(conf.num_atoms(), 0, "number of atoms is not zero-initialized");
    ASSERT_EQ(conf.xs.size(), 0, "number of position elements is not 0");
    ASSERT_EQ(conf.vs.size(), 0, "number of velocity elements is not 0");
    ASSERT_EQ(conf.fs.size(), 0, "number of force elements is not 0");
    ASSERT_EQ(conf.xs.capacity(), 9, "capacity of vector is not 3*3 elems as input");

    array<double, 3> box { 0.0, 0.0, 0.0 };
    ASSERT_EQ_VEC(conf.box, box, "box is not zero-initialized");

    ASSERT_EQ(conf.title, "", "title is not initialized as empty");

    finalize_test();
}

void test_system_add_atom()
{
    SystemConf conf {3};
    auto num = conf.add_atom(1.0, 2.0, 3.0);
    ASSERT_EQ(num, 1, "new atom number was not correctly returned after add");
    num = conf.add_atom(1.0, 2.0, 3.0);
    ASSERT_EQ(num, 2, "new atom number was not correctly returned after add");

    ASSERT_EQ(conf.num_atoms(), 2, "atom number was not updated after add");
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

void test_system_set_box()
{
    SystemConf conf {1};
    conf.set_box(1.0, 2.0, 3.0);
    array<double, 3> box {1.0, 2.0, 3.0};
    ASSERT_EQ_VEC(conf.box, box, "box was not set correctly");

    finalize_test();
}

void test_system_read_grofile()
{
    const string path = "tests/include/grofile_small.gro";
    const auto conf = read_conf_from_grofile(path);

    ASSERT_EQ(conf.title, "small test grofile", "incorrect title read");

    // The file has three atoms and only the positions should be read
    ASSERT_EQ(conf.num_atoms(), 3, "incorrect number of atoms read");

    vector<double> xs {
        0.129, 0.079, 0.464,
        0.119, 0.061, 0.314,
        0.109, 0.042, 0.165
    };
    ASSERT_EQ_VEC(conf.xs, xs, "positions were not read correctly");

    vector<double> zs (9, 0.0);
    ASSERT_EQ_VEC(conf.vs, zs, "velocities are not zero-initialized");
    ASSERT_EQ_VEC(conf.fs, zs, "forces are not zero-initialized");

    array<double, 3> box { 216.00000, 4.67650, 110.00000 };
    ASSERT_EQ_VEC(conf.box, box, "incorrect box size read");

    finalize_test();
}

int main(int argc, char* argv[])
{
    cout << "Running tests in '" << __FILE__ << "'.\n";

    test_system_init();
    test_system_add_atom();
    test_system_set_box();
    test_system_read_grofile();

    if (num_failed == 0) {
        cout << "All tests passed.\n\n";
    }
    else {
        cout << num_failed << " test(s) failed.\n\n";
    }

    return 0;
}
