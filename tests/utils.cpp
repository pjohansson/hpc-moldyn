#include <iostream>
#include <sstream>
#include <string>

//#define ASSERT(cond, err_string) if (!cond) { print_test_fail(err_string, __func__, __LINE__); }

// Compare to input objects and see that they are equivalent. They may
// be of different types, in which case a static_cast must be possible
// to perform.
#define ASSERT_EQ(a, b, err_string) assert_eq(a, b, err_string, __func__, __LINE__);

// Compare two input arrays and assert that they are identical.
// The objects must be of identical types implementing size(), cbegin()
// and cend() methods for iteration over the values. Their size and
// values are both controlled.
#define ASSERT_EQ_VEC(a, b, err_string) assert_eq_vec(a, b, err_string, __func__, __LINE__);

// Macro which adds a main function which runs the input functions.
#define RUN_TESTS(...) \
    int main(int argc, char* argv[]) { \
        std::cout << "Running tests in '" << __FILE__ << "'.\n"; \
        __VA_ARGS__ \
        if (num_failed == 0) { \
            cout << "All tests passed.\n\n"; \
        } \
        else if (num_failed == 1) { \
            cout << num_failed << " test failed.\n\n"; \
        } \
        else { \
            cout << num_failed << " tests failed.\n\n"; \
        } \
        return num_failed; \
    }

// Macro which wraps a test function with a void declaration
// and updates the test status afterwards.
#define ADD_TEST(name, ...) \
    void name() { \
        __VA_ARGS__ \
        finalize_test(); \
    }

int num_failed = 0;
bool test_failed = false;

void print_test_fail(const std::string err_string, const std::string file, const int line)
{
    std::cerr << "  [FAIL] " << file << ':' << line << ' ' << err_string << '\n';
    test_failed = true;
}

template<typename T1, typename T2>
void assert_eq(const T1 recv, const T2 expt, const std::string err_string, const std::string file, const int line)
{
    std::ostringstream oss;
    oss << err_string << " (expected: " << expt << ", received: " << recv << ')';
    if (recv != static_cast<T1>(expt)) {
        print_test_fail(oss.str(), file, line);
    }
}

template<typename T>
void assert_eq_vec(const T recv, const T expt, const std::string err_string, const std::string file, const int line)
{
    bool not_equal = false;

    if (recv.size() != expt.size()) {
        not_equal = true;
    }
    else {
        auto r = recv.cbegin();
        for (auto e : expt) {
            if (e != *r++) {
                not_equal = true;
                break;
            }
        }
    }

    if (not_equal) {
        std::ostringstream oss;
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
