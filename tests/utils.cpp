#include <iostream>
#include <sstream>
#include <string>

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
