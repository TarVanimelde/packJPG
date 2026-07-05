// test_framework.h - a minimal, single-header C++ test framework for packJPG.
//
// Why home-grown instead of googletest/doctest/Catch2?
//   The prompt's hard constraint is "single-header, no new build system". This
//   header satisfies that in ~150 lines with zero external downloads, which also
//   matters because the build sandbox has no package-manager network access to
//   fetch a vendored framework. It provides the small slice of functionality the
//   suite needs: auto-registered test cases, CHECK/REQUIRE assertions, and a
//   runner. If the project later wants richer tooling, swapping in doctest is a
//   drop-in change (same TEST_CASE / CHECK / REQUIRE spelling).
//
// Usage:
//   #include "test_framework.h"
//   TEST_CASE("bitreader reads bits") { CHECK(x == y); REQUIRE(ptr != nullptr); }
//   // exactly one translation unit must #define TEST_FRAMEWORK_MAIN before include.

#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <cstdio>
#include <exception>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace tf {

// Thrown by REQUIRE on failure to abort the current test case immediately.
struct RequireFailure : std::exception {};

struct TestCase {
    std::string name;
    std::function<void()> fn;
};

// Global registry (function-local static avoids static-init-order issues).
inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> cases;
    return cases;
}

// Per-run counters.
inline int& checks_failed() { static int n = 0; return n; }
inline int& checks_total()  { static int n = 0; return n; }
// Set true by any failing CHECK/REQUIRE within the current case.
inline bool& current_failed() { static bool b = false; return b; }

struct Registrar {
    Registrar(const std::string& name, std::function<void()> fn) {
        registry().push_back(TestCase{name, std::move(fn)});
    }
};

// Records a single assertion outcome; prints on failure.
inline void report(bool ok, const char* expr, const char* file, int line,
                   const std::string& extra = std::string()) {
    ++checks_total();
    if (!ok) {
        ++checks_failed();
        current_failed() = true;
        std::printf("    FAILED: %s\n      at %s:%d\n", expr, file, line);
        if (!extra.empty()) std::printf("      %s\n", extra.c_str());
    }
}

template <typename A, typename B>
inline std::string eq_detail(const A& a, const B& b) {
    std::ostringstream os;
    os << "expected equal: lhs=" << a << " rhs=" << b;
    return os.str();
}

inline int run_all() {
    int failed_cases = 0;
    std::printf("Running %zu test case(s)...\n", registry().size());
    for (auto& tc : registry()) {
        current_failed() = false;
        int fails_before = checks_failed();
        try {
            tc.fn();
        } catch (const RequireFailure&) {
            // REQUIRE already reported; case is marked failed.
        } catch (const std::exception& e) {
            current_failed() = true;
            std::printf("    FAILED: uncaught std::exception: %s\n", e.what());
        } catch (...) {
            current_failed() = true;
            std::printf("    FAILED: uncaught non-standard exception\n");
        }
        bool ok = !current_failed();
        std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", tc.name.c_str());
        if (!ok) ++failed_cases;
        (void)fails_before;
    }
    std::printf("\n%d/%zu cases passed, %d/%d checks passed.\n",
                (int)registry().size() - failed_cases, registry().size(),
                checks_total() - checks_failed(), checks_total());
    return failed_cases == 0 ? 0 : 1;
}

} // namespace tf

// ---- Macros -------------------------------------------------------------

#define TF_CONCAT_INNER(a, b) a##b
#define TF_CONCAT(a, b) TF_CONCAT_INNER(a, b)

#define TEST_CASE(NAME)                                                       \
    static void TF_CONCAT(tf_test_, __LINE__)();                             \
    static ::tf::Registrar TF_CONCAT(tf_reg_, __LINE__)(                     \
        NAME, &TF_CONCAT(tf_test_, __LINE__));                              \
    static void TF_CONCAT(tf_test_, __LINE__)()

// Non-fatal: records failure, continues the case.
#define CHECK(expr) ::tf::report((expr), #expr, __FILE__, __LINE__)
#define CHECK_FALSE(expr) ::tf::report(!(expr), "!(" #expr ")", __FILE__, __LINE__)
// NOTE: arguments are evaluated exactly ONCE (captured into locals) so that
// side-effecting expressions like reader.read_byte() behave correctly.
#define CHECK_EQ(a, b)                                                        \
    do {                                                                      \
        auto tf_a_ = (a);                                                     \
        auto tf_b_ = (b);                                                     \
        ::tf::report(tf_a_ == tf_b_, #a " == " #b, __FILE__, __LINE__,        \
                     ::tf::eq_detail(tf_a_, tf_b_));                          \
    } while (0)

// Fatal: records failure and aborts the current case.
#define REQUIRE(expr)                                                         \
    do {                                                                      \
        bool tf_ok_ = (expr);                                                 \
        ::tf::report(tf_ok_, #expr, __FILE__, __LINE__);                     \
        if (!tf_ok_) throw ::tf::RequireFailure{};                            \
    } while (0)

#ifdef TEST_FRAMEWORK_MAIN
int main() { return ::tf::run_all(); }
#endif

#endif // TEST_FRAMEWORK_H
