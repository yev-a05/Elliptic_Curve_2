#pragma once

// Elliptic curve in short Weierstrass form: y^2 = x^3 + ax + b (mod p)

template <typename T>
class EllipticCurve {
public:
    T a, b, p;
    T n;           // curve order (#E(F_p), including O_E)
    bool verbose;  // step-by-step output

    EllipticCurve(const T& a, const T& b, const T& p, const T& n, bool verbose = false)
        : a(a), b(b), p(p), n(n), verbose(verbose) {
    }

    static std::string toString(const T& val) {
        std::ostringstream oss;
        oss << val;
        return oss.str();
    }

    // Check that the curve is non-singular: Delta = -16(4a^3 + 27b^2) != 0 (mod p).
    bool isNonSingular() const {
        T disc_core = mod(T(4) * a * a * a + T(27) * b * b, p);
        if (verbose) {
            std::cout << "[isNonSingular] Delta check: " << toString(disc_core)
                << (disc_core != T(0) ? " != 0 -> non-singular" : " == 0 -> SINGULAR!")
                << std::endl;
        }
        return disc_core != T(0);
    }

    // Print the curve equation
    void print() const {
        std::cout << "Elliptic Curve: y^2 = x^3";
        if (a != T(0)) std::cout << " + " << toString(a) << "*x";
        if (b != T(0)) std::cout << " + " << toString(b);
        std::cout << "  (mod " << toString(p) << ")" << std::endl;
    }
};
