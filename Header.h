#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <random>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4146)   // unary minus on unsigned (GMP internal)
#pragma warning(disable: 26812)  // unscoped enum in GMP headers (gmp_randalg_t)
#endif
#include <gmpxx.h>
#ifdef _MSC_VER
#pragma warning(pop)
#endif

// ======================== Modular Arithmetic ========================

// --- mod: always returns a non-negative result in [0, modulus) ---
inline mpz_class mod(const mpz_class& value, const mpz_class& modulus) {
    mpz_class result = value % modulus;
    if (result < 0) result += modulus;
    return result;
}

inline long long mod(long long value, long long modulus) {
    long long result = value % modulus;
    if (result < 0) result += modulus;
    return result;
}

// --- modInverse: modular multiplicative inverse (a^{-1} mod p) ---
inline mpz_class modInverse(const mpz_class& a, const mpz_class& p) {
    mpz_class result;
    if (mpz_invert(result.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t()) == 0)
        throw std::runtime_error("modInverse: inverse does not exist");
    return result;
}

inline long long modInverse(long long a, long long p) {
    a = ((a % p) + p) % p;
    if (a == 0)
        throw std::runtime_error("modInverse: inverse of 0 does not exist");
    long long old_r = a, r = p;
    long long old_s = 1, s = 0;
    while (r != 0) {
        long long q = old_r / r;
        long long temp_r = r; r = old_r - q * r; old_r = temp_r;
        long long temp_s = s; s = old_s - q * s; old_s = temp_s;
    }
    if (old_r != 1)
        throw std::runtime_error("modInverse: inverse does not exist (gcd != 1)");
    return ((old_s % p) + p) % p;
}

// --- modPow: modular exponentiation (base^exp mod p) ---
inline mpz_class modPow(const mpz_class& base, const mpz_class& exp, const mpz_class& p) {
    mpz_class result;
    mpz_powm(result.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    return result;
}

inline long long modPow(long long base, long long exp, long long p) {
    base = ((base % p) + p) % p;
    long long result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) result = (result * base) % p;
        base = (base * base) % p;
        exp /= 2;
    }
    return result;
}

// Include curve and point class definitions
#include "EllipticCurve.h"
#include "EllipticCurvePoint.h"