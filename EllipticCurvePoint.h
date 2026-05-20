#pragma once

// Elliptic curve point in projective coordinates (X, Y, Z)

template <typename T>
class EllipticCurvePoint {
public:
    T X, Y, Z;
    const EllipticCurve<T>* curve;

    // --- Constructors ---
    EllipticCurvePoint(const T& X, const T& Y, const T& Z, const EllipticCurve<T>* curve)
        : X(X), Y(Y), Z(Z), curve(curve) {
    }

    static EllipticCurvePoint infinity(const EllipticCurve<T>* curve) {
        return EllipticCurvePoint(T(0), T(1), T(0), curve);
    }

    static EllipticCurvePoint fromAffine(const T& x, const T& y, const EllipticCurve<T>* curve) {
        EllipticCurvePoint pt(x, y, T(1), curve);
        if (!pt.isOnCurve())
            throw std::runtime_error("fromAffine: point is NOT on the curve");
        return pt;
    }

    // --- Checks ---
    bool isInfinity() const {
        return Z == 0;
    }

    void assertOnCurve(const std::string& context) const {
        if (isInfinity()) return;
        const T& p = curve->p;
        const T& a = curve->a;
        const T& b = curve->b;

        T Y2 = mod(Y * Y, p);
        T lhs = mod(Y2 * Z, p);
        T X2 = mod(X * X, p);
        T X3 = mod(X2 * X, p);
        T Z2 = mod(Z * Z, p);
        T Z3 = mod(Z2 * Z, p);
        T rhs = mod(X3 + mod(mod(a * X, p) * Z2, p) + mod(b * Z3, p), p);

        if (lhs != rhs)
            throw std::runtime_error(context + ": point is NOT on the curve");
    }

    bool isOnCurve() const {
        if (isInfinity()) return true;
        const T& p = curve->p;
        const T& a = curve->a;
        const T& b = curve->b;

        T Y2 = mod(Y * Y, p);
        T lhs = mod(Y2 * Z, p);
        T X2 = mod(X * X, p);
        T X3 = mod(X2 * X, p);
        T Z2 = mod(Z * Z, p);
        T Z3 = mod(Z2 * Z, p);
        T aXZ2 = mod(mod(a * X, p) * Z2, p);
        T bZ3 = mod(b * Z3, p);
        T rhs = mod(X3 + aXZ2 + bZ3, p);

        return lhs == rhs;
    }

    // --- Conversion ---
    std::pair<T, T> toAffine() const {
        if (isInfinity())
            throw std::runtime_error("toAffine: point at infinity has no affine coordinates");
        const T& p = curve->p;
        T zInv = modInverse(Z, p);
        return { mod(X * zInv, p), mod(Y * zInv, p) };
    }

    // --- Arithmetic ---
    EllipticCurvePoint pointDouble() const {
        if (isInfinity()) return infinity(curve);
        assertOnCurve("pointDouble");

        const T& p = curve->p;
        const T& a = curve->a;

        if (mod(Y, p) == 0) return infinity(curve);

        T W = mod(mod(a * mod(Z * Z, p), p) + mod(3 * mod(X * X, p), p), p);
        T S = mod(Y * Z, p);
        T B = mod(mod(X * Y, p) * S, p);
        T H = mod(mod(W * W, p) - mod(8 * B, p), p);
        T Xr = mod(mod(2 * H, p) * S, p);
        T Yr = mod(mod(W * mod(4 * B - H, p), p) - mod(mod(8 * mod(Y * Y, p), p) * mod(S * S, p), p), p);
        T Zr = mod(8 * mod(mod(S * S, p) * S, p), p);

        return EllipticCurvePoint(Xr, Yr, Zr, curve);
    }

    EllipticCurvePoint pointAdd(const EllipticCurvePoint& other) const {
        const T& p = curve->p;

        if (isInfinity()) return other;
        if (other.isInfinity()) return *this;

        assertOnCurve("pointAdd (P)");
        other.assertOnCurve("pointAdd (Q)");

        T U1 = mod(other.Y * Z, p);
        T U2 = mod(Y * other.Z, p);
        T V1 = mod(other.X * Z, p);
        T V2 = mod(X * other.Z, p);

        if (V1 == V2) {
            if (U1 != U2) return infinity(curve);
            else return pointDouble();
        }

        T U = mod(U1 - U2, p);
        T V = mod(V1 - V2, p);
        T W = mod(Z * other.Z, p);
        T Vsq = mod(V * V, p);
        T Vcb = mod(Vsq * V, p);
        T Usq = mod(U * U, p);
        T VsqV2 = mod(Vsq * V2, p);
        T A = mod(mod(Usq * W, p) - Vcb - mod(2 * VsqV2, p), p);
        T X3 = mod(V * A, p);
        T Y3 = mod(mod(U * mod(VsqV2 - A, p), p) - mod(Vcb * U2, p), p);
        T Z3 = mod(Vcb * W, p);

        return EllipticCurvePoint(X3, Y3, Z3, curve);
    }

    // --- Scalar Multiplication ---
    EllipticCurvePoint scalarMul(const T& k) const {
        if (k == 0 || isInfinity()) return infinity(curve);
        assertOnCurve("scalarMul");

        EllipticCurvePoint res = infinity(curve);
        EllipticCurvePoint temp = *this;
        std::vector<int> bits = getBits(k);

        for (size_t i = 0; i < bits.size(); i++) {
            if (bits[i] == 1) res = res.pointAdd(temp);
            temp = temp.pointDouble();
        }
        return res;
    }

private:
    static std::vector<int> getBits(const T& k) {
        std::vector<int> bits;
        T val = k;
        while (val > 0) {
            T remainder = val % 2;
            bits.push_back(remainder == 0 ? 0 : 1);
            val /= 2;
        }
        return bits;
    }

public:
    // --- Output ---
    static std::string toString(const T& val) {
        std::ostringstream oss;
        oss << val;
        return oss.str();
    }

    void print() const {
        if (isInfinity()) {
            std::cout << "O_E (0 : 1 : 0)" << std::endl;
        }
        else {
            const T& p = curve->p;
            T zInv = modInverse(Z, p);
            T ax = mod(X * zInv, p);
            T ay = mod(Y * zInv, p);
            std::cout << "(" << toString(ax) << ", " << toString(ay) << ")" << std::endl;
        }
    }

    // --- Comparison ---
    bool equals(const EllipticCurvePoint& other) const {
        if (isInfinity() && other.isInfinity()) return true;
        if (isInfinity() || other.isInfinity()) return false;
        const T& p = curve->p;
        return mod(X * other.Z, p) == mod(other.X * Z, p) &&
            mod(Y * other.Z, p) == mod(other.Y * Z, p);
    }
};
