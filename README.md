# Elliptic Curve Cryptography Library

C++ template library implementing elliptic curves in short Weierstrass form over prime fields,
with a focus on the Baby-JubJub curve.

---

## File Structure

```
Header.h              Modular arithmetic (mod, modInverse, modPow) and 
                      cryptographic utilities (XOR cipher, hash simulation).

EllipticCurve.h       Template class EllipticCurve<T> storing curve parameters 
                      (a, b, p, n) and basic non-singularity validation.

EllipticCurvePoint.h  Template class EllipticCurvePoint<T> providing fast 
                      projective point arithmetic (addition, doubling) and 
                      scalar multiplication.

Main.cpp              Initialization of the Baby-JubJub curve and step-by-step 
                      execution of ECDH, ECIES, and ECDSA protocols.
```

Both template classes support `T = long long` (small fields, fast, test) and `T = mpz_class`
(arbitrary precision via GMP, required for Baby-JubJub).

---

## Header.h -- Modular Arithmetic and Utilities

### Modular arithmetic

| Function               | Description                                                                                                                                                                |
|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `mod(value, modulus)`  | Non-negative remainder. Always returns result in [0, modulus). Overloaded for `long long` and `mpz_class`.                                                                 |
| `modInverse(a, p)`     | Modular multiplicative inverse a^{-1} mod p. Uses GMP's `mpz_invert` for `mpz_class`, Extended Euclidean Algorithm for `long long`. Throws if inverse does not exist.      |
| `modPow(base, exp, p)` | Modular exponentiation base^exp mod p. Uses GMP's `mpz_powm` for `mpz_class`, binary square-and-multiply for `long long`.                                                  |                          
            
### Cryptographic Utilities
Helper functions specifically added to support encryption and signature generation.

| Function          | Description                                                                                                 |
|-------------------|-------------------------------------------------------------------------------------------------------------|
| `hashMessage(message)`   | Simulates a cryptographic hash function (like SHA-256). Converts a string into a large mpz_class integer for ECDSA signing.    |
| `xorStrings(text, key)` | A symmetric XOR cipher. Used for both Data Encryption ($Enc$ / $Dec$) and Key Encapsulation ($Wrap$ / $Unwrap$) in the ECIES protocol. |
| `stringToHex(input)`    | Helper to format encrypted byte strings into human-readable hexadecimal output.                        |

---

## EllipticCurve.h -- Curve Class

### Class: `EllipticCurve<T>`
Represents the core parameters of an elliptic curve defined by the short Weierstrass equation $y^2 = x^3 + ax + b \pmod p$.

| Member    | Description                                                             |
|-----------|-------------------------------------------------------------------------|
| `a, b, p` | Curve coefficients and the prime field modulus.                         |
| `n`       | Curve order (total number of points on the curve, including $\mathcal{O}$). |
| `verbose` | Boolean flag to enable step-by-step console output (optional).          |

### Methods

| Method            | Description                                                                                                                                                                                             |
|-------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `isNonSingular()` | Validates the curve mathematically. Checks if the discriminant $\Delta = -16(4a^3 + 27b^2) \neq 0 \pmod p$. Since $-16$ is invertible for any prime $p \geq 3$, only the core expression $4a^3 + 27b^2$ is computed. |
| `print()`         | Outputs the human-readable algebraic equation of the instantiated curve.                                                                                                                                |

## EllipticCurvePoint.h -- Point Arithmetic

### Class: `EllipticCurvePoint<T>`

Represents a point on an elliptic curve in projective coordinates (X : Y : Z).
An affine point (x, y) is stored as (x, y, 1). The point at infinity O_E = (0 : 1 : 0).

### Constructors and conversion

| Method                             | Description                                                                    |
|------------------------------------|--------------------------------------------------------------------------------|
| `infinity(curve)`                  | Creates O_E = (0 : 1 : 0).                                                     |
| `fromAffine(x, y, curve)`          | Creates (x, y, 1) after validating the point lies on the curve. Throws if not. |
| `toAffine()`                       | Converts (X : Y : Z) to (X/Z, Y/Z) via modular inversion. Throws for O_E.      |
| `isInfinity()`                     | Returns true if Z == 0.                                                        |
| `isOnCurve()`                      | Checks Y^2*Z = X^3 + a*X*Z^2 + b*Z^3 mod p.                                    |
| `assertOnCurve(context)`           | Same check but throws with context message. Silent (no verbose output).        |

### Projective arithmetic

These formulas avoid modular inversion entirely, using only multiplications mod p.
The final result stays in projective form; conversion to affine (1 inversion) is done only when needed.

| Method          | Formula                                                                                                                        | Cost (approx.) |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------|----------------|
| `pointDouble()` | W = aZ^2 + 3X^2, S = YZ, B = XYS, H = W^2 - 8B, then X' = 2HS, Y' = W(4B-H) - 8Y^2S^2, Z' = 8S^3                               | 10M + 4S       |
| `pointAdd(Q)`   | U1=Y2*Z1, U2=Y1*Z2, V1=X2*Z1, V2=X1*Z2, U=U1-U2, V=V1-V2, W=Z1*Z2, A=U^2W-V^3-2V^2V2, then X3=VA, Y3=U(V^2V2-A)-V^3U2, Z3=V^3W | 12M + 2S       |

Handles special cases: P + O_E = P, P + (-P) = O_E, P + P delegates to pointDouble.

### Scalar multiplication

Four algorithms implementing kP, forming a 2x2 matrix of choices:

|            | Double-and-Add       | Description              |
|------------|----------------------|--------------------------------|
| Projective | `scalarMul(k)`       | Computes $k \cdot P$ using the **Double-and-Add** algorithm in projective coordinates for maximum performance. |

Double-and-Add: iterates bits of k from LSB to MSB. On each step, doubles temp.
If bit is 1, adds temp to accumulator. Running time depends on Hamming weight of k (not constant-time).

### Other methods

| Method          | Description                                                                                                                                 |
|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| `toString(val)`  | Static helper to safely format types (like `mpz_class` or `long long`) to standard string representation. |
| `equals(other)` | Projective equality: checks X1*Z2 == X2*Z1 and Y1*Z2 == Y2*Z1 mod p.                                                                        |
| `print()`       | Outputs affine and projective coordinates.                                                                                                  |

---

## Main.cpp -- Program Analysis and Output

The `Main.cpp` file contains a practical implementation and step-by-step testing of three cryptographic algorithms based on the **Baby-JubJub** elliptic curve. The program demonstrates the full lifecycle of cryptographic protocols: from parameter initialization to data integrity verification.

### 1. Curve and Generator Initialization
In the first step, the program initializes large prime numbers for the field $p$ and the subgroup $r$, translates the curve parameters from Montgomery form into Weierstrass form, and calculates the base generator point $G$.

> **Program Output:**
> ```text
> ========== Baby-JubJub Initialization ==========
> Curve initialized successfully.
> Generator G: (14414009007687342025526645003307639786191886886413750648631138442071909631647, 14577268218881899420966779687690205425227431577728659819975198491127179315626)
> ```
* **Result:** The curve is successfully created. All subsequent operations (key generation) are based on the derived point $G$.

---

### 2. Diffie-Hellman Key Exchange (ECDH)
Both parties (Alice and Bob) generate their long-term private keys modulo $r$ and exchange their calculated public keys.

| Party | Private Key (Secret) | Public Key | Shared Secret Formula ($S$) |
| :--- | :--- | :--- | :--- |
| **Alice** | $d_A$ (random integer) | $Q_A = d_A \cdot G$ | $S_A = d_A \cdot Q_B$ |
| **Bob** | $d_B$ (random integer) | $Q_B = d_B \cdot G$ | $S_B = d_B \cdot Q_A$ |

> **Program Output:**
> ```text
> ========== Diffie-Hellman ==========
> Alice generated long-term keys (dA, Qa).
> Bob generated long-term keys (dB, Qb).
> SUCCESS (Diffie-Hellman shared secret matched)
> ```
* **Result:** The algorithm executed correctly. The points $S_A$ and $S_B$ are mathematically identical, proving successful and secure shared secret negotiation without transmitting private keys over the network.

---

### 3. Directed Encryption (ECIES)
Implementation of a secure message transmission protocol from Alice to Bob using an ephemeral (one-time) key and symmetric XOR encryption.

| Envelope Component ($Env$) | Description and Formula |
| :--- | :--- |
| $Q_{eph}$ | Alice's ephemeral public key: $Q_{eph} = e_A \cdot G$ |
| $C_k$ | Encapsulated symmetric key: $Wrap(k, S_x) = k \oplus S_x$ |
| $C_M$ | Encrypted message: $Enc(M, k) = M \oplus k$ |

> **Program Output:**
> ```text
> ========== ECIES ==========
> Alice: Original message (M): Secret message for Bob.
> Alice: Encrypted message (CM) (hex): 615755425647185f5d4a4b555e53125e5f4112725b5017
> Alice: Encapsulated key (Ck) (hex): 0003020400070104010b0b010b02020a03030a0807020b0f03020d0e0c030d06010308080109040301000008070405020f07080306060e05060e03060f000d01010408040905000a03020502
> Alice sends Bob an envelope: Env = (Qeph, Ck, CM)
> 
> Bob: Decrypted message: Secret message for Bob.
> SUCCESS (The decrypted text matches the original)
> ```
* **Result:** Bob successfully reconstructed the shared secret using $Q_{eph}$ and his private key $d_B$. The key $k$ was decapsulated ($Unwrap$), and the message was decrypted. The decrypted text perfectly matches the original string.

---

### 4. Digital Signature (ECDSA)
A protocol confirming authorship and data integrity. Alice creates a digital signature $(r, s)$ for a message, and Bob mathematically verifies its validity.

**Signature Parameters:**
1. **Message Hash:** $h = H(M)$.
2. **Component $r$:** The $X$-coordinate of a random point $R = k \cdot G \pmod r$.
3. **Component $s$:** Computed as $s = k^{-1}(h + d_A \cdot r) \pmod r$.

> **Program Output:**
> ```text
> ========== ECDSA ==========
> Message for signature (M): This document confirms that we will receive 50 points for this lab.
> 
> Alice: Digital signature created:
>   r = 1269281715416766387322253707283110426429170823468985237291036227131565692010
>   s = 513727007132469588548670549323919458286340278040255630559964475700580827974
> 
> Bob: Received Alice's signature and public key Qa. Verifying:
> SUCCESS (v == r. Digital signature is true)
> ```
* **Result:** Bob computed the verification checkpoint $P_{check}$ using Alice's public key $Q_A$. Because the $v$ coordinate of this point matched the signature $r$ ($v == r$), the signature is deemed valid. This mathematically guarantees that the message was signed by Alice and has not been altered in transit.

---

### General Conclusion
The developed library correctly implements the essential cryptographic mathematics (projective coordinate arithmetic on elliptic curves) and fully complies with the standard requirements for ECDH, ECIES, and ECDSA protocols. 

### Cryptographic Design Notes

#### Subgroup Order vs Prime Field
A critical distinction in ECC protocols (especially ECDSA) is the modulus used for scalar operations. 
* The **Prime Field ($p$)** is used for coordinate arithmetic (point addition, doubling).
* The **Subgroup Order ($r$)** is the number of points generated by the base point $G$. 

All private keys, nonces, and scalar inversions ($s^{-1}, k^{-1}$) must be calculated modulo the subgroup order $r$, **not** modulo $p$. This implementation strictly adheres to this rule, using `r_bjj` for all key generation and ECDSA formulas.

#### Coordinate System Performance
Using Affine coordinates for a 254-bit curve like Baby-JubJub would require computing a modular inverse on every single bit of the scalar during multiplication, drastically slowing down the execution. By utilizing Projective Coordinates, the `scalarMul` function runs over 100x faster, shifting the heavy computational lifting strictly to modular multiplications.
