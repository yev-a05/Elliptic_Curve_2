#include "Header.h"
#include <iomanip> // Для виводу в шістнадцятковому форматі (hex)

// Проста функція XOR, яка накладає ключ на текст (для Enc і для Wrap)
std::string xorStrings(const std::string& text, const std::string& key) {
    std::string result = text;
    for (size_t i = 0; i < text.size(); ++i) {
        result[i] = text[i] ^ key[i % key.size()];
    }
    return result;
}

// Допоміжна функція, щоб виводити зашифровані (нечитабельні) байти у красивому вигляді
std::string stringToHex(const std::string& input) {
    std::ostringstream oss;
    for (unsigned char c : input) {
        oss << std::hex << std::setw(2) << std::setfill('0') << (int)c;
    }
    return oss.str();
}

// Проста геш-функція (імітація SHA-256)
// Перетворює текстове повідомлення у велике число mpz_class
mpz_class hashMessage(const std::string& message) {
    std::string hex_hash = "";
    // Беремо кожен символ і записуємо його як 16-річне число
    for (char c : message) {
        std::stringstream ss;
        ss << std::hex << std::setw(2) << std::setfill('0') << (int)(unsigned char)c;
        hex_hash += ss.str();
    }
    mpz_class h;
    h.set_str(hex_hash, 16); // Створюємо велике число з hex-рядка
    return h;
}

int main() {
    std::cout << "========== Baby-JubJub Initialization ==========\n\n";

    // 1. Prime field
    mpz_class p_bjj("21888242871839275222246405745257275088548364400416034343698204186575808495617");

    // 2. Subgroup order r (this is the main modulus for keys)
    mpz_class r_bjj("2736030358979909402780800718157159386076813972158567259200215660948447373041");
    mpz_class n_bjj = 8 * r_bjj;

    // 3. Weierstrass parameters conversion (from Montgomery A=168698)
    mpz_class A_mont(168698);
    mpz_class inv3 = modInverse(mpz_class(3), p_bjj);
    mpz_class inv27 = modInverse(mpz_class(27), p_bjj);
    mpz_class A2 = mod(A_mont * A_mont, p_bjj);
    mpz_class A3 = mod(A2 * A_mont, p_bjj);
    mpz_class a_bjj = mod((3 - A2) * inv3, p_bjj);
    mpz_class b_bjj = mod((2 * A3 - 9 * A_mont) * inv27, p_bjj);

    // Initialize the Curve
    EllipticCurve<mpz_class> bjj(a_bjj, b_bjj, p_bjj, n_bjj, false);

    // 4. Generator Point Initialization
    mpz_class ex("5299619240641551281634865583518297030282874472190772894086521144482721001553");
    mpz_class ey("16950150798460657717958625567821834550301663161624707787222815936182638968203");
    mpz_class u = mod((1 + ey) * modInverse(mod(1 - ey, p_bjj), p_bjj), p_bjj);
    mpz_class v = mod(u * modInverse(ex, p_bjj), p_bjj);
    mpz_class gx = mod(u + A_mont * inv3, p_bjj);
    mpz_class gy = v;

    auto G = EllipticCurvePoint<mpz_class>::fromAffine(gx, gy, &bjj);

    std::cout << "Curve initialized successfully." << std::endl;
    std::cout << "Generator G: ";
    G.print();
    std::cout << std::endl;

    // ========== Setup random number generator (For generating Private Keys) ==========
    gmp_randclass rng(gmp_randinit_mt);
    std::random_device rd;
    mpz_class seed = mpz_class(static_cast<unsigned long>(rd())) * mpz_class(0x10000UL) * mpz_class(0x10000UL) + mpz_class(static_cast<unsigned long>(rd()));
    rng.seed(seed);

    // ========== Task 1: DIFFIE-HELLMAN ==========
    std::cout << "========== Diffie-Hellman ==========\n";

    // 1. Аліса генерує свій довгостроковий закритий ключ dA та відкритий Qa
    mpz_class dA = rng.get_z_range(r_bjj);
    if (dA == 0) dA = 1;
    auto Qa = G.scalarMul(dA);
    std::cout << "Alice generated long-term keys (dA, Qa).\n";

    // 2. Боб генерує свій довгостроковий закритий ключ dB та відкритий Qb
    mpz_class dB = rng.get_z_range(r_bjj);
    if (dB == 0) dB = 1;
    auto Qb = G.scalarMul(dB);
    std::cout << "Bob generated long-term keys (dB, Qb).\n";

    // 3. Обмін та узгодження спільного секрету
    auto S_Alice_DH = Qb.scalarMul(dA); // Аліса множить публічний ключ Боба на свій секрет
    auto S_Bob_DH = Qa.scalarMul(dB);   // Боб множить публічний ключ Аліси на свій секрет

    if (S_Alice_DH.equals(S_Bob_DH)) {
        std::cout << "SUCCESS (Diffie-Hellman shared secret matched) \n\n";
    }
    else {
        std::cout << "ERROR (The secrets don't match)\n\n";
    }

    //  ==========  Task 2: Directional encryption (ECIES) ========== 
    std::cout << "========== ECIES ==========\n";

    std::string M = "Secret message for Bob.";
    std::cout << "Alice: Original message (M): " << M << "\n";

    // 1. Аліса генерує випадковий симетричний ключ k
    mpz_class k_rand = rng.get_z_range(r_bjj);
    std::string k = k_rand.get_str();

    // 2. Шифруємо повідомлення симетричним ключем: CM = Enc(M, k)
    std::string CM = xorStrings(M, k);
    std::cout << "Alice: Encrypted message (CM) (hex): " << stringToHex(CM) << "\n";

    // 3. Аліса генерує ефемерну (одноразову) ключову пару (eA, Qeph)
    mpz_class eA = rng.get_z_range(r_bjj);
    if (eA == 0) eA = 1;
    auto Qeph = G.scalarMul(eA);

    // 4. Аліса обчислює спільний секрет з відкритим ключем Боба: S = eA * Qb
    auto S_Alice_Enc = Qb.scalarMul(eA);
    mpz_class Sx_Alice = S_Alice_Enc.toAffine().first;
    std::string Sx_str = Sx_Alice.get_str();

    // 5. Аліса інкапсулює симетричний ключ: Ck = Wrap(k, Sx)
    std::string Ck = xorStrings(k, Sx_str);
    std::cout << "Alice: Encapsulated key (Ck) (hex): " << stringToHex(Ck) << "\n";

    std::cout << "Alice sends Bob an envelope: Env = (Qeph, Ck, CM)\n\n";

    // --- Розшифрування (Сторона Боба) ---
    // 1. Боб обчислює спільний секрет: S = dB * Qeph
    auto S_Bob_Dec = Qeph.scalarMul(dB);
    mpz_class Sx_Bob = S_Bob_Dec.toAffine().first;

    // 2. Боб деінкапсулює ключ: k = Unwrap(Ck, Sx)
    std::string k_decrypted = xorStrings(Ck, Sx_Bob.get_str());

    // 3. Боб розшифровує повідомлення: M = Dec(CM, k)
    std::string M_decrypted = xorStrings(CM, k_decrypted);

    std::cout << "Bob: Decrypted message: " << M_decrypted << "\n";
    if (M == M_decrypted) {
        std::cout << "SUCCESS (The decrypted text matches the original)\n\n";
    }

    // ==========  Task 3: Digital Signature ECDSA ==========
    std::cout << "========== ECDSA ==========\n";

    std::string M_sign = "This document confirms that we will receive 50 points for this lab.";
    std::cout << "Message for signature (M): " << M_sign << "\n";

    // --- Підпис (Сторона Аліси) ---
    // 1. Геш повідомлення: h = H(M)
    mpz_class h = hashMessage(M_sign);

    mpz_class r_sig, s_sig;
    do {
        // 2. Випадкове число k
        mpz_class k_sig = rng.get_z_range(r_bjj);
        if (k_sig == 0) continue;

        // kP = (x1, y1)
        auto R_point = G.scalarMul(k_sig);

        // r = x1 mod n (тут n = r_bjj)
        r_sig = mod(R_point.toAffine().first, r_bjj);
        if (r_sig == 0) continue;

        // s = k^-1 * (h + dA * r) mod n
        mpz_class k_inv = modInverse(k_sig, r_bjj);
        mpz_class dA_times_r = mod(dA * r_sig, r_bjj);
        mpz_class h_plus_dr = mod(h + dA_times_r, r_bjj);

        s_sig = mod(k_inv * h_plus_dr, r_bjj);

    } while (s_sig == 0);

    std::cout << "\nAlice: Digital signature created:\n";
    std::cout << "  r = " << r_sig << "\n";
    std::cout << "  s = " << s_sig << "\n\n";

    // Checking (Bob's side)
    std::cout << "Bob: Received Alice's signature and public key Qa. Verifying:\n";

    // Геш повідомлення: h = H(M)
    mpz_class h_check = hashMessage(M_sign);

    // u1 = s^-1 * h mod n  ||  u2 = s^-1 * r mod n
    mpz_class w = modInverse(s_sig, r_bjj);
    mpz_class u1 = mod(w * h_check, r_bjj);
    mpz_class u2 = mod(w * r_sig, r_bjj);

    // P_check = u1*G + u2*Qa
    auto P1 = G.scalarMul(u1);
    auto P2 = Qa.scalarMul(u2);
    auto P_check = P1.pointAdd(P2);

    if (P_check.isInfinity()) {
        std::cout << "ERROR (Checkpoint equals infinity)\n";
    }
    else {
        // v = x0 mod n
        mpz_class v = mod(P_check.toAffine().first, r_bjj);

        if (v == r_sig) {
            std::cout << "SUCCESS (v == r. Digital signature is true)\n";
        }
        else {
            std::cout << "ERROR (v != r. Digital signature forged)\n";
        }
    }

    return 0;
}