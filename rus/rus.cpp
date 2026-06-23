#include "rus.h"
#include <mpfr.h>

RUS::RUS(int debug_level, const mpf_class &epsilon, int effort_level, int pslq_iters, CRITERION criterion,
         bool pslq_iters_from_cli)
    : idb(debug_level),
      eps(epsilon), effort(effort_level), pslq_max_iter(pslq_iters), pslq_iters_from_cli(pslq_iters_from_cli),
      crit(criterion)
{
    mpf_set_default_prec(2048);
    mpfr_set_default_prec(2048);
}

void RUS::run(const mpf_class &theta, circuit &best_cir)
{
    if (idb >= 1)
        std::cout << "Stage 1: PSLQ\n";

    if (!pslq_iters_from_cli)
        pslq_max_iter = DEFAULT_PSLQ_MAX_ITER;

    std::vector<std::complex<mpf_class>> x = create_vector(theta);
    auto get_error = gen_error_function(theta);
    PslqComplex pslq(x, eps, pslq_max_iter, idb, effort, get_error);
    pslq.run();
    auto results = pslq.get_results();

    int r_count = 0;
    for (auto r : results)
    {
        ++r_count;
        Normalization::OmegaRing z = create_omega_ring(r);

        auto cz = z.toComplex(0);

        if (idb >= 2)
        {
            print_OmegaRing(z);

            // todo: use mpf_class instead of double (modify the .toComplex() function in OmegaRing)
            std::cout << "real epsilon = " << 2 * ((decltype(cz)(cos(theta / 2).get_d(), sin(theta / 2).get_d())) * cz).imag() / abs(cz) << std::endl;
        }
        if (idb >= 1)
            std::cout << "Stage 2-" << r_count << ": Normalization\n";

        Normalization normSolver(z);
        normSolver.solve();
        matrix2x2<mpz_class> result2 = normSolver.get_result();
        if (idb >= 1)
            std::cout << "Stage 3-" << r_count << ": Decomposition\n";

        circuit cir = exactDecomposer::decompose(result2);
        switch (crit)
        {
        case G_COUNT:
            if (r_count == 1)
                best_cir = cir;
            else if (count_gate(cir) < count_gate(best_cir))
                best_cir = cir;
            break;
        case T_COUNT:
        default:
            if (r_count == 1)
                best_cir = cir;
            else if (count_t(cir) < count_t(best_cir))
                best_cir = cir;
            break;
        }
    }
}

mpf_class RUS::parse_theta(const std::string &theta_str)
{
    if (theta_str.empty())
        return mpf_class(0);
    std::size_t found = theta_str.find("pi");
    if (found != std::string::npos)
    {
        std::string str_without_pi = theta_str.substr(0, found);
        if (found < theta_str.size() - 2)
            str_without_pi += theta_str.substr(found + 2);
        std::size_t found_div = str_without_pi.find("/");
        if (found_div != std::string::npos)
        {
            auto num_str = str_without_pi.substr(0, found_div);
            auto den_str = str_without_pi.substr(found_div + 1);
            mpf_class numerator = num_str.size() == 0 ? mpf_class(1) : mpf_class(num_str);
            mpf_class denominator = den_str.size() == 0 ? mpf_class(1) : mpf_class(den_str);
            return (numerator / denominator) * M_PI;
        }
        else if (str_without_pi.empty())
            return mpf_class(1) * M_PI;
        else
            return mpf_class(str_without_pi) * M_PI;
    }
    return mpf_class(theta_str);
}

std::vector<std::complex<mpf_class>> RUS::create_vector(const double &theta)
{
    using std::sin;
    std::vector<std::complex<mpf_class>> v;
    for (size_t i = 0; i < 4; i++)
    {
        v.push_back(mpf_class(sin(theta / 2 + i * M_PI / 4)));
    }
    return v;
}

std::vector<std::complex<mpf_class>> RUS::create_vector(const mpf_class &theta)
{
    std::vector<std::complex<mpf_class>> v;
    for (size_t i = 0; i < 4; i++)
    {
        mpf_class x = theta / 2 + i * M_PI / 4;
        mpf_class y;

        y = sin(x);
        v.push_back(y);
    }
    return v;
}

Normalization::OmegaRing RUS::create_omega_ring(const std::vector<std::complex<mpf_class>> &v)
{
    std::vector<mpz_class> vz;
    for (auto &z : v)
    {
        vz.push_back(mpz_class(z.real()));
    }
    return Normalization::OmegaRing(vz[0], vz[1], vz[2], vz[3]);
}

std::function<real_t(const PslqComplex::Complex &min_y_val, const PslqComplex::ComplexVector &b_col)> RUS::gen_error_function(const mpf_class &theta)
{
    // todo: use mpf_class instead of double (modify the .toComplex() function in OmegaRing)
    return [theta](const std::complex<mpf_class> &min_y_val, const std::vector<std::complex<mpf_class>> &b_col) -> mpf_class
    {
        Normalization::OmegaRing z = create_omega_ring(b_col);
        auto cz = z.toComplex(0);
        return abs(2 * ((decltype(cz)(cos(theta / 2).get_d(), sin(theta / 2).get_d())) * cz).imag() / abs(cz));
    };
}

mpf_class RUS::sin(const mpf_class &theta)
{
    mpf_class x = theta;
    mpf_class y;

    mpfr_t xx;
    mpfr_t yy;
    mpfr_init(xx);
    mpfr_init(yy);

    mpfr_set_f(xx, x.get_mpf_t(), MPFR_RNDN);
    mpfr_sin(yy, xx, MPFR_RNDN);

    mpfr_get_f(y.get_mpf_t(), yy, MPFR_RNDN);

    mpfr_clear(xx);
    mpfr_clear(yy);

    return y;
}

mpf_class RUS::cos(const mpf_class &theta)
{
    mpf_class x = theta;
    mpf_class y;

    mpfr_t xx;
    mpfr_t yy;
    mpfr_init(xx);
    mpfr_init(yy);

    mpfr_set_f(xx, x.get_mpf_t(), MPFR_RNDN);
    mpfr_cos(yy, xx, MPFR_RNDN);

    mpfr_get_f(y.get_mpf_t(), yy, MPFR_RNDN);

    mpfr_clear(xx);
    mpfr_clear(yy);

    return y;
}

int RUS::count_t(const circuit &cir)
{
    std::ostringstream oss;
    cir.toStream(oss);
    const std::string str = oss.str();
    return std::count(str.cbegin(), str.cend(), 'T');
}

int RUS::count_gate(const circuit &cir)
{
    std::ostringstream oss;
    cir.toStream(oss);
    const std::string str = oss.str();
    return std::count(str.cbegin(), str.cend(), '\n');
}

void QasmGenerator::to_qasm(circuit &cir, std::ostream &os, std::string qbit_def, bool without_header)
{
    if (!without_header)
    {
        os << "OPENQASM 2.0;\n";
        os << "include \"qelib1.inc\";\n\n";
    }
    if (qbit_def.empty())
    {
        os << "qreg q[1];\n";
        qbit_def = "q[0]";
    }
    if (cir.size() == 0 || (cir.size() == 1 && cir[0] == gateLibrary::Id))
    {
        os << "id " << qbit_def << ";\n";
        return;
    }
    auto name_qasm = gen_name_qasm();
    fill_in_gatename(name_qasm, qbit_def);

    for (size_t i = 0; i < cir.size(); i++)
    {
        if (cir[i] == gateLibrary::Id)
            continue;
        if (cir[i] >= gateLibrary::GLw1 && cir[i] <= gateLibrary::GLw7)
            throw std::runtime_error("Unsupported gate in QASM: " + std::to_string(cir[i]));
        os << name_qasm[cir[i]];
    }
}

void QasmGenerator::to_rus_qasm(circuit &cir, std::ostream &os, std::string qbit_def, std::string ancil_def, bool without_header, bool with_tail_syntax)
{
    if (!without_header)
    {
        os << "OPENQASM 2.0;\n";
        os << "include \"qelib1.inc\";\n\n";
    }
    if (qbit_def.empty())
    {
        os << "qreg q[1];\n";
        qbit_def = "q[0]";
    }
    if (ancil_def.empty())
    {
        os << "qreg a[1];\n";
        ancil_def = "a[0]";
    }
    if (cir.size() == 0 || (cir.size() == 1 && cir[0] == gateLibrary::Id))
    {
        os << "id " << qbit_def << ";\n";
        return;
    }

    os << "cx " << qbit_def << "," << ancil_def << ";\n";
    to_qasm(cir, os, ancil_def, true);
    os << "cx " << qbit_def << "," << ancil_def << ";\n";
    if (with_tail_syntax)
        os << "rus " << ancil_def << "==0;";
    else
        os << "// only valid when measure " << ancil_def << " == 0";
}

std::vector<std::string> QasmGenerator::gen_name_qasm()
{
    std::vector<std::string> name_qasm(gateLibrary::GLw7 + 10);

    name_qasm[gateLibrary::Id] = "id __gate_name__;";
    name_qasm[gateLibrary::T] = "t __gate_name__;";
    name_qasm[gateLibrary::P] = "s __gate_name__;";
    name_qasm[gateLibrary::TP] = "tdg __gate_name__;z __gate_name__;";
    name_qasm[gateLibrary::Z] = "z __gate_name__;";
    name_qasm[gateLibrary::TZ] = "t __gate_name__;z __gate_name__;";
    name_qasm[gateLibrary::Pd] = "sdg __gate_name__;";
    name_qasm[gateLibrary::Td] = "tdg __gate_name__;";
    name_qasm[gateLibrary::H] = "h __gate_name__;";
    name_qasm[gateLibrary::X] = "x __gate_name__;";
    name_qasm[gateLibrary::Y] = "y __gate_name__;";
    for (int i = 1; i < 8; i++)
    {
        std::stringstream ss;
        ss << "GLw" << i << " __gate_name__"
           << ";";
        name_qasm[(gateLibrary::GLw1 - 1 + i)] = ss.str();
    }
    return name_qasm;
}

void QasmGenerator::fill_in_gatename(std::vector<std::string> &name_qasm, const std::string &gate_name)
{
    for (size_t i = 0; i < name_qasm.size(); i++)
    {
        std::string &str = name_qasm[i];

        size_t pos;
        while (pos = str.find("__gate_name__"), pos != std::string::npos)
            str.replace(pos, 13, gate_name);

        str = split_by_semicolon(str);
    }
}

std::string QasmGenerator::split_by_semicolon(const std::string &str)
{
    std::string res;
    std::string s = str;
    size_t pos = 0;
    while ((pos = s.find(";")) != std::string::npos)
    {
        std::string token = s.substr(0, pos);
        res = res + token + ";\n";
        s.erase(0, pos + 1);
    }
    return res;
}