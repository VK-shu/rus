#include "rus.h"
#include <mpfr.h>
#include <cmath>

namespace
{
mpf_class mpf_abs(const mpf_class &x)
{
    return x >= 0 ? x : -x;
}

mpf_class mpf_hypot(const mpf_class &re, const mpf_class &im)
{
    mpf_class sum = re * re + im * im;
    mpf_class r;
    mpf_sqrt(r.get_mpf_t(), sum.get_mpf_t());
    return r;
}
} // namespace

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

    std::vector<std::complex<mpf_class>> x = create_vector(theta);
    auto get_error = gen_error_function(theta);

    unsigned pslq_limit = static_cast<unsigned>(pslq_max_iter);
    const unsigned pslq_cap = pslq_iters_from_cli ? pslq_limit : static_cast<unsigned>(DEFAULT_PSLQ_MAX_ITER);
    if (!pslq_iters_from_cli)
        pslq_limit = static_cast<unsigned>(INITIAL_PSLQ_MAX_ITER);

    std::vector<PslqComplex::ComplexVector> results;
    while (true)
    {
        PslqComplex pslq(x, eps, pslq_limit, idb, effort, get_error);
        pslq.run();
        results = pslq.get_results();

        if (!results.empty() || pslq_iters_from_cli)
            break;
        if (pslq_limit >= pslq_cap)
            break;

        const unsigned next_limit = std::min(pslq_limit * 2u, pslq_cap);
        if (idb >= 1)
            std::cerr << "PSLQ: no candidates at limit " << pslq_limit << ", retrying with " << next_limit
                      << std::endl;
        pslq_limit = next_limit;
    }

    if (results.empty())
        throw std::runtime_error("RUS: PSLQ produced no candidates (try increasing -F, -I, or relaxing -E)");

    int r_count = 0;
    int success_count = 0;
    for (auto r : results)
    {
        ++r_count;
        Normalization::OmegaRing z = create_omega_ring(r);

        if (idb >= 2)
        {
            print_OmegaRing(z);
            std::cout << "real epsilon = " << compute_phase_error(z, theta) << std::endl;
        }
        if (idb >= 1)
            std::cout << "Stage 2-" << r_count << ": Normalization\n";

        matrix2x2<mpz_class> result2;
        try
        {
            Normalization normSolver(z);
            normSolver.solve();
            result2 = normSolver.get_result();
        }
        catch (const std::exception &e)
        {
            if (idb >= 1)
                std::cerr << "RUS: normalization failed for candidate " << r_count << ": " << e.what() << std::endl;
            continue;
        }

        if (idb >= 1)
            std::cout << "Stage 3-" << r_count << ": Decomposition\n";

        circuit cir;
        if (!exactDecomposer::instance().decompose(result2, cir))
        {
            if (idb >= 1)
                std::cerr << "RUS: decomposition failed for candidate " << r_count << std::endl;
            continue;
        }

        switch (crit)
        {
        case G_COUNT:
            if (success_count == 0)
                best_cir = cir;
            else if (count_gate(cir) < count_gate(best_cir))
                best_cir = cir;
            break;
        case T_COUNT:
        default:
            if (success_count == 0)
                best_cir = cir;
            else if (count_t(cir) < count_t(best_cir))
                best_cir = cir;
            break;
        }
        ++success_count;
    }

    if (success_count == 0)
        throw std::runtime_error("RUS: all candidates failed normalization or decomposition");

    if (best_cir.empty())
        throw std::runtime_error("RUS: no circuit selected after synthesis");
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

std::complex<mpf_class> RUS::omega_ring_to_mpf_complex(const Normalization::OmegaRing &z, int de)
{
    const mpf_class sqrt2ov2 = SQRT2 / 2;
    mpf_class denom(1);
    if (de > 0)
    {
        int denom_power2 = de / 2;
        denom = 1;
        for (int i = 0; i < denom_power2; ++i)
            denom *= 2;
    }

    mpf_class arr[4];
    for (int i = 0; i < 4; ++i)
        arr[i] = mpf_class(z[i]) / denom;

    if (de % 2 == 1)
        for (int i = 0; i < 4; ++i)
            arr[i] *= sqrt2ov2;

    const mpf_class real = arr[0] + sqrt2ov2 * (arr[1] - arr[3]);
    const mpf_class imag = arr[2] + sqrt2ov2 * (arr[1] + arr[3]);
    return std::complex<mpf_class>(real, imag);
}

mpf_class RUS::compute_phase_error(const Normalization::OmegaRing &z, const mpf_class &theta)
{
    const std::complex<mpf_class> cz = omega_ring_to_mpf_complex(z, 0);
    const mpf_class cz_abs = mpf_hypot(cz.real(), cz.imag());
    if (cz_abs == 0)
        return mpf_class(0);

    const std::complex<mpf_class> target(cos(theta / 2), sin(theta / 2));
    const std::complex<mpf_class> prod = target * cz;
    return mpf_abs(2 * prod.imag() / cz_abs);
}

std::function<real_t(const PslqComplex::Complex &min_y_val, const PslqComplex::ComplexVector &b_col)> RUS::gen_error_function(const mpf_class &theta)
{
    return [theta](const std::complex<mpf_class> &, const std::vector<std::complex<mpf_class>> &b_col) -> mpf_class
    {
        Normalization::OmegaRing z = create_omega_ring(b_col);
        return compute_phase_error(z, theta);
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
    int n = 0;
    for (int g : cir)
    {
        switch (g)
        {
        case gateLibrary::T:
        case gateLibrary::Td:
        case gateLibrary::TP:
        case gateLibrary::TZ:
            ++n;
            break;
        default:
            break;
        }
    }
    return n;
}

int RUS::count_gate(const circuit &cir)
{
    if (cir.empty())
        return 1;
    return static_cast<int>(cir.size());
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