/**
 * @file
 * @copyright Copyright 2021. Tom de Geus. All rights reserved.
 * @par license This project is released under the MIT License.
 */

#ifndef APTS_H
#define APTS_H

#include <cmath>
#include <prrng.h>
#include <xtensor/xmath.hpp>
#include <xtensor/xtensor.hpp>

/**
 * \cond
 */
#define APTS_QUOTEHELPER(x) #x
#define APTS_QUOTE(x) APTS_QUOTEHELPER(x)

#define APTS_ASSERT_IMPL(expr, assertion, file, line, function) \
    if (!(expr)) { \
        throw assertion( \
            std::string(file) + ":" + std::to_string(line) + " (" + std::string(function) + ")" + \
            ": assertion failed (" #expr ") \n\t"); \
    }

/**
 * \endcond
 */

/**
 * Library version.
 *
 * Either:
 *
 * -   Configure using CMake at install time. Internally uses:
 *
 *         python -c "from setuptools_scm import get_version; print(get_version())"
 *
 * -   Define externally using:
 *
 *         -DAPTS_VERSION="`python -c "from setuptools_scm import get_version;
 * print(get_version())"`"
 *
 *     From the root of this project. This is what ``setup.py`` does.
 *
 * Note that both ``CMakeLists.txt`` and ``setup.py`` will construct the version using
 * ``setuptools_scm``. Tip: use the environment variable ``SETUPTOOLS_SCM_PRETEND_VERSION`` to
 * overwrite the automatic version.
 */
#ifndef APTS_VERSION
#define APTS_VERSION "@PROJECT_VERSION@"
#endif

/**
 * All assertions are implementation as:
 *
 *     APTS_ASSERT(...)
 *
 * They can be enabled by:
 *
 *     #define APTS_ENABLE_ASSERT
 *
 * (before including prrng).
 * The advantage is that:
 *
 * -   File and line-number are displayed if the assertion fails.
 * -   prrng's assertions can be enabled/disabled independently from those of other libraries.
 *
 * \throw std::runtime_error
 */
#ifdef APTS_ENABLE_ASSERT
#define APTS_ASSERT(expr, assertion) \
    APTS_ASSERT_IMPL(expr, assertion, __FILE__, __LINE__, __FUNCTION__)
#else
#define APTS_ASSERT(expr, assertion)
#endif

/**
 * Assertion that cannot be switched off. Implement assertion by::
 *
 *     APTS_REQUIRE(...)
 *
 * \throw std::runtime_error
 */
#define APTS_REQUIRE(expr, assertion) \
    APTS_ASSERT_IMPL(expr, assertion, __FILE__, __LINE__, __FUNCTION__)

/**
 * @brief Analytical Prandtl-Tomlison Solutions
 */
namespace apts {

namespace detail {

/**
 * Remove " from string.
 *
 * @param arg Input string.
 * @return String without "
 */
inline std::string unquote(const std::string& arg)
{
    std::string ret = arg;
    ret.erase(std::remove(ret.begin(), ret.end(), '\"'), ret.end());
    return ret;
}

} // namespace detail

/**
 * Return version string, for example `"0.1.0"`
 *
 * @return std::string
 */
inline std::string version()
{
    return detail::unquote(std::string(APTS_QUOTE(APTS_VERSION)));
}

/**
 * @brief Quadratic potential.
 */
class Quadratic {
protected:
    double m_w;
    double m_m;
    double m_eta;
    double m_kappa;
    double m_F;

    double m_r0;
    double m_v0;
    double m_r0prime;
    double m_v0prime;
    double m_rmax;
    double m_delta_r;

    double m_lambda;
    double m_omega;
    double m_chi;
    double m_phi;
    double m_L;
    double m_Q;
    double m_phase;
    double m_pi;

    double m_vc;
    bool m_vc_computed;

    double m_tauexit;
    bool m_tauexit_computed;

public:
    /**
     * @param v0 Initial velocity.
     * @param w Width of the potential.
     * @param m Mass.
     * @param eta Damping.
     * @param mu Stiffness.
     * @param f External force.
     */
    Quadratic(
        double v0 = 0,
        double w = 1,
        double m = 1,
        double eta = 0.1,
        double mu = 1,
        double f = 0)
        : m_m(m), m_eta(eta), m_kappa(mu), m_F(f)
    {
        APTS_REQUIRE(std::pow(eta, 2.0) <= 4.0 * m, std::out_of_range);

        m_lambda = m_eta / (2.0 * m_m);
        m_omega = std::sqrt(m_kappa / m_m - std::pow(m_lambda, 2.0));
        m_phase = std::atan(m_omega / m_lambda);
        m_pi = xt::numeric_constants<double>::PI;
        m_delta_r = m_F / m_kappa;

        this->set_w(w);
        this->set_v0(v0);
    }

    /**
     * @brief Change the width of the potential.
     * @param value Width.
     */
    void set_w(double value)
    {
        m_w = value;
        m_r0 = -0.5 * m_w;
        m_rmax = 0.5 * m_w;
        m_r0prime = m_r0 - m_delta_r;
        m_vc_computed = false;
        m_tauexit_computed = false;
    }

    /**
     * @brief Set the initial velocity.
     * @note Calling this allows to avoid the constructor and gives a minor speed-up.
     * @param value Velocity.
     */
    void set_v0(double value)
    {
        m_v0 = value;
        m_v0prime = value;

        double re_alpha = 0.5 * m_r0prime;
        double im_alpha = 0.5 * (m_lambda * m_r0prime + m_v0prime) / m_omega;

        if (re_alpha >= 0.0) {
            m_chi = 0.0;
        }
        else {
            if (im_alpha >= 0.0) {
                m_chi = 1.0;
            }
            else {
                m_chi = -1.0;
            }
        }

        m_phi = m_chi * m_pi - std::atan(m_lambda / m_omega + m_v0prime / (m_omega * m_r0prime));
        m_L = std::sqrt(
            (std::pow(m_omega * m_r0prime, 2.0) + std::pow(m_lambda * m_r0prime + m_v0prime, 2.0)) /
            std::pow(m_omega, 2.0));
        m_Q = m_L * m_lambda * std::sqrt(1 + std::pow(m_omega / m_lambda, 2.0));

        m_vc_computed = false;
        m_tauexit_computed = false;
    }

    /**
     * @brief Parameter: width
     * @return double
     */
    double w() const
    {
        return m_w;
    }

    /**
     * @brief Parameter: mass
     * @return double
     */
    double m() const
    {
        return m_m;
    }

    /**
     * @brief Parameter: damping coefficient
     * @return double
     */
    double eta() const
    {
        return m_eta;
    }

    /**
     * @brief Parameter: elastic modulus
     * @return double
     */
    double kappa() const
    {
        return m_kappa;
    }

    /**
     * @brief Parameter: external force
     * @return double
     */
    double F() const
    {
        return m_F;
    }

    /**
     * @brief Parameter: initial position
     * @return double
     */
    double r0() const
    {
        return m_r0;
    }

    /**
     * @brief Parameter: initial velocity
     * @return double
     */
    double v0() const
    {
        return m_v0;
    }

    /**
     * @brief Parameter: delta r
     * @return double
     */
    double delta_r() const
    {
        return m_delta_r;
    }

    /**
     * @brief Parameter: delta v
     * @return double
     */
    double delta_v() const
    {
        return m_v0 - m_v0prime;
    }

    /**
     * @brief Parameter: lambda
     * @return double
     */
    double lambda() const
    {
        return m_lambda;
    }

    /**
     * @brief Parameter: omega
     * @return double
     */
    double omega() const
    {
        return m_omega;
    }

    /**
     * @brief Parameter: chi
     * @return double
     */
    double chi() const
    {
        return m_chi;
    }

    /**
     * @brief Parameter: phi
     * @return double
     */
    double phi() const
    {
        return m_phi;
    }

    /**
     * @brief Parameter: L
     * @return double
     */
    double L() const
    {
        return m_L;
    }

    /**
     * @brief Parameter: Q
     * @return double
     */
    double Q() const
    {
        return m_Q;
    }

    /**
     * @brief Position at certain time.
     * @param tau Time.
     * @return double
     */
    double r_scalar(double tau) const
    {
        return m_L * std::exp(-m_lambda * tau) * std::cos(m_omega * tau + m_phi) + m_delta_r;
    }

    /**
     * @brief Velocity at certain time.
     * @param tau Time.
     * @return double
     */
    double v_scalar(double tau) const
    {
        return -m_Q * std::exp(-m_lambda * tau) * std::cos(m_omega * tau + m_phi - m_phase);
    }

    /**
     * @brief Position at different times.
     * @param tau Array of time.
     * @return Array of position.
     */
    template <class T>
    T r(const T& tau) const
    {
        return m_L * xt::exp(-m_lambda * tau) * xt::cos(m_omega * tau + m_phi) + m_delta_r;
    }

    /**
     * @brief Velocity at different times.
     * @param tau Array of time.
     * @return Array of velocity.
     */
    template <class T>
    T v(const T& tau) const
    {
        return -m_Q * xt::exp(-m_lambda * tau) * xt::cos(m_omega * tau + m_phi - m_phase);
    }

    /**
     * @brief Time at the maximum position take could be reached if the potential would be infinite.
     * @return double
     */
    double tau_at_rmax() const
    {
        double tp = (1.5 * m_pi + m_phase - m_phi) / m_omega;
        double tm = (-0.5 * m_pi + m_phase - m_phi) / m_omega;

        if (tp > 0 && tm > 0) {
            return std::min(tp, tm);
        }

        if (tp > 0) {
            return tp;
        }

        return tm;
    }

    /**
     * @brief Return `true` if the particle exits the well.
     * @return bool
     */
    bool exits() const
    {
        return this->r_scalar(this->tau_at_rmax()) > m_rmax;
    }

    /**
     * @brief Time when the particle exists the well.
     * This requires a minimisation, that is here performed with a simple Newton-Raphson method.
     * @return double
     */
    double tau_exit()
    {
        if (m_tauexit_computed) {
            return m_tauexit;
        }

        APTS_REQUIRE(this->exits(), std::out_of_range);

        double tau = 0.5 * this->tau_at_rmax();

        for (size_t i = 0; i < 100; ++i) {

            double x = this->r_scalar(tau);
            double v = this->v_scalar(tau);
            double f = x - m_rmax;

            if (std::abs(f / m_rmax) < 1e-6) {
                m_tauexit = tau;
                m_tauexit_computed = true;
                return m_tauexit;
            }

            tau -= f / v;
        }

        throw std::runtime_error("tau_exit: failure to converge");
    }

    /**
     * @brief Initial velocity below which the particle does not exit the well.
     * @return double
     */
    double vc()
    {
        if (m_vc_computed) {
            return m_vc;
        }

        double v0_bak = m_v0;

        double v0 = 0.0;
        double v1 = 10.0;

        // initial search
        for (size_t i = 0; i < 101; ++i) {
            this->set_v0(v1);

            if (this->exits()) {
                break;
            }

            v1 *= 2.0;

            if (i >= 100) {
                throw std::runtime_error("vc: failure to find an initial guess");
            }
        }

        for (size_t i = 0; i < 1000; ++i) {
            double vm = 0.5 * (v0 + v1);
            this->set_v0(vm);

            if (this->exits()) {
                v1 = vm;
            }
            else {
                v0 = vm;
            }

            if ((1.0 - v0 / v1) < 1e-6) {
                this->set_v0(v0_bak);
                m_vc = 0.5 * (v0 + v1);
                m_vc_computed = true;
                return m_vc;
            }
        }

        throw std::runtime_error("vc: failure to converge");
    }
};

/**
 * @brief Search final well of a thrown particle.
 * @param particle Particle (`w` overriden from the start, `v0` is important).
 * @param distribution Type of distribution for `w`, see prrng.
 * @param parameters Parameters for the distribution (defaults appended), see prrng.
 * @param seed `initstate` of the random generator.
 */
inline double throw_particle_Quadratic(
    Quadratic particle,
    enum prrng::distribution distribution,
    const std::vector<double>& parameters,
    uint64_t seed = 0)
{
    prrng::pcg32 gen(seed);
    double w;
    auto param = prrng::default_parameters(distribution, parameters);

    for (size_t i = 0; i < 1e8; ++i) {

        while (true) {
            switch (distribution) {
            case prrng::distribution::random:
                w = gen.random() * param[0] + param[1];
                break;
            case prrng::distribution::exponential:
                w = gen.exponential(param[0]) + param[1];
                break;
            case prrng::distribution::power:
                w = gen.power(param[0]) + param[1];
                break;
            case prrng::distribution::delta:
                w = gen.delta(param[0]) + param[1];
                break;
            case prrng::distribution::pareto:
                w = gen.pareto(param[0], param[1]) + param[2];
                break;
            case prrng::distribution::weibull:
                w = gen.weibull(param[0], param[1]) + param[2];
                break;
            case prrng::distribution::gamma:
                w = gen.gamma(param[0], param[1]) + param[2];
                break;
            case prrng::distribution::normal:
                w = gen.normal(param[0], param[1]) + param[2];
                break;
            case prrng::distribution::custom:
                throw std::runtime_error("Unknow distribution");
            }
        }

        particle.set_w(w);

        if (!particle.exits()) {
            return particle.w();
        }

        particle.set_v0(particle.v_scalar(particle.tau_exit()));
    }

    throw std::runtime_error("throw_particle_Quadratic: failure to stop");
}

/**
 * @brief Search final well of a thrown particles.
 * @param v0 Initial velocity per particle.
 * @param particle Particle (`w` and `v0` are overriden).
 * @param distribution Type of distribution for `w`, see prrng.
 * @param parameters Parameters for the distribution (defaults appended), see prrng.
 * @param seed `initstate` of the first random generator.
 */
template <class T>
inline T throw_particles_Quadratic(
    const T& v0,
    Quadratic particle,
    enum prrng::distribution distribution,
    const std::vector<double>& parameters,
    uint64_t seed = 0)
{
    uint64_t n = static_cast<uint64_t>(v0.size());
    T ret = v0;

    for (uint64_t i = 0; i < n; ++i) {
        particle.set_v0(v0[i]);
        ret[i] = throw_particle_Quadratic(particle, distribution, parameters, seed + i);
    }

    return ret;
}

} // namespace apts

#endif
