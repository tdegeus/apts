#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>
#include <xtensor-python/xtensor_python_config.hpp> // todo: remove for xtensor-python >0.26.1

#include <apts.h>

namespace py = pybind11;

/**
 * Overrides the `__name__` of a module.
 * Classes defined by pybind11 use the `__name__` of the module as of the time they are defined,
 * which affects the `__repr__` of the class type objects.
 */
class ScopedModuleNameOverride {
public:
    explicit ScopedModuleNameOverride(py::module m, std::string name) : module_(std::move(m))
    {
        original_name_ = module_.attr("__name__");
        module_.attr("__name__") = name;
    }
    ~ScopedModuleNameOverride()
    {
        module_.attr("__name__") = original_name_;
    }

private:
    py::module module_;
    py::object original_name_;
};

PYBIND11_MODULE(_apts, m)
{
    // Ensure members to display as `apts.X` rather than `apts._apts.X`
    ScopedModuleNameOverride name_override(m, "apts");

    xt::import_numpy();

    m.doc() = "Portable Reconstructible Random Number Generator";

    m.def(
        "version",
        &apts::version,
        "Return version string. "
        "See :cpp:func:`apts::version`.");

    {
        using S = apts::Quadratic;

        py::class_<S> cls(m, "Quadratic");

        cls.def(
            py::init<double, double, double, double, double>(),
            "Constructor.",
            py::arg("v0") = 0,
            py::arg("w") = 1,
            py::arg("m") = 1,
            py::arg("eta") = 0.1,
            py::arg("mu") = 1);

        cls.def("__repr__", [](const S&) { return "<apts.Quadratic>"; });

        cls.def_property("v0", &S::v0, &S::set_v0);
        cls.def_property_readonly("exits", &S::exits);
        cls.def_property_readonly("tau_exit", &S::tau_exit);
        cls.def_property_readonly("vc", &S::vc);
        cls.def("r", &S::r<xt::pyarray<double>>, py::arg("tau"));
        cls.def("v", &S::v<xt::pyarray<double>>, py::arg("tau"));
    }

} // PYBIND11_MODULE
