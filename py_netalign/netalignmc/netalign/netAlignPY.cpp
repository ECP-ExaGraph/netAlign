#include <pybind11/pybind11.h>
#include "netAlignKernel.h"

namespace py = pybind11;

PYBIND11_MODULE(netAlignPY, m) {
    m.doc() = "Network Alignment plugin";
    m.def("netalign", &netAlign, "A function which aligns two netowrks",py::arg("argc") = 2, py::arg("alpha")= 1);
}
