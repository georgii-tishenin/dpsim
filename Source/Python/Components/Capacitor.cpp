/** Python binding for Capacitors.
 *
 * @author Steffen Vogel <stvogel@eonerc.rwth-aachen.de>
 * @copyright 2017, Institute for Automation of Complex Power Systems, EONERC
 *
 * CPowerSystems
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************************/

#include "Python/Components/Capacitor.h"

const char *DPsim::Python::Components::DocCapacitor =
"Capacitor(name, node1, node2, capacitance)\n"
"Construct a new Capacitor.\n"
"\n"
"Attributes: ``capacitance``.\n"
"\n"
":param capacitance: Capacitance in Farad.\n"
":returns: A new `Component` representing this Capacitor.\n";

template<>
PyObject* DPsim::Python::Components::Capacitor<CPS::EMT::Ph1::Capacitor>(PyObject* self, PyObject* args)
{
    const char *name;
    double capacitance;

    PyObject *pyNodes;

    if (!PyArg_ParseTuple(args, "sOd", &name, &pyNodes, &capacitance))
        return nullptr;

    try {
        CPS::Node<CPS::Real>::List nodes = Python::Node<CPS::Real>::fromPython(pyNodes);

        Component *pyComp = PyObject_New(Component, &DPsim::Python::ComponentType);
        Component::init(pyComp);
        pyComp->comp = std::make_shared<CPS::EMT::Ph1::Capacitor>(name, nodes, capacitance);

        return (PyObject*) pyComp;
    } catch (...) {
        return nullptr;
    }
}

template<>
PyObject* DPsim::Python::Components::Capacitor<CPS::DP::Ph1::Capacitor>(PyObject* self, PyObject* args)
{
    const char *name;
    double capacitance;

    PyObject *pyNodes;

    if (!PyArg_ParseTuple(args, "sOd", &name, &pyNodes, &capacitance))
        return nullptr;

    try {
        CPS::Node<CPS::Complex>::List nodes = Python::Node<CPS::Complex>::fromPython(pyNodes);

        Component *pyComp = PyObject_New(Component, &DPsim::Python::ComponentType);
        Component::init(pyComp);
        pyComp->comp = std::make_shared<CPS::DP::Ph1::Capacitor>(name, nodes, capacitance);

        return (PyObject*) pyComp;
    } catch (...) {
        return nullptr;
    }
}