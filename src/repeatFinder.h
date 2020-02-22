

#ifndef PHISPYREPEATFINDER_H
#define PHISPYREPEATFINDER_H

static PyObject * python_input(PyObject *self, PyObject *args);

static PyMethodDef PhiSpyRepeatFinderMethods[] = {
    {"repeatFinder", python_input, METH_VARARGS, "Python interface for C++ repeat finder for PhiSpy"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef PhiSpyRepeatFinderModule = {
    PyModuleDef_HEAD_INIT,
    "repeatFinder",
    "Python for a C++ repeat finder used by PhiSpy to identify potential prophage ends",
    -1,
    PhiSpyRepeatFinderMethods
};

#endif //PHISPYREPEATFINDER_H