#include "DCViz.h"

#include <Python.h>

DCViz::DCViz(string fileOrDir) :
    launched(false),
    fileOrDir(fileOrDir)
{

}

void DCViz::launchPython(std::string cmd)
{

    Py_Initialize();

    PyRun_SimpleString(cmd.c_str());

    Py_Finalize();

}
