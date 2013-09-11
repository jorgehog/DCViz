#include <Python.h>
#include <thread>
#include <sstream>

using namespace std;

class DCViz {
public:

    DCViz(string fileOrDir) : launched(false), fileOrDir(fileOrDir) {}

    inline ~DCViz() {

        if (DCVizThread.joinable()){
            DCVizThread.join();
        }

        DCVizThread.detach();
    }

    inline void launchGUI(); //Not implemented
    inline void launch(const bool dynamic=false, const double delay=3);
    inline void finalize();


private:

    bool launched;
    string fileOrDir;

    thread DCVizThread;

    inline void launchThread(string cmd);

    inline static void launchPython(string cmd);

};

inline void DCViz::launchGUI() {

    assert(!launched);

    launchThread("import DCVizGUI; DCVizGUI.main('" + (fileOrDir + "')"));

}


inline void DCViz::launch(const bool dynamic, const double delay) {

    assert(delay >= 0);
    assert(!launched);

    stringstream cmd;

    cmd << "import DCVizWrapper; DCVizWrapper.main('" <<fileOrDir << "', ";

    if (dynamic) {
        cmd << "True, ";
    } else {
        cmd << "False, ";
    }

    cmd << delay << ")";

    launchThread(cmd.str());

}

//This function is standalone because it's called by both the
//GUI and the non-GUI methods.
inline void DCViz::launchThread(string cmd) {
    DCVizThread = thread(DCViz::launchPython, cmd);
    launched = true;
}

inline void DCViz::finalize() {
    DCVizThread.join();
}

inline void DCViz::launchPython(string cmd) {

    Py_Initialize();

    PyRun_SimpleString(cmd.c_str());

    Py_Finalize();

}

















