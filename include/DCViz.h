#include <stdlib.h>
#include <thread>
#include <Python.h>
#include <sstream>

using namespace std;

class DCViz {
public:

    DCViz(string fileOrDir) : launched(false), fileOrDir(fileOrDir) {}
    inline ~DCViz() {
        DCVizThread.detach();
    }

    inline void launchGUI() {assert(false);} //Not implemented
    inline void launch(const bool dynamic=false, const int delay=3);
    inline void join();


private:

    bool launched;
    string fileOrDir;

    thread DCVizThread;

    inline void launchThread(const string cmd);
    inline static void pyLaunch(const string cmd);
};


inline void DCViz::launch(const bool dynamic, const int delay) {

    assert(delay >= 0);
    assert(!launched);

    stringstream s;

    s << "import DCVizWrapper; DCVizWrapper.main('" <<fileOrDir << "', ";

    if (dynamic) {
        s << "True, ";
    } else {
        s << "False, ";
    }

    s << delay << ")";

    launchThread(s.str());
}


inline void DCViz::launchThread(const string cmd) {
    DCVizThread = thread(DCViz::pyLaunch, cmd);
    launched = true;
}

inline void DCViz::join() {
    DCVizThread.join();
}

inline void DCViz::pyLaunch(const string cmd) {
    Py_Initialize();
    PyRun_SimpleString(cmd.c_str());
    Py_Finalize();
}
