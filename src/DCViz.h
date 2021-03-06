#pragma once

#include <assert.h>
#include <thread>
#include <sstream>

using namespace std;

class DCViz {
public:

    DCViz(string fileOrDir);

    inline ~DCViz() {

        if (DCVizThread.joinable()){
            DCVizThread.join();
        }

    }

    inline void launchGUI();
    inline void launch(const bool dynamic=false,
                       const double delay=3,
                       const int figSize_x = -1,
                       const int figSize_y = -1);
    inline void finalize();


private:

    bool launched;
    string fileOrDir;

    thread DCVizThread;

    inline void launchThread(string cmd);

    static void launchPython(string cmd);

};

inline void DCViz::launchGUI() {

    assert(!launched);

    launchThread("import DCVizGUI; DCVizGUI.main('" + (fileOrDir + "')"));

}

inline void
DCViz::launch(const bool dynamic,
              const double delay,
              const int figSize_x,
              const int figSize_y)
{

    assert(delay >= 0);
    assert(!launched);

    stringstream cmd;

    cmd << "import sys; sys.path.append('.'); ";

    cmd << "import DCVizWrapper; DCVizWrapper.main('" << fileOrDir << "'";

    if (dynamic) {
        cmd << ", True";
    } else {
        cmd << ", False";
    }

    cmd << ", delay=" <<delay;

    if ((figSize_x != -1) && (figSize_y != -1)) {
        cmd << ", fs=[" << figSize_x << "," << figSize_y << "]";
    }


    cmd << ")";
    launchThread(cmd.str());

}

//This function is standalone because it's called by both the
//GUI and the non-GUI methods.
inline void DCViz::launchThread(string cmd)
{
    DCVizThread = thread(DCViz::launchPython, cmd);
    launched = true;
}

inline void DCViz::finalize()
{
    DCVizThread.join();
}
