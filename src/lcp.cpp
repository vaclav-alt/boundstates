#include <iostream>
#include "lcp-app.hpp"

int main(int argc, const char* argv[]) {
    try {
        LcpApp app;
        app.Exec(argc, argv);
        return 0;
    }
    catch(std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    catch(std::exception * e) {
        std::cout << "ERROR: " << e->what() << std::endl;
        return 1;
    }

    catch(...) {
        std::cout << "ERROR" << std::endl;
        return 1;
    }

	return 0;
};
