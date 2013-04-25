#include <sys/stat.h> // mkdir()
#include "kmerizer.h"

int main(int argc, char *argv[])
{
    // parse args
    if (argc != 4) {
        fprintf(stdout, "Usage: %s <k> <threads> <input dir>\n", argv[0]);
        return 1;
    }
    size_t k = atoi(argv[1]);
    size_t threads = atoi(argv[2]);
    char* inprefix = argv[3];

    Kmerizer *counter = new Kmerizer(k, threads, inprefix, CANONICAL);

    cout << "Load" << endl;
    counter->load();
    counter->histogram();
    return 0;
}