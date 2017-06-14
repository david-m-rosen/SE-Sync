#include "SESync.h"
#include "SESync_utils.h"

using namespace std;
using namespace SESync;

int main(int argc, char** argv)
{
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " [input .g2o file]" << endl;
        exit(1);
    }

    size_t num_poses;
    vector<SESync::RelativePoseMeasurement> measurements = read_g2o_file(argv[1], num_poses);
    cout << "Loaded " << measurements.size() << " measurements between "
         << num_poses << " poses from file " << argv[1] << endl
         << endl;

    SESyncOpts opts;
    //opts.tolgradnorm = 0;
    //opts.rel_func_decrease_tol = 1e-10;
    //opts.use_chordal_initialization = false;
    //opts.use_Cholesky = false;
    //opts.eig_comp_tol = 1e-12;
    //opts.num_Lanczos_vectors = 50;
    opts.verbose = true; // Print output to stdout

    SESyncResult results = SESync::SESync(measurements, opts);
}
