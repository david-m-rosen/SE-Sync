#include "SESync/SESync.h"
#include "SESync/SESync_utils.h"

#include <fstream>

using namespace std;
using namespace SESync;

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " [input .g2o file]" << endl;
    exit(1);
  }

  size_t num_poses;
  vector<SESync::RelativePoseMeasurement> measurements =
      read_g2o_file(argv[1], num_poses);
  cout << "Loaded " << measurements.size() << " measurements between "
       << num_poses << " poses from file " << argv[1] << endl
       << endl;
  if (measurements.size() == 0) {
    cout << "Error: No measurements were read!"
         << " Are you sure the file exists?" << endl;
    exit(1);
  }

  SESyncOpts opts;
  opts.verbose = true; // Print output to stdout
  opts.num_threads = 4;

  SESyncResult results = SESync::SESync(measurements, opts);

  string filename = "poses.txt";

  cout << "Saving final poses to file: " << filename << endl;

  ofstream poses_file(filename);
  poses_file << results.xhat;
  poses_file.close();
}
