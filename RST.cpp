#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

int main() {
    std::ifstream fin("files/M6Tw025_Stat.dat");   // <-- your input file
    if (!fin) {
        std::cerr << "Error: could not open input.dat\n";
        return 1;
    }

    std::string line;
    int vel_file_offset = 142;   // example: how many header lines to skip
    int Ny = 0;                  // will count how many lines we keep

    // Skip header lines
    for (int i = 0; i < vel_file_offset; ++i) {
        std::getline(fin, line);
    }

    // Prepare storage
    std::vector<double> z, zd, up, vp, wp, uvp;

    // Read the rest
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::vector<double> values;
        double val;
        while (iss >> val) {
            values.push_back(val);
        }

        // Grab only the needed columns:
        // y_d[count]    = values[1]  → call this zd
        // urms_us       = values[8]  → up
        // wrms_us       = values[9]  → wp
        // vrms_us       = values[10] → vp
        // uvrms_us      = values[15] → uvp
        if (values.size() > 15) {
            z.push_back(values[0]); // just index as "z"
            zd.push_back(values[1]);
            up.push_back(values[8]);
            vp.push_back(values[10]);
            wp.push_back(values[9]);
            uvp.push_back(values[15]);
        }
    }

    Ny = static_cast<int>(z.size());

    // Write output
    std::ofstream fout("RST.dat");
    if (!fout) {
        std::cerr << "Error: could not write RST.dat\n";
        return 1;
    }

    fout << "VARIABLES= y, yd, up, vp, wp, uvp\n";
    fout << "ZONE f=point, i=" << Ny << "\n";

    for (int i = 0; i < Ny; ++i) {
        fout << z[i]   << " "
             << zd[i]  << " "
             << up[i]  << " "
             << vp[i]  << " "
             << wp[i]  << " "
             << uvp[i] << "\n";
    }

    std::cout << "Wrote " << Ny << " rows to RST.dat\n";
    return 0;
}
