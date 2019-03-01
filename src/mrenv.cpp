#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "getkw/Getkw.hpp"

#include <algorithm>

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

using namespace std;
using namespace mrcpp;
using namespace mrchem;

namespace mrenv {

template <int D> auto *initialize_world(int min_scale, std::array<int, D> &corner, std::array<int, D> &boxes,
                                        std::array<double, D> sf, bool periodic) {

    if (periodic) return new BoundingBox<D>(sf, periodic);
    return new BoundingBox<D>(min_scale, corner, boxes);
}

void initialize(int argc, char **argv) {
    const char *infile = nullptr;
    if (argc == 1) {
        infile = "STDIN";
    } else if (argc == 2) {
        infile = argv[1];
    } else {
        MSG_ERROR("Ivalid number of arguments!");
    }

    // Initialize input
    input = Getkw(infile, false, true);

    // Initialize printing
    int printlevel = input.get<int>("printlevel");
    bool teletype = input.get<bool>("teletype");
    if (teletype) {
        Printer::init(printlevel, mpi::orb_rank, mpi::orb_size, "mrchem");
    } else {
        Printer::init(printlevel, mpi::orb_rank, mpi::orb_size);
    }
    Printer::setPrecision(15);

    mpi::numerically_exact = input.get<bool>("mpi.numerically_exact");
    mpi::share_nuc_pot = input.get<bool>("mpi.share_nuclear_potential");
    mpi::share_coul_dens = input.get<bool>("mpi.share_coulomb_density");
    mpi::share_coul_pot = input.get<bool>("mpi.share_coulomb_potential");
    mpi::share_xc_dens = input.get<bool>("mpi.share_xc_density");
    mpi::share_xc_pot = input.get<bool>("mpi.share_xc_potential");
    mpi::shared_memory_size = input.get<int>("mpi.shared_memory_size");

    // Initialize world box
    int min_scale = input.get<int>("mra.min_scale");
    int max_scale = input.get<int>("mra.max_scale");
    vector<int> corner = input.getIntVec("mra.corner");
    vector<int> boxes = input.getIntVec("mra.boxes");
    auto scaling_factor = input.getDblVec("mra.scaling_factor");
    auto periodic = input.get<bool>("mra.periodic");
    std::array<int, 3> c_idx;
    std::array<int, 3> n_bxs;
    std::array<double, 3> sf;
    std::copy_n(corner.begin(), 3, c_idx.begin());
    std::copy_n(boxes.begin(), 3, n_bxs.begin());
    std::copy_n(scaling_factor.begin(), 3, sf.begin());
    BoundingBox<3> *world = initialize_world<3>(min_scale, c_idx, n_bxs, sf, periodic);

    // Initialize scaling basis
    int order = input.get<int>("mra.order");
    string btype = input.get<string>("mra.basis_type");

    int max_depth = max_scale - min_scale;
    if (min_scale < MinScale) MSG_FATAL("Root scale too large");
    if (max_scale > MaxScale) MSG_FATAL("Max scale too large");
    if (max_depth > MaxDepth) MSG_FATAL("Max depth too large");

    // Initialize global MRA
    if (btype == "i") {
        InterpolatingBasis basis(order);
        MRA = new MultiResolutionAnalysis<3>(*world, basis, max_depth);
    } else if (btype == "l") {
        LegendreBasis basis(order);
        MRA = new MultiResolutionAnalysis<3>(*world, basis, max_depth);
    } else {
        MSG_FATAL("Invalid basis type!");
    }

    println(0, endl << endl);
    println(0, "************************************************************");
    println(0, "***     __  __ ____   ____ _                             ***");
    println(0, "***    |  \\/  |  _ \\ / ___| |__   ___ _ __ ___           ***");
    println(0, "***    | |\\/| | |_) | |   | '_ \\ / _ \\ '_ ` _ \\          ***");
    println(0, "***    | |  | |  _ <| |___| | | |  __/ | | | | |         ***");
    println(0, "***    |_|  |_|_| \\_\\\\____|_| |_|\\___|_| |_| |_|         ***");
    println(0, "***                                                      ***");
    println(0, "***    VERSION " << PROGRAM_VERSION << " (rev. " << GIT_REVISION << ")                     ***");
    println(0, "***                                                      ***");
    println(0, "***    Stig Rune Jensen <stig.r.jensen@uit.no>           ***");
    println(0, "***    Luca Frediani    <luca.frediani@uit.no>           ***");
    println(0, "***    Peter Wind       <peter.wind@uit.no>              ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, endl);

    if (mpi::orb_size > 1 or omp::n_threads > 1) {
        println(0, "+++ Parallel execution: ");
        println(0, "  MPI hosts available     : " << mpi::orb_size);
        println(0, "  Threads/host            : " << omp::n_threads);
        println(0, "  Total used CPUs         : " << mpi::orb_size * omp::n_threads);
        println(0, "");
    } else {
        println(0, "+++ Serial execution" << endl);
    }

    Printer::printEnvironment();

    int oldprec = Printer::setPrecision(6);
    MRA->print();
    Printer::setPrecision(oldprec);
}

void finalize(double wt) {
    // Delete global MRA
    if (MRA != nullptr) delete MRA;
    MRA = nullptr;

    Printer::setPrecision(6);
    println(0, endl);
    println(0, "************************************************************");
    println(0, "***                                                      ***");
    println(0, "***                    Exiting MRChem                    ***");
    println(0, "***                                                      ***");
    println(0, "***               Wall time: " << wt << "                ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, endl);
}

} // namespace mrenv
