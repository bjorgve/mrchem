#include "MRCPP/Printer"

#include "parallel.h"
#include "qmfunctions/ComplexFunction.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"

using mrcpp::Printer;

namespace mrchem {

namespace omp {

int n_threads = omp_get_max_threads();

} // namespace omp

namespace mpi {

bool numerically_exact = false;
int shared_memory_size = 1000;

int world_size = 0;
int world_rank = 1;
int orb_rank = 0;
int orb_size = 1;
int share_rank = 0;
int share_size = 1;
int sh_group_rank = 0;
int is_bank = 0;
int is_bankclient = 1;
int is_bankmaster = 0; // only one bankmaster is_bankmaster
int bank_size = 0;
std::vector<int> bankmaster;

MPI_Comm comm_orb;
MPI_Comm comm_share;
MPI_Comm comm_sh_group;
MPI_Comm comm_bank;

Bank orb_bank;

} // namespace mpi

int const func_id_shift = 10000; // to ensure that orbitals and function do not collide

void mpi::initialize() {
    omp_set_dynamic(0);

#ifdef HAVE_MPI
    MPI_Init(nullptr, nullptr);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // divide the world into groups
    // each group has its own group communicator definition

    // define independent group of MPI processes, that are not part of comm_orb
    // for now the new group does not include comm_share
    comm_bank = MPI_COMM_WORLD; // clients and master
    MPI_Comm comm_remainder;    // clients only

    bankmaster.resize(bank_size);
    for (int i = 0; i < bank_size; i++) {
        bankmaster[i] = world_size - i - 1; // rank of the bankmasters
    }
    if (world_rank < world_size - bank_size) {
        // everything which is left
        is_bank = 0;
        is_bankclient = 1;
    } else {
        // special group of bankmasters
        is_bank = 1;
        is_bankclient = 0;
        if (world_rank == world_size - bank_size) is_bankmaster = 1;
    }
    MPI_Comm_split(MPI_COMM_WORLD, is_bankclient, world_rank, &comm_remainder);

    // split world into groups that can share memory
    MPI_Comm_split_type(comm_remainder, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_share);

    MPI_Comm_rank(comm_share, &share_rank);
    MPI_Comm_size(comm_share, &share_size);

    // define a rank of the group
    MPI_Comm_split(comm_remainder, share_rank, world_rank, &comm_sh_group);
    // mpiShRank is color (same color->in same group)
    // MPI_worldrank is key (orders rank within the groups)

    // we define a new orbital rank, so that the orbitals within
    // a shared memory group, have consecutive ranks
    MPI_Comm_rank(comm_sh_group, &sh_group_rank);

    orb_rank = share_rank + sh_group_rank * world_size;
    MPI_Comm_split(comm_remainder, 0, orb_rank, &comm_orb);
    // 0 is color (same color->in same group)
    // mpiOrbRank is key (orders rank in the group)

    MPI_Comm_rank(comm_orb, &orb_rank);
    MPI_Comm_size(comm_orb, &orb_size);
#endif
}

void mpi::finalize() {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

void mpi::barrier(MPI_Comm comm) {
#ifdef HAVE_MPI
    MPI_Barrier(comm);
#endif
}

/*********************************
 * Orbital related MPI functions *
 *********************************/

bool mpi::grand_master() {
    return (mpi::orb_rank == 0) ? true : false;
}

bool mpi::share_master() {
    return (mpi::share_rank == 0) ? true : false;
}

/** @brief Test if orbital belongs to this MPI rank (or is common)*/
bool mpi::my_orb(const Orbital &orb) {
    return (orb.rankID() < 0 or orb.rankID() == mpi::orb_rank) ? true : false;
}

/** @brief Test if orbital belongs to this MPI rank */
bool mpi::my_unique_orb(const Orbital &orb) {
    return (orb.rankID() == mpi::orb_rank) ? true : false;
}

/** @brief Distribute orbitals in vector round robin. Orbitals should be empty.*/
void mpi::distribute(OrbitalVector &Phi) {
    for (int i = 0; i < Phi.size(); i++) Phi[i].setRankID(i % mpi::orb_size);
}

/** @brief Free all function pointers not belonging to this MPI rank */
void mpi::free_foreign(OrbitalVector &Phi) {
    for (auto &i : Phi) {
        if (not mpi::my_orb(i)) i.free(NUMBER::Total);
    }
}

/** @brief Return the subset of an OrbitalVector that belongs to this MPI rank */
OrbitalChunk mpi::get_my_chunk(OrbitalVector &Phi) {
    OrbitalChunk chunk;
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) chunk.push_back(std::make_tuple(i, Phi[i]));
    }
    return chunk;
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(IntVector &vec, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(DoubleVector &vec, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the vector with contributions from all MPI ranks */
void mpi::allreduce_vector(ComplexVector &vec, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = vec.size();
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), N, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(IntMatrix &mat, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_INT, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_DOUBLE, MPI_SUM, comm);
#endif
}

/** @brief Add up each entry of the matrix with contributions from all MPI ranks */
void mpi::allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm) {
#ifdef HAVE_MPI
    int N = mat.size();
    MPI_Allreduce(MPI_IN_PLACE, mat.data(), N, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
#endif
}

// send an orbital with MPI, includes orbital meta data
void mpi::send_orbital(Orbital &orb, int dst, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    mpi::send_function(orb, dst, tag, comm);
    OrbitalData &orbinfo = orb.getOrbitalData();
    MPI_Send(&orbinfo, sizeof(OrbitalData), MPI_BYTE, dst, 0, comm);
#endif
}

// receive an orbital with MPI, includes orbital meta data
void mpi::recv_orbital(Orbital &orb, int src, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    mpi::recv_function(orb, src, tag, comm);

    MPI_Status status;
    OrbitalData &orbinfo = orb.getOrbitalData();
    MPI_Recv(&orbinfo, sizeof(OrbitalData), MPI_BYTE, src, 0, comm, &status);
#endif
}

// send a function with MPI
void mpi::send_function(QMFunction &func, int dst, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    if (func.isShared()) MSG_WARN("Sending a shared function is not recommended");
    FunctionData &funcinfo = func.getFunctionData();
    MPI_Send(&funcinfo, sizeof(FunctionData), MPI_BYTE, dst, 0, comm);
    if (func.hasReal()) mrcpp::send_tree(func.real(), dst, tag, comm, funcinfo.real_size);
    if (func.hasImag()) mrcpp::send_tree(func.imag(), dst, tag + 10000, comm, funcinfo.imag_size);
#endif
}

// receive a function with MPI
void mpi::recv_function(QMFunction &func, int src, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    if (func.isShared()) MSG_WARN("Receiving a shared function is not recommended");
    MPI_Status status;

    FunctionData &funcinfo = func.getFunctionData();
    MPI_Recv(&funcinfo, sizeof(FunctionData), MPI_BYTE, src, 0, comm, &status);
    if (funcinfo.real_size > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasReal()) func.alloc(NUMBER::Real);
        mrcpp::recv_tree(func.real(), src, tag, comm, funcinfo.real_size);
    }

    if (funcinfo.imag_size > 0) {
        // We must have a tree defined for receiving nodes. Define one:
        if (not func.hasImag()) func.alloc(NUMBER::Imag);
        mrcpp::recv_tree(func.imag(), src, tag + 10000, comm, funcinfo.imag_size);
    }
#endif
}

/** Update a shared function after it has been changed by one of the MPI ranks. */
void mpi::share_function(QMFunction &func, int src, int tag, MPI_Comm comm) {
#ifdef HAVE_MPI
    if (func.isShared()) {
        if (func.hasReal()) mrcpp::share_tree(func.real(), src, tag, comm);
        if (func.hasImag()) mrcpp::share_tree(func.imag(), src, 2 * tag, comm);
    }
#endif
}

/** @brief Add all mpi function into rank zero */
void mpi::reduce_function(double prec, QMFunction &func, MPI_Comm comm) {
/* 1) Each odd rank send to the left rank
   2) All odd ranks are "deleted" (can exit routine)
   3) new "effective" ranks are defined within the non-deleted ranks
      effective rank = rank/fac , where fac are powers of 2
   4) repeat
 */
#ifdef HAVE_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) {
        if ((comm_rank / fac) % 2 == 0) {
            // receive
            int src = comm_rank + fac;
            if (src < comm_size) {
                QMFunction func_i(false);
                int tag = 3333 + src;
                mpi::recv_function(func_i, src, tag, comm);
                func.add(1.0, func_i); // add in place using union grid
                func.crop(prec);
            }
        }
        if ((comm_rank / fac) % 2 == 1) {
            // send
            int dest = comm_rank - fac;
            if (dest >= 0) {
                int tag = 3333 + comm_rank;
                mpi::send_function(func, dest, tag, comm);
                break; // once data is sent we are done
            }
        }
        fac *= 2;
    }
    MPI_Barrier(comm);
#endif
}

/** @brief Distribute rank zero function to all ranks */
void mpi::broadcast_function(QMFunction &func, MPI_Comm comm) {
/* use same strategy as a reduce, but in reverse order */
#ifdef HAVE_MPI
    int comm_size, comm_rank;
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    if (comm_size == 1) return;

    int fac = 1; // powers of 2
    while (fac < comm_size) fac *= 2;
    fac /= 2;

    while (fac > 0) {
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 1) {
            // receive
            int src = comm_rank - fac;
            int tag = 4334 + comm_rank;
            mpi::recv_function(func, src, tag, comm);
        }
        if (comm_rank % fac == 0 and (comm_rank / fac) % 2 == 0) {
            // send
            int dst = comm_rank + fac;
            int tag = 4334 + dst;
            if (dst < comm_size) mpi::send_function(func, dst, tag, comm);
        }
        fac /= 2;
    }
    MPI_Barrier(comm);
#endif
}

Bank::Bank() {}
Bank::~Bank() {
    for (int ix = 1; ix < this->deposits.size(); ix++) this->clear(ix);
}

int Bank::open() {
#ifdef HAVE_MPI
    MPI_Status status;
    char safe_data1;
    int deposit_size = sizeof(deposit);
    int n_chunks, ix;
    int message;

    deposits.resize(1); // we reserve 0, since it is the value returned for undefined key
    queue.resize(1);    // we reserve 0, since it is the value returned for undefined key

    // The bank never goes out of this loop as long as it receives a close signal!
    while (true) {
        MPI_Recv(&message, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, mpi::comm_bank, &status);
        if (message == CLOSE_BANK) {
            if (mpi::is_bankmaster) std::cout << "Bank is closing" << std::endl;
            this->clear_bank();
            break; // close bank, i.e stop listening for incoming messages
        }
        if (message == CLEAR_BANK) {
            this->clear_bank();
            // send message that it is ready (value of message is not used)
            MPI_Send(&message, 1, MPI_INTEGER, status.MPI_SOURCE, 77, mpi::comm_bank);
        }
        if (message == GET_ORBITAL or message == GET_FUNCTION) {
            // withdrawal
            int ix = id2ix[status.MPI_TAG];
            if (ix == 0) {
                // the id does not exist. Put in queue and Wait until it is defined
                if (id2qu[status.MPI_TAG] == 0) {
                    queue.push_back({status.MPI_TAG, {status.MPI_SOURCE}});
                    id2qu[status.MPI_TAG] = queue.size() - 1;
                } else {
                    // somebody is already waiting for this id. queue in queue
                    queue[id2qu[status.MPI_TAG]].clients.push_back(status.MPI_SOURCE);
                }
            } else {
                if (deposits[ix].id != status.MPI_TAG) std::cout << ix << " Bank acounting error " << std::endl;
                if (message == GET_ORBITAL) {
                    mpi::send_orbital(*deposits[ix].orb, status.MPI_SOURCE, deposits[ix].id, mpi::comm_bank);
                }
                if (message == GET_FUNCTION) {
                    mpi::send_function(*deposits[ix].orb, status.MPI_SOURCE, deposits[ix].id, mpi::comm_bank);
                }
            }
        }
        if (message == SAVE_ORBITAL or message == SAVE_FUNCTION) {
            // make a new deposit
            if (id2ix[status.MPI_TAG]) {
                ix = id2ix[status.MPI_TAG];
            } else {
                ix = deposits.size(); // NB: ix is now index of last element + 1
                deposits.resize(ix + 1);
            }
            deposits[ix].orb = new Orbital(0);
            deposits[ix].id = status.MPI_TAG;
            id2ix[deposits[ix].id] = ix;
            deposits[ix].source = status.MPI_SOURCE;
            if (message == SAVE_ORBITAL) {
                mpi::recv_orbital(*deposits[ix].orb, deposits[ix].source, status.MPI_TAG, mpi::comm_bank);
            }
            if (message == SAVE_FUNCTION) {
                mpi::recv_function(*deposits[ix].orb, deposits[ix].source, status.MPI_TAG, mpi::comm_bank);
            }
            if (id2qu[deposits[ix].id] != 0) {
                // someone is waiting for those data. Send to them
                int iq = id2qu[deposits[ix].id];
                if (deposits[ix].id != queue[iq].id) std::cout << ix << " Bank queue acounting error " << std::endl;
                for (int iqq : queue[iq].clients) {
                    if (message == SAVE_ORBITAL) {
                        mpi::send_orbital(*deposits[ix].orb, iqq, queue[iq].id, mpi::comm_bank);
                    }
                    if (message == SAVE_FUNCTION) {
                        mpi::send_function(*deposits[ix].orb, iqq, queue[iq].id, mpi::comm_bank);
                    }
                }
                queue.erase(queue.begin() + iq);
                id2qu.erase(deposits[ix].id);
            }
        }
    }
#endif
}

// save orbital in Bank with identity id
int Bank::put_orb(int id, Orbital &orb) {
#ifdef HAVE_MPI
    // for now we distribute according to id
    if (id > func_id_shift) MSG_WARN("Bank id should be less than func_id_shift (10000)");
    MPI_Send(&SAVE_ORBITAL, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    mpi::send_orbital(orb, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
#endif
}

// get orbital with identity id
int Bank::get_orb(int id, Orbital &orb) {
#ifdef HAVE_MPI
    MPI_Status status;
    MPI_Send(&GET_ORBITAL, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    mpi::recv_orbital(orb, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
#endif
}

// save function in Bank with identity id . NB:not tested
int Bank::put_func(int id, QMFunction &func) {
#ifdef HAVE_MPI
    // for now we distribute according to id
    id += func_id_shift;
    MPI_Send(&SAVE_FUNCTION, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    mpi::send_function(func, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
#endif
}

// get function with identity id
int Bank::get_func(int id, QMFunction &func) {
#ifdef HAVE_MPI
    MPI_Status status;
    MPI_Send(&GET_FUNCTION, 1, MPI_INTEGER, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
    id += func_id_shift;
    mpi::recv_function(func, mpi::bankmaster[id % mpi::bank_size], id, mpi::comm_bank);
#endif
}

// Ask to close the Bank
int Bank::close() {
#ifdef HAVE_MPI
    for (int i = 0; i < mpi::bank_size; i++) {
        MPI_Send(&CLOSE_BANK, 1, MPI_INTEGER, mpi::bankmaster[i], 0, mpi::comm_bank);
    }
#endif
}

// remove all deposits
// NB:: collective call. All clients must call this
int Bank::clear_all(int i, MPI_Comm comm) {
#ifdef HAVE_MPI
    // 1) wait until all clients are ready
    mpi::barrier(comm);
    // master send signal to bank
    if (i == 0) {
        for (int i = 0; i < mpi::bank_size; i++) {
            MPI_Send(&CLEAR_BANK, 1, MPI_INTEGER, mpi::bankmaster[i], 0, mpi::comm_bank);
            ;
        }
        for (int i = 0; i < mpi::bank_size; i++) {
            // wait until Bank is finished and has sent signal
            MPI_Status status;
            int message;
            MPI_Recv(&message, 1, MPI_INTEGER, mpi::bankmaster[i], 77, mpi::comm_bank, &status);
        }
    }
    mpi::barrier(comm);
#endif
}

int Bank::clear_bank() {
#ifdef HAVE_MPI
    for (int ix = 1; ix < this->deposits.size(); ix++) this->clear(ix);
    this->deposits.resize(1);
    this->queue.resize(1);
    this->id2ix.clear();
    this->id2qu.clear();
#endif
}

int Bank::clear(int ix) {
#ifdef HAVE_MPI
    delete deposits[ix].orb;
#endif
}

} // namespace mrchem
