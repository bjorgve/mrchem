#pragma once

#include "MRCPP/Parallel"
#include <map>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "mrchem.h"
#include "qmfunctions/qmfunction_fwd.h"

namespace mrchem {

namespace omp {
extern int n_threads;
} // namespace omp

namespace mpi {
extern bool numerically_exact;
extern int shared_memory_size;

extern int world_rank;
extern int world_size;
extern int orb_rank;
extern int orb_size;
extern int share_rank;
extern int share_size;
extern int sh_group_rank;
extern int is_bank;
extern int is_bankclient;
extern int bank_size;
extern std::vector<int> bankmaster;

extern MPI_Comm comm_orb;
extern MPI_Comm comm_share;
extern MPI_Comm comm_sh_group;
extern MPI_Comm comm_bank;

void initialize();
void finalize();
void barrier(MPI_Comm comm);

bool grand_master();
bool share_master();
bool my_orb(const Orbital &orb);
bool my_unique_orb(const Orbital &orb);
void distribute(OrbitalVector &Phi);
void free_foreign(OrbitalVector &Phi);
OrbitalChunk get_my_chunk(OrbitalVector &Phi);

void send_orbital(Orbital &orb, int dst, int tag, MPI_Comm comm = mpi::comm_orb);
void recv_orbital(Orbital &orb, int src, int tag, MPI_Comm comm = mpi::comm_orb);

void send_function(QMFunction &func, int dst, int tag, MPI_Comm comm);
void recv_function(QMFunction &func, int src, int tag, MPI_Comm comm);
void share_function(QMFunction &func, int src, int tag, MPI_Comm comm);

void reduce_function(double prec, QMFunction &func, MPI_Comm comm);
void broadcast_function(QMFunction &func, MPI_Comm comm);

void allreduce_vector(IntVector &vec, MPI_Comm comm);
void allreduce_vector(DoubleVector &vec, MPI_Comm comm);
void allreduce_vector(ComplexVector &vec, MPI_Comm comm);
void allreduce_matrix(IntMatrix &vec, MPI_Comm comm);
void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm);
void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm);

} // namespace mpi

namespace bank {

struct deposit {
    Orbital *orb;
    double *data; // for pure data arrays
    bool hasdata;
    int datasize;
    int id;     // to identify what is deposited
    int source; // mpi rank from the source of the data
};

struct queue_struct {
    int id;
    std::vector<int> clients;
};

} // namespace bank

class Bank {
public:
    Bank() = default;
    ~Bank();
    void open();
    void close();
    void clear_all(int i, MPI_Comm comm);
    void clear(int ix);
    int put_orb(int id, Orbital &orb);
    int get_orb(int id, Orbital &orb, int wait = 0);
    int get_orb_del(int id, Orbital &orb);
    int put_func(int id, QMFunction &func);
    int get_func(int id, QMFunction &func);
    int set_datasize(int datasize);
    int put_data(int id, int size, double *data);
    int get_data(int id, int size, double *data);

private:
    int const CLOSE_BANK = 1;
    int const CLEAR_BANK = 2;
    int const GET_ORBITAL = 3;
    int const GET_ORBITAL_AND_WAIT = 4;
    int const GET_ORBITAL_AND_DELETE = 5;
    int const SAVE_ORBITAL = 6;
    int const GET_FUNCTION = 7;
    int const SAVE_FUNCTION = 8;
    int const SET_DATASIZE = 9;
    int const GET_DATA = 10;
    int const SAVE_DATA = 11;
    std::map<int, int> id2ix;
    std::vector<bank::deposit> deposits;
    std::map<int, int> id2qu;
    std::vector<bank::queue_struct> queue;

    void clear_bank();
};

namespace mpi {
extern Bank orb_bank;
} // namespace mpi

} // namespace mrchem
