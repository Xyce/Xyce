//-------------------------------------------------------------------------
//   Copyright 2002-2022 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//
// Purpose        : 
//                  
//                  
//
// Special Notes  : 
//                  
//
// Creator        : David Baur
//
// Creation Date  : 
//
//
//
//
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#ifdef Xyce_PARALLEL_MPI

#include <sstream>

#include <N_PDS_MPI.h>
#include <N_UTL_Op.h>
#include <N_UTL_Marshal.h>

namespace Xyce {
namespace Parallel {

template struct Loc<int>;
template struct Loc<double>;
template struct Loc<float>;

MPI_Datatype
double_complex_type()
{
  static MPI_Datatype s_mpi_double_complex;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Type_contiguous(2, MPI_DOUBLE, &s_mpi_double_complex);
    MPI_Type_commit(&s_mpi_double_complex);
  }
  return s_mpi_double_complex;
}

MPI_Datatype
float_complex_type()
{
  static MPI_Datatype s_mpi_float_complex;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Type_contiguous(2, MPI_FLOAT, &s_mpi_float_complex);
    MPI_Type_commit(&s_mpi_float_complex);
  }
  return s_mpi_float_complex;
}


// #ifdef MPI_LONG_LONG_INT
// MPI_Datatype
// long_long_int_int_type()
// {
//   static MPI_Datatype s_mpi_long_long_int_int;
//   static bool initialized = false;

//   int B[] = {2, 1};
//   MPI_Aint D[] = {0, 8};
//   MPI_Datatype T[] = {MPI_LONG_LONG_INT, MPI_INT};
  
//   if (!initialized) {
//     initialized = true;

//     MPI_Type_struct(2, B, D, T, &s_mpi_long_long_int_int);
//     MPI_Type_commit(&s_mpi_long_long_int_int);
//   }
//   return s_mpi_long_long_int_int;
// }
// #endif


MPI_Datatype
double_double_int_type()
{
  static MPI_Datatype s_mpi_double_double_int;
  static bool initialized = false;

  int B[] = {2, 1};
  MPI_Aint D[] = {0, 16};
  MPI_Datatype T[] = {MPI_DOUBLE, MPI_INT};

  if (!initialized) {
    initialized = true;

    MPI_Type_create_struct(2, B, D, T, &s_mpi_double_double_int);
    MPI_Type_commit(&s_mpi_double_double_int);
  }
  return s_mpi_double_double_int;
}


namespace {

extern "C" {
  void
  mpi_double_complex_sum(
    void *		invec,
    void *		inoutvec,
    int *		len,
    MPI_Datatype *	datatype)
  {
    std::complex<double> *complex_in = static_cast<std::complex<double> *>(invec);
    std::complex<double> *complex_inout = static_cast<std::complex<double> *>(inoutvec);

    for (int i = 0; i < *len; ++i)
      complex_inout[i] += complex_in[i];
  }

  void
  mpi_op_identifier_compare(
    void *                invec,
    void *                inoutvec,
    int *                 len,
    MPI_Datatype *        datatype)
  {
    Util::Op::Identifier *identifier_in = static_cast<Util::Op::Identifier *>(invec);
    Util::Op::Identifier *identifier_inout = static_cast<Util::Op::Identifier *>(inoutvec);

    // if (*datatype != Parallel::Datatype<Util::Op::Identifier>()::type())
    //   return;

    for (int i = 0; i < *len; ++i) {
      if (identifier_in[i] != identifier_inout[i]) { // if equal do nothing
        if (identifier_inout[i] == Util::Op::identifier<Util::Op::UndefinedOp>()) // if inout is undefined, set it to in
          identifier_inout[i] = identifier_in[i];
        else if (identifier_in[i] != Util::Op::identifier<Util::Op::UndefinedOp>()) // if in is not undefined, the set inout to 0
          identifier_inout[i] = 0;
      }
    }
  }
} // extern "C"

} // namespace <unnamed>


MPI_Op
double_complex_sum_op()
{
  static MPI_Op s_mpi_double_complex_sum;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Op_create(mpi_double_complex_sum, true, &s_mpi_double_complex_sum);
  }
  return s_mpi_double_complex_sum;
}

MPI_Op
op_identifier_compare_op()
{
  static MPI_Op s_op_identifier_op;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Op_create(mpi_op_identifier_compare, true, &s_op_identifier_op);
  }
  return s_op_identifier_op;
}

namespace {

const Parallel::ReduceSet *s_currentReduceSet = 0;

extern "C" {
  typedef void (*ParallelReduceOp)
  (void * inv, void * outv, int *, MPI_Datatype *);
}

void
all_reduce(
  MPI_Comm		arg_comm,
  ParallelReduceOp	arg_op,
  void *		arg_in,
  void *		arg_out,
  unsigned		arg_len)
{
  MPI_Op mpi_op = MPI_OP_NULL ;

  MPI_Op_create(arg_op, 0, & mpi_op);

  // The SUN was buggy when combining an
  // MPI_Allreduce with a user defined operator,
  // use reduce/broadcast instead.

  const int result = MPI_Allreduce(arg_in,arg_out,arg_len,MPI_BYTE,mpi_op,arg_comm);

  MPI_Op_free(& mpi_op);

  if (MPI_SUCCESS != result) {
    std::ostringstream msg ;
    msg << "Xyce::MPI::all_reduce FAILED: MPI_Allreduce = " << result;
    throw std::runtime_error(msg.str());
  }
}

struct ReduceCheck : public ReduceInterface
{
  ReduceCheck()
  {}

  void setSize(unsigned size) {
    m_size = size;
  }

  virtual void size(void *&inbuf) const {
    unsigned *t = align_cast<unsigned>(inbuf);
    t += sizeof(unsigned);
    inbuf = t;
  }

  virtual void copyin(void *&inbuf) const {
    unsigned *t = align_cast<unsigned>(inbuf);
    *t++ = m_size;
    inbuf = t;
  }

  virtual void copyout(void *&outbuf) const {
    unsigned *t = align_cast<unsigned>(outbuf);

    unsigned size = *t++;
    if (m_size != size)
      throw std::runtime_error("size mismatch");

    outbuf = t;
  }

  virtual void op(void *&inbuf, void *&outbuf) const {
    unsigned *tin = align_cast<unsigned>(inbuf);
    unsigned *tout = align_cast<unsigned>(outbuf);

    *tout = std::min(*tout, *tin);

    inbuf = ++tin;
    outbuf = ++tout;
  }

private:
  unsigned	m_size;
};

} // namespace <unnamed>


ReduceSet::ReduceSet()
{
  add(new ReduceCheck);
}


ReduceSet::~ReduceSet()
{
  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    delete (*it);
}


size_t
ReduceSet::size() const {
  void *buffer_end = 0;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->size(buffer_end);

  ReduceCheck *reduce_check = static_cast<ReduceCheck *>(m_reduceVector.front());
  reduce_check->setSize(reinterpret_cast<char *>(buffer_end) - (char *) 0);

  return reinterpret_cast<char *>(buffer_end) - (char *) 0;
}

void
ReduceSet::copyin(void * const buffer_in) const {
  void *inbuf = buffer_in;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->copyin(inbuf);
}

void
ReduceSet::copyout(void * const buffer_out) const {
  void *outbuf = buffer_out;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->copyout(outbuf);
}

void
ReduceSet::op(void * const buffer_in, void * const buffer_out) const {
  void *inbuf = buffer_in;
  void *outbuf = buffer_out;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->op(inbuf, outbuf);
}

void ReduceSet::void_op(void * inv, void * outv, int *, MPI_Datatype *) {
  s_currentReduceSet->op(inv, outv);
}


void
ReduceSet::add(
  ReduceInterface *	reduce_interface)
{
  m_reduceVector.push_back(reduce_interface);
}


void
AllReduce(
  MPI_Comm		comm,
  const ReduceSet &	reduce_set)
{
  size_t size = reduce_set.size();

  if (size) {
    char *input_buffer  = new char[size];
    char *output_buffer = new char[size];
    void *inbuf = (void *) input_buffer;
    void *outbuf = (void *) output_buffer;

    s_currentReduceSet = &reduce_set;

    ParallelReduceOp f = reinterpret_cast<ParallelReduceOp>(& ReduceSet::void_op);

    reduce_set.copyin(inbuf);
    all_reduce(comm, f, inbuf, outbuf, size);
    reduce_set.copyout(outbuf);
    delete [] output_buffer;
    delete [] input_buffer;
  }
}

void
AllWriteString(
  MPI_Comm              mpi_comm,
  std::ostream &        os,
  const std::string &   message)
{
  if (mpi_parallel_run(mpi_comm)) {
    const int root = 0;
    const unsigned size = Parallel::size(mpi_comm);
    const unsigned rank = Parallel::rank(mpi_comm);

    int result;

    // Gather the send counts on root processor
    int send_count = message.size();

    std::vector<int> recv_count(size, 0);
    int * const recv_count_ptr = &recv_count[0];

    result = MPI_Gather(& send_count, 1, MPI_INT,
                        recv_count_ptr, 1, MPI_INT,
                        root, mpi_comm);

    if (MPI_SUCCESS != result) {
      std::ostringstream message;
      message << "Parallel::AllWriteString FAILED: MPI_Gather = " << result;
      throw std::runtime_error(message.str());
    }

    // Receive counts are only non-zero on the root processor:
    std::vector<int> recv_displ(size + 1, 0);
    for (unsigned i = 0; i < size; ++i) {
      recv_displ[i + 1] = recv_displ[i] + recv_count[i];
    }

    const unsigned recv_size = (unsigned) recv_displ[ size ];
    std::vector<char> buffer(recv_size);

    {
      const char * const send_ptr = message.c_str();
      char * const recv_ptr = recv_size ? & buffer[0] : (char *) NULL;
      int * const recv_displ_ptr = & recv_displ[0];

      result = MPI_Gatherv((void*) send_ptr, send_count, MPI_CHAR,
                           recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                           root, mpi_comm);
    }

    if (MPI_SUCCESS != result) {
      std::ostringstream message;
      message << "Parallel::AllWriteString FAILED: MPI_Gatherv = " << result;
      throw std::runtime_error(message.str());
    }

    if (root == (int) rank) {
      for (unsigned i = 0; i < size; ++i) {
        if (recv_count[i]) {
          char * const ptr = & buffer[ recv_displ[i] ];
          os.write(ptr, recv_count[i]);
          os << std::endl;
        }
      }
      os.flush();
    }
  }
  else
    os << message;
}

void
GatherV(
  MPI_Comm                      mpi_comm,
  unsigned                      root,
  const std::string &           src,
  std::vector<std::string> &    dest)
{
  if (mpi_parallel_run(mpi_comm)) {
    unsigned size = Parallel::size(mpi_comm);
    unsigned rank = Parallel::rank(mpi_comm);

    int send_size = src.size();
    std::vector<int> receive_count;
    std::vector<int> receive_displacement;
    int *receive_count_ptr = 0;
    int *receive_displacement_ptr = 0;

    if (rank == root) {
      receive_count.resize(size);
      receive_displacement.resize(size);
      receive_count_ptr = &receive_count[0];
      receive_displacement_ptr = &receive_displacement[0];
    }

    MPI_Gather(&send_size, 1, MPI_INT, receive_count_ptr, 1, MPI_INT, root, mpi_comm);

    std::vector<char> receive;
    char *receive_ptr = 0;
    if (rank == root) {
      int displacement = 0;
      for (unsigned i = 0; i < size; ++i) {
        receive_displacement[i] = displacement;
        displacement += receive_count[i];
      }

      receive.resize(displacement);
      receive_ptr = &receive[0];
    }

    MPI_Gatherv(const_cast<char *>(src.data()), send_size, MPI_BYTE,
                receive_ptr, receive_count_ptr, receive_displacement_ptr, MPI_BYTE,
                root, mpi_comm);

    if (rank == root) {
      dest.resize(size);
      for (unsigned i = 0; i < size; ++i)
        dest[i] = std::string(&receive[receive_displacement[i]], &receive[receive_displacement[i] + receive_count[i]]);
    }
  }
  else {
    dest.resize(1);
    dest[0] = src;
  }
}

void
AllGatherV(
  MPI_Comm                      mpi_comm,
  const std::string &           src,
  std::vector<std::string> &    dest)
{
  if (mpi_parallel_run(mpi_comm)) {
    unsigned size = Parallel::size(mpi_comm);
    unsigned rank = Parallel::rank(mpi_comm);

    int send_size = src.size();
    std::vector<int> receive_count;
    std::vector<int> receive_displacement;
    int *receive_count_ptr = 0;
    int *receive_displacement_ptr = 0;

    receive_count.resize(size);
    receive_displacement.resize(size);
    receive_count_ptr = &receive_count[0];
    receive_displacement_ptr = &receive_displacement[0];

    MPI_Allgather(&send_size, 1, MPI_INT, receive_count_ptr, 1, MPI_INT, mpi_comm);

    std::vector<char> receive;
    char *receive_ptr = 0;
    int displacement = 0;
    for (unsigned i = 0; i < size; ++i) {
      receive_displacement[i] = displacement;
      displacement += receive_count[i];
    }

    receive.resize(displacement);
    receive_ptr = &receive[0];

    MPI_Allgatherv(const_cast<char *>(src.data()), send_size, MPI_BYTE,
                   receive_ptr, receive_count_ptr, receive_displacement_ptr, MPI_BYTE,
                   mpi_comm);

    dest.resize(size);
    for (unsigned i = 0; i < size; ++i)
      dest[i] = std::string(&receive[receive_displacement[i]], &receive[receive_displacement[i] + receive_count[i]]);
  }
  else {
    dest.resize(1);
    dest[0] = src;
  }
}

void
Broadcast(
  MPI_Comm      mpi_comm,
  std::string & s,
  int           root)
{
  if (mpi_parallel_run(mpi_comm)) {
    std::string::size_type size = s.length();
    if (MPI_Bcast(&size, 1, Datatype<std::string::size_type>::type(), root, mpi_comm) != MPI_SUCCESS)
      throw std::runtime_error("MPI_Bcast failed");
    if (Parallel::rank(mpi_comm) == root) {
      if (MPI_Bcast(const_cast<std::string::value_type *>(s.data()), size, Datatype<std::string::value_type>::type(), root, mpi_comm) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Bcast failed");
    }
    else {
      std::vector<std::string::value_type> t(size);

      if (MPI_Bcast(&t[0], size, Datatype<std::string::value_type>::type(), root, mpi_comm) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Bcast failed");
      s.assign(t.begin(), t.end());
    }
  }
}

void
Broadcast(
  MPI_Comm              mpi_comm,
  Util::Marshal &       m,
  int                   root)
{
  if (mpi_parallel_run(mpi_comm)) {
    std::string s;
    if (Parallel::rank(mpi_comm) == root)
      s = m.str();
    std::string::size_type size = s.length();
    if (MPI_Bcast(&size, 1, Datatype<std::string::size_type>::type(), root, mpi_comm) != MPI_SUCCESS)
      throw std::runtime_error("MPI_Bcast failed");
    if (Parallel::rank(mpi_comm) == root) {
      if (MPI_Bcast(const_cast<std::string::value_type *>(s.data()), size, Datatype<std::string::value_type>::type(), root, mpi_comm) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Bcast failed");
    }
    else {
      std::vector<std::string::value_type> t(size);

      if (MPI_Bcast(&t[0], size, Datatype<std::string::value_type>::type(), root, mpi_comm) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Bcast failed");
      s.assign(t.begin(), t.end());
      m.str(s);
    }
  }
}

} // namespace Parallel
} // namespace Xyced

#endif // Xyce_PARALLEL_MPI
