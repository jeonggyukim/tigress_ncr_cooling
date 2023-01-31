#ifndef IO_WRAPPER_HPP_
#define IO_WRAPPER_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file io_wrapper.hpp
//! \brief defines a set of small wrapper functions for MPI versus Serial Output.

// C headers

// C++ headers
#include <cstdio>

// Athena++ headers
// #include "../athena.hpp"
#include "defs.hpp"

using IOWrapperFile = FILE*;

using IOWrapperSizeT = std::uint64_t;

class IOWrapper {
 public:
  IOWrapper() {fh_=nullptr;}
  ~IOWrapper() {}
  // nested type definition of strongly typed/scoped enum in class definition
  enum class FileMode {read, write};

  // wrapper functions for basic I/O tasks
  int Open(const char* fname, FileMode rw);
  std::size_t Read(void *buf, IOWrapperSizeT size, IOWrapperSizeT count);
  std::size_t Read_all(void *buf, IOWrapperSizeT size, IOWrapperSizeT count);
  std::size_t Read_at_all(void *buf, IOWrapperSizeT size,
                          IOWrapperSizeT count, IOWrapperSizeT offset);
  std::size_t Write(const void *buf, IOWrapperSizeT size, IOWrapperSizeT count);
  std::size_t Write_at_all(const void *buf, IOWrapperSizeT size,
                           IOWrapperSizeT cnt, IOWrapperSizeT offset);
  int Close();
  int Seek(IOWrapperSizeT offset);
  IOWrapperSizeT GetPosition();

 private:
  IOWrapperFile fh_;

};
#endif // IO_WRAPPER_HPP_
