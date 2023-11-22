//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file io_wrapper.cpp

// C headers

// C++ headers
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
// #include "../athena.hpp"
#include "io_wrapper.hpp"

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname, FileMode rw)

int IOWrapper::Open(const char* fname, FileMode rw) {
  std::stringstream msg;

  if (rw == FileMode::read) {
      if ((fh_ = std::fopen(fname,"rb")) == nullptr) // NOLINT
      {
        msg << "### FATAL ERROR in function [IOWrapper:Open]"
            <<std::endl<< "Input file '" << fname << "' could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
        return false;
      }

  } else if (rw == FileMode::write) {
      if ((fh_ = std::fopen(fname,"wb")) == nullptr) // NOLINT
      {
        msg << "### FATAL ERROR in function [IOWrapper:Open]"
            << std::endl<< "Output file '" << fname
            << "' could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
        return false;
      }
  } else {
    return false;
  }

  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read(void *buf, IOWrapperSizeT size, IOWrapperSizeT count)

std::size_t IOWrapper::Read(void *buf, IOWrapperSizeT size, IOWrapperSizeT count) {
  return std::fread(buf,size,count,fh_);
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_all(void *buf, IOWrapperSizeT size, IOWrapperSizeT count)

std::size_t IOWrapper::Read_all(void *buf, IOWrapperSizeT size, IOWrapperSizeT count) {
  return std::fread(buf,size,count,fh_);
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_at_all(void *buf, IOWrapperSizeT size,
//!                            IOWrapperSizeT count, IOWrapperSizeT offset)

std::size_t IOWrapper::Read_at_all(void *buf, IOWrapperSizeT size,
                                   IOWrapperSizeT count, IOWrapperSizeT offset) {
  std::fseek(fh_, offset, SEEK_SET);
  return std::fread(buf,size,count,fh_);
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write(const void *buf, IOWrapperSizeT size, IOWrapperSizeT cnt)

std::size_t IOWrapper::Write(const void *buf, IOWrapperSizeT size, IOWrapperSizeT cnt) {
  return std::fwrite(buf,size,cnt,fh_);
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write_at_all(const void *buf, IOWrapperSizeT size,
//!                                 IOWrapperSizeT cnt, IOWrapperSizeT offset)

std::size_t IOWrapper::Write_at_all(const void *buf, IOWrapperSizeT size,
                                    IOWrapperSizeT cnt, IOWrapperSizeT offset) {
  std::fseek(fh_, offset, SEEK_SET);
  return std::fwrite(buf,size,cnt,fh_);
}


//----------------------------------------------------------------------------------------
//! \fn void IOWrapper::Close()

int IOWrapper::Close() {
  return std::fclose(fh_);
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Seek(IOWrapperSizeT offset)

int IOWrapper::Seek(IOWrapperSizeT offset) {
  return std::fseek(fh_, offset, SEEK_SET);
}

//----------------------------------------------------------------------------------------
//! \fn IOWrapperSizeT IOWrapper::GetPosition()

IOWrapperSizeT IOWrapper::GetPosition() {
  return ftell(fh_);
}
