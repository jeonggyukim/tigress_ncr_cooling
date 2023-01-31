//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file parameter_input.cpp
//! \brief implementation of functions in class ParameterInput
//!
//! PURPOSE
//! =======
//! Member functions of this class are used to read and parse the input file.
//! Functionality is loosely modeled after FORTRAN namelist.
//!
//! EXAMPLE
//! =======
//! Input file in 'Athena++' format:
//!
//!     <blockname1>      # block name; must be on a line by itself
//!                       # everything after a hash symbol is a comment and is ignored
//!     name1=value       # each parameter name must be on a line by itself
//!     name2 = value1    # whitespace around the = is optional
//!                       # blank lines are OK
//!     # my comment here   comment lines are OK
//!     # name3 = value3    values (and blocks) that are commented out are ignored
//!
//!     <blockname2>      # start new block
//!     name1 = value1    # note that same parameter names can appear in different blocks
//!     name2 = value2    # empty lines (like following) are OK
//!
//!     <blockname1>      # same blockname can re-appear, although NOT recommended
//!     name3 = value3    # this would be the 3rd parameter name in blockname1
//!     name1 = value4    # if parameter name is repeated, previous value is overwritten!
//!
//! LIMITATIONS
//! ===========
//!  - parameter specification (name=val # comment) must all be on a single line
//!
//! HISTORY
//! =======
//!  - Nov 2002:  Created for Athena1.0/Cambridge release by Peter Teuben
//!  - 2003-2008: Many improvements and extensions by T. Gardiner and J.M. Stone
//!  - Jan 2014:  Rewritten in C++ for the Athena++ code by J.M. Stone
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // transform
#include <cmath>      // std::fmod()
#include <cstdlib>    // atoi(), atof(), nullptr, std::size_t
#include <fstream>    // ifstream
#include <iostream>   // endl, ostream
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
// ##include "athena.hpp"
#include "defs.hpp"
//#include "globals.hpp"
#include "parameter_input.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn ParameterInput::ParameterInput()
//! \brief ParameterInput constructor

ParameterInput::ParameterInput() :pfirst_block{}, last_filename_{} {
#ifdef OPENMP_PARALLEL
  omp_init_lock(&lock_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn ParameterInput::~ParameterInput()
//! \brief ParameterInput destructor
//!
//! iterates through nested singly linked lists of blocks/lines and deletes each
//! InputBlock node (whose destructor below deletes linked list "line" nodes)

ParameterInput::~ParameterInput() {
  InputBlock *pib = pfirst_block;
  while (pib != nullptr) {
    InputBlock *pold_block = pib;
    pib = pib->pnext;
    delete pold_block;
  }
#ifdef OPENMP_PARALLEL
  omp_destroy_lock(&lock_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn InputBlock::~InputBlock()
//! \brief InputBlock destructor
//!
//! iterates through singly linked list of "line" nodes and deletes them

InputBlock::~InputBlock() {
  InputLine *pil = pline;
  while (pil != nullptr) {
    InputLine *pold_line = pil;
    pil = pil->pnext;
    delete pold_line;
  }
}

//----------------------------------------------------------------------------------------
//! \fn  void ParameterInput::LoadFromStream(std::istream &is)
//! \brief Load input parameters from a stream
//!
//! Input block names are allocated and stored in a singly linked list of InputBlocks.
//! Within each InputBlock the names, values, and comments of each parameter are allocated
//! and stored in a singly linked list of InputLines.

void ParameterInput::LoadFromStream(std::istream &is) {
  std::string line, block_name, param_name, param_value, param_comment;
  std::size_t first_char, last_char;
  std::stringstream msg;
  InputBlock *pib{};
  int line_num{-1}, blocks_found{0};

  while (is.good()) {
    std::getline(is, line);
    line_num++;
    if (line.find('\t') != std::string::npos) {
      line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
      // msg << "### FATAL ERROR in function [ParameterInput::LoadFromStream]"
      //     << std::endl << "Tab characters are forbidden in input files";
      // ATHENA_ERROR(msg);
    }
    if (line.empty()) continue;                             // skip blank line
    first_char = line.find_first_not_of(" ");               // skip white space
    if (first_char == std::string::npos) continue;          // line is all white space
    if (line.compare(first_char, 1, "#") == 0) continue;      // skip comments
    if (line.compare(first_char, 9, "<par_end>") == 0) break; // stop on <par_end>

    if (line.compare(first_char, 1, "<") == 0) {              // a new block
      first_char++;
      last_char = (line.find_first_of(">", first_char));
      block_name.assign(line, first_char, last_char-1);       // extract block name

      if (last_char == std::string::npos) {
        msg << "### FATAL ERROR in function [ParameterInput::LoadFromStream]"
            << std::endl << "Block name '" << block_name
            << "' in the input stream'" << "' not properly ended";
        ATHENA_ERROR(msg);
      }

      pib = FindOrAddBlock(block_name);  // find or add block to singly linked list

      if (pib == nullptr) {
        msg << "### FATAL ERROR in function [ParameterInput::LoadFromStream]"
            << std::endl << "Block name '" << block_name
            << "' could not be found/added";
        ATHENA_ERROR(msg);
      }
      blocks_found++;
      continue;  // skip to next line if block name was found
    } // end "a new block was found"

    // if line does not contain a block name or skippable information (comments,
    // whitespace), it must contain a parameter value
    if (blocks_found == 0) {
        msg << "### FATAL ERROR in function [ParameterInput::LoadFromStream]"
            << std::endl << "Input file must specify a block name before the first"
            << " parameter = value line";
        ATHENA_ERROR(msg);
    }
    // parse line and add name/value/comment strings (if found) to current block name
    ParseLine(pib, line, param_name, param_value, param_comment);
    AddParameter(pib, param_name, param_value, param_comment);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void ParameterInput::LoadFromFile(IOWrapper &input)
//! \brief Read the parameters from an input file or restarting file.
//!
//! Return the position at the end of the header, which is used in restarting

void ParameterInput::LoadFromFile(IOWrapper &input) {
  std::stringstream par, msg;
  constexpr int kBufSize = 4096;
  char buf[kBufSize];
  IOWrapperSizeT header = 0, ret, loc;

  // search <par_end> or EOF.
  do {
    // if (Globals::my_rank == 0) // only the master process reads the header from the file
      ret = input.Read(buf, sizeof(char), kBufSize);
// #ifdef MPI_PARALLEL
//     // then broadcasts it
//     MPI_Bcast(&ret, sizeof(IOWrapperSizeT), MPI_BYTE, 0, MPI_COMM_WORLD);
//     MPI_Bcast(buf, ret, MPI_BYTE, 0, MPI_COMM_WORLD);
// #endif
    par.write(buf, ret); // add the buffer into the stream
    header += ret;
    std::string sbuf = par.str(); // create string for search
    loc = sbuf.find("<par_end>", 0); // search from the top of the stream
    if (loc != std::string::npos) { // found <par_end>
      header = loc + 10; // store the header length
      break;
    }
    if (header > kBufSize*10) {
      msg << "### FATAL ERROR in function [ParameterInput::LoadFromFile]"
          << "<par_end> is not found in the first 40KBytes." << std::endl
          << "Probably the file is broken or a wrong file is specified" << std::endl;
      ATHENA_ERROR(msg);
    }
  } while (ret == kBufSize); // till EOF (or par_end is found)

  // Now par contains the parameter inputs + some additional including <par_end>
  // Read the stream and load the parameters
  LoadFromStream(par);
  // Seek the file to the end of the header
  input.Seek(header);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn InputBlock* ParameterInput::FindOrAddBlock(std::string name)
//! \brief find or add specified InputBlock.  Returns pointer to block.

InputBlock* ParameterInput::FindOrAddBlock(std::string name) {
  InputBlock *pib, *plast;
  plast = pfirst_block;
  pib = pfirst_block;

  // Search singly linked list of InputBlocks to see if name exists, return if found.
  while (pib != nullptr) {
    if (name.compare(pib->block_name) == 0) return pib;
    plast = pib;
    pib = pib->pnext;
  }

  // Create new block in list if not found above
  pib = new InputBlock;
  pib->block_name.assign(name);  // store the new block name
  pib->pline = nullptr;             // Terminate the InputLine list
  pib->pnext = nullptr;             // Terminate the InputBlock list

  // if this is the first block in list, save pointer to it in class
  if (pfirst_block == nullptr) {
    pfirst_block = pib;
  } else {
    plast->pnext = pib;      // link new node into list
  }

  return pib;
}

//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::ParseLine(InputBlock *pib, std::string line,
//!          std::string& name, std::string& value, std::string& comment)
//! \brief parse "name = value # comment" format, return name/value/comment strings.

void ParameterInput::ParseLine(InputBlock *pib, std::string line,
                               std::string& name, std::string& value,
                               std::string& comment) {
  std::size_t first_char, last_char, equal_char, hash_char, len;

  first_char = line.find_first_not_of(" ");   // find first non-white space
  equal_char = line.find_first_of("=");       // find "=" char
  hash_char  = line.find_first_of("#");       // find "#" (optional)

  // copy substring into name, remove white space at end of name
  len = equal_char - first_char;
  name.assign(line, first_char, len);

  last_char = name.find_last_not_of(" ");
  name.erase(last_char+1, std::string::npos);

  // copy substring into value, remove white space at start and end
  len = hash_char - equal_char - 1;
  value.assign(line, equal_char+1, len);

  first_char = value.find_first_not_of(" ");
  value.erase(0, first_char);

  last_char = value.find_last_not_of(" ");
  value.erase(last_char+1, std::string::npos);

  // copy substring into comment, if present
  if (hash_char != std::string::npos) {
    comment = line.substr(hash_char);
  } else {
    comment = "";
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::AddParameter(InputBlock *pb, std::string name,
//!  std::string value, std::string comment)
//! \brief add name/value/comment tuple to the InputLine singly linked list in block *pb.
//!
//! If a parameter with the same name already exists, the value and comment strings
//! are replaced (overwritten).

void ParameterInput::AddParameter(InputBlock *pb, std::string name,
                                  std::string value, std::string comment) {
  InputLine *pl, *plast;
  // Search singly linked list of InputLines to see if name exists.  This also sets *plast
  // to point to the tail node (but not storing a pointer to the tail node in InputBlock)
  pl = pb->pline;
  plast = pb->pline;
  while (pl != nullptr) {
    if (name.compare(pl->param_name) == 0) {   // param name already exists
      pl->param_value.assign(value);           // replace existing param value
      pl->param_comment.assign(comment);       // replace exisiting param comment
      if (value.length() > pb->max_len_parvalue) pb->max_len_parvalue = value.length();
      return;
    }
    plast = pl;
    pl = pl->pnext;
  }

  // Create new node in singly linked list if name does not already exist
  pl = new InputLine;
  pl->param_name.assign(name);
  pl->param_value.assign(value);
  pl->param_comment.assign(comment);
  pl->pnext = nullptr;

  // if this is the first parameter in list, save pointer to it in block.
  if (pb->pline == nullptr) {
    pb->pline = pl;
    pb->max_len_parname = name.length();
    pb->max_len_parvalue = value.length();
  } else {
    plast->pnext = pl;  // link new node into list
    if (name.length() > pb->max_len_parname) pb->max_len_parname = name.length();
    if (value.length() > pb->max_len_parvalue) pb->max_len_parvalue = value.length();
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::ModifyFromCmdline(int argc, char *argv[])
//! \brief parse commandline for changes to input parameters
//!
//! Note this function is very forgiving (no warnings!) if there is an error in format

void ParameterInput::ModifyFromCmdline(int argc, char *argv[]) {
  std::string input_text, block,name, value;
  std::stringstream msg;
  InputBlock *pb;
  InputLine *pl;

  for (int i=1; i<argc; i++) {
    input_text = argv[i];
    std::size_t slash_posn = input_text.find_first_of("/");   // find "/" character
    std::size_t equal_posn = input_text.find_first_of("=");   // find "=" character

    // skip if either "/" or "=" do not exist in input
    if ((slash_posn == std::string::npos) || (equal_posn == std::string::npos)) continue;

    // extract block/name/value strings
    block = input_text.substr(0, slash_posn);
    name  = input_text.substr(slash_posn+1, (equal_posn - slash_posn - 1));
    value = input_text.substr(equal_posn+1, std::string::npos);

    // get pointer to node with same block name in singly linked list of InputBlocks
    pb = GetPtrToBlock(block);
    if (pb == nullptr) {
      msg << "### FATAL ERROR in function [ParameterInput::ModifyFromCmdline]"
          << std::endl << "Block name '" << block << "' on command line not found";
      ATHENA_ERROR(msg);
    }

    // get pointer to node with same parameter name in singly linked list of InputLines
    pl = pb->GetPtrToLine(name);
    if (pl == nullptr) {
      msg << "### FATAL ERROR in function [ParameterInput::ModifyFromCmdline]"
          << std::endl << "Parameter '" << name << "' in block '" << block
          << "' on command line not found";
      ATHENA_ERROR(msg);
    }
    pl->param_value.assign(value);   // replace existing value

    if (value.length() > pb->max_len_parvalue) pb->max_len_parvalue = value.length();
  }
}

//----------------------------------------------------------------------------------------
//! \fn InputBlock* ParameterInput::GetPtrToBlock(std::string name)
//! \brief return pointer to specified InputBlock if it exists

InputBlock* ParameterInput::GetPtrToBlock(std::string name) {
  InputBlock *pb;
  for (pb = pfirst_block; pb != nullptr; pb = pb->pnext) {
    if (name.compare(pb->block_name) == 0) return pb;
  }
  return nullptr;
}

//----------------------------------------------------------------------------------------
//! \fn int ParameterInput::DoesParameterExist(std::string block, std::string name)
//! \brief check whether parameter of given name in given block exists

int ParameterInput::DoesParameterExist(std::string block, std::string name) {
  InputLine *pl;
  InputBlock *pb;
  pb = GetPtrToBlock(block);
  if (pb == nullptr) return 0;
  pl = pb->GetPtrToLine(name);
  return (pl == nullptr ? 0 : 1);
}

//----------------------------------------------------------------------------------------
//! \fn int ParameterInput::GetInteger(std::string block, std::string name)
//! \brief returns integer value of string stored in block/name

int ParameterInput::GetInteger(std::string block, std::string name) {
  InputBlock* pb;
  InputLine* pl;
  std::stringstream msg;

  Lock();

  // get pointer to node with same block name in singly linked list of InputBlocks
  pb = GetPtrToBlock(block);
  if (pb == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetInteger]" << std::endl
        << "Block name '" << block << "' not found when trying to set value "
        << "for parameter '" << name << "'";
    ATHENA_ERROR(msg);
  }

  // get pointer to node with same parameter name in singly linked list of InputLines
  pl = pb->GetPtrToLine(name);
  if (pl == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetInteger]" << std::endl
        << "Parameter name '" << name << "' not found in block '" << block << "'";
    ATHENA_ERROR(msg);
  }

  std::string val=pl->param_value;
  Unlock();

  // Convert string to integer and return value
  return atoi(val.c_str());
}

//----------------------------------------------------------------------------------------
//! \fn Real ParameterInput::GetReal(std::string block, std::string name)
//! \brief returns real value of string stored in block/name

Real ParameterInput::GetReal(std::string block, std::string name) {
  InputBlock* pb;
  InputLine* pl;
  std::stringstream msg;

  Lock();

  // get pointer to node with same block name in singly linked list of InputBlocks
  pb = GetPtrToBlock(block);
  if (pb == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Block name '" << block << "' not found when trying to set value "
        << "for parameter '" << name << "'";
    ATHENA_ERROR(msg);
  }

  // get pointer to node with same parameter name in singly linked list of InputLines
  pl = pb->GetPtrToLine(name);
  if (pl == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Parameter name '" << name << "' not found in block '" << block << "'";
    ATHENA_ERROR(msg);
  }

  std::string val=pl->param_value;
  Unlock();

  // Convert string to real and return value
  return static_cast<Real>(atof(val.c_str()));
}

//----------------------------------------------------------------------------------------
//! \fn bool ParameterInput::GetBoolean(std::string block, std::string name)
//! \brief returns boolean value of string stored in block/name

bool ParameterInput::GetBoolean(std::string block, std::string name) {
  InputBlock* pb;
  InputLine* pl;
  std::stringstream msg;

  Lock();

  // get pointer to node with same block name in singly linked list of InputBlocks
  pb = GetPtrToBlock(block);
  if (pb == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Block name '" << block << "' not found when trying to set value "
        << "for parameter '" << name << "'";
    ATHENA_ERROR(msg);
  }

  // get pointer to node with same parameter name in singly linked list of InputLines
  pl = pb->GetPtrToLine(name);
  if (pl == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Parameter name '" << name << "' not found in block '" << block << "'";
    ATHENA_ERROR(msg);
  }

  std::string val=pl->param_value;
  Unlock();

  // check is string contains integers 0 or 1 (instead of true or false) and return
  if (val.compare(0, 1, "0")==0 || val.compare(0, 1, "1")==0) {
    return static_cast<bool>(atoi(val.c_str()));
  }

  // convert string to all lower case
  std::transform(val.begin(), val.end(), val.begin(), ::tolower);
  // Convert string to bool and return value
  bool b;
  std::istringstream is(val);
  is >> std::boolalpha >> b;

  return (b);
}

//----------------------------------------------------------------------------------------
//! \fn std::string ParameterInput::GetString(std::string block, std::string name)
//! \brief returns string stored in block/name

std::string ParameterInput::GetString(std::string block, std::string name) {
  InputBlock* pb;
  InputLine* pl;
  std::stringstream msg;

  Lock();

  // get pointer to node with same block name in singly linked list of InputBlocks
  pb = GetPtrToBlock(block);
  if (pb == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Block name '" << block << "' not found when trying to set value "
        << "for parameter '" << name << "'";
    ATHENA_ERROR(msg);
  }

  // get pointer to node with same parameter name in singly linked list of InputLines
  pl = pb->GetPtrToLine(name);
  if (pl == nullptr) {
    msg << "### FATAL ERROR in function [ParameterInput::GetReal]" << std::endl
        << "Parameter name '" << name << "' not found in block '" << block << "'";
    ATHENA_ERROR(msg);
  }

  std::string val=pl->param_value;
  Unlock();

  // return value
  return val;
}

//----------------------------------------------------------------------------------------
//! \fn int ParameterInput::GetOrAddInteger(std::string block, std::string name,
//!   int default_value)
//! \brief returns integer value stored in block/name if it exists, or creates and sets
//!   value to def_value if it does not exist

int ParameterInput::GetOrAddInteger(std::string block, std::string name, int def_value) {
  InputBlock* pb;
  InputLine *pl;
  std::stringstream ss_value;
  int ret;

  Lock();
  if (DoesParameterExist(block, name)) {
    pb = GetPtrToBlock(block);
    pl = pb->GetPtrToLine(name);
    std::string val = pl->param_value;
    ret = atoi(val.c_str());
  } else {
    pb = FindOrAddBlock(block);
    ss_value << def_value;
    AddParameter(pb, name, ss_value.str(), "# Default value added at run time");
    ret = def_value;
  }
  Unlock();
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn Real ParameterInput::GetOrAddReal(std::string block, std::string name,
//!   Real def_value)
//! \brief returns real value stored in block/name if it exists, or creates and sets
//!   value to def_value if it does not exist

Real ParameterInput::GetOrAddReal(std::string block, std::string name, Real def_value) {
  InputBlock* pb;
  InputLine *pl;
  std::stringstream ss_value;
  Real ret;

  Lock();
  if (DoesParameterExist(block, name)) {
    pb = GetPtrToBlock(block);
    pl = pb->GetPtrToLine(name);
    std::string val = pl->param_value;
    ret = static_cast<Real>(atof(val.c_str()));
  } else {
    pb = FindOrAddBlock(block);
    ss_value << def_value;
    AddParameter(pb, name, ss_value.str(), "# Default value added at run time");
    ret = def_value;
  }
  Unlock();
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn bool ParameterInput::GetOrAddBoolean(std::string block, std::string name,
//!   bool def_value)
//! \brief returns boolean value stored in block/name if it exists, or creates and sets
//!   value to def_value if it does not exist

bool ParameterInput::GetOrAddBoolean(std::string block,std::string name, bool def_value) {
  InputBlock* pb;
  InputLine *pl;
  std::stringstream ss_value;
  bool ret;

  Lock();
  if (DoesParameterExist(block, name)) {
    pb = GetPtrToBlock(block);
    pl = pb->GetPtrToLine(name);
    std::string val = pl->param_value;
    if (val.compare(0, 1, "0")==0 || val.compare(0, 1, "1")==0) {
      ret = static_cast<bool>(atoi(val.c_str()));
    } else {
      std::transform(val.begin(), val.end(), val.begin(), ::tolower);
      std::istringstream is(val);
      is >> std::boolalpha >> ret;
    }
  } else {
    pb = FindOrAddBlock(block);
    ss_value << def_value;
    AddParameter(pb, name, ss_value.str(), "# Default value added at run time");
    ret = def_value;
  }
  Unlock();
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn std::string ParameterInput::GetOrAddString(std::string block, std::string name,
//!                                                std::string def_value)
//! \brief returns string value stored in block/name if it exists, or creates and sets
//!   value to def_value if it does not exist

std::string ParameterInput::GetOrAddString(std::string block, std::string name,
                                           std::string def_value) {
  InputBlock* pb;
  InputLine *pl;
  std::stringstream ss_value;
  std::string ret;

  Lock();
  if (DoesParameterExist(block, name)) {
    pb = GetPtrToBlock(block);
    pl = pb->GetPtrToLine(name);
    ret = pl->param_value;
  } else {
    pb = FindOrAddBlock(block);
    AddParameter(pb, name, def_value, "# Default value added at run time");
    ret = def_value;
  }
  Unlock();
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn int ParameterInput::SetInteger(std::string block, std::string name, int value)
//! \brief updates an integer parameter; creates it if it does not exist

int ParameterInput::SetInteger(std::string block, std::string name, int value) {
  InputBlock* pb;
  std::stringstream ss_value;

  Lock();
  pb = FindOrAddBlock(block);
  ss_value << value;
  AddParameter(pb, name, ss_value.str(), "# Updated during run time");
  Unlock();
  return value;
}

//----------------------------------------------------------------------------------------
//! \fn Real ParameterInput::SetReal(std::string block, std::string name, Real value)
//! \brief updates a real parameter; creates it if it does not exist

Real ParameterInput::SetReal(std::string block, std::string name, Real value) {
  InputBlock* pb;
  std::stringstream ss_value;

  Lock();
  pb = FindOrAddBlock(block);
  ss_value << value;
  AddParameter(pb, name, ss_value.str(), "# Updated during run time");
  Unlock();
  return value;
}

//----------------------------------------------------------------------------------------
//! \fn bool ParameterInput::SetBoolean(std::string block, std::string name, bool value)
//! \brief updates a boolean parameter; creates it if it does not exist

bool ParameterInput::SetBoolean(std::string block, std::string name, bool value) {
  InputBlock* pb;
  std::stringstream ss_value;

  Lock();
  pb = FindOrAddBlock(block);
  ss_value << value;
  AddParameter(pb, name, ss_value.str(), "# Updated during run time");
  Unlock();
  return value;
}

//----------------------------------------------------------------------------------------
//! \fn std::string ParameterInput::SetString(std::string block, std::string name,
//!                                           std::string  value)
//! \brief updates a string parameter; creates it if it does not exist

std::string ParameterInput::SetString(std::string block, std::string name,
                                      std::string value) {
  InputBlock* pb;

  Lock();
  pb = FindOrAddBlock(block);
  AddParameter(pb, name, value, "# Updated during run time");
  Unlock();
  return value;
}

//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::RollbackNextTime()
//! \brief rollback next_time by dt for each output block

void ParameterInput::RollbackNextTime() {
  InputBlock *pb = pfirst_block;
  InputLine* pl;
  std::stringstream msg;
  Real next_time;

  while (pb != nullptr) {
    if (pb->block_name.compare(0, 6, "output") == 0) {
      pl = pb->GetPtrToLine("next_time");
      if (pl == nullptr) {
        msg << "### FATAL ERROR in function [ParameterInput::RollbackNextTime]"
            << std::endl << "Parameter name 'next_time' not found in block '"
            << pb->block_name << "'";
        ATHENA_ERROR(msg);
      }
      next_time = static_cast<Real>(atof(pl->param_value.c_str()));
      pl = pb->GetPtrToLine("dt");
      if (pl == nullptr) {
        msg << "### FATAL ERROR in function [ParameterInput::RollbackNextTime]"
            << std::endl << "Parameter name 'dt' not found in block '"
            << pb->block_name << "'";
        ATHENA_ERROR(msg);
      }
      next_time -= static_cast<Real>(atof(pl->param_value.c_str()));
      msg << next_time;
      //AddParameter(pb, "next_time", msg.str().c_str(), "# Updated during run time");
      SetReal(pb->block_name, "next_time", next_time);
    }
    pb = pb->pnext;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::ForwardNextTime()
//! \brief add dt to next_time until next_time >  mesh_time - dt for each output block

void ParameterInput::ForwardNextTime(Real mesh_time) {
  InputBlock *pb = pfirst_block;
  InputLine* pl;
  Real next_time;
  Real dt0, dt;
  bool fresh = false;

  while (pb != nullptr) {
    if (pb->block_name.compare(0, 6, "output") == 0) {
      std::stringstream msg;
      pl = pb->GetPtrToLine("next_time");
      if (pl == nullptr) {
        next_time = mesh_time;
        // This is a freshly added output
        fresh = true;
      } else {
        next_time = static_cast<Real>(atof(pl->param_value.c_str()));
      }
      pl = pb->GetPtrToLine("dt");
      if (pl == nullptr) {
        msg << "### FATAL ERROR in function [ParameterInput::ForwardNextTime]"
            << std::endl << "Parameter name 'dt' not found in block '"
            << pb->block_name << "'";
        ATHENA_ERROR(msg);
      }
      dt0 = static_cast<Real>(atof(pl->param_value.c_str()));
      dt = dt0 * static_cast<int>((mesh_time - next_time) / dt0) + dt0;
      if (dt > 0) {
        next_time += dt;
        // If the user has added a new/fresh output round to multiple of dt0,
        // and make sure that mesh_time - dt0 < next_time < mesh_time,
        // to ensure immediate writing
        if (fresh) next_time -= std::fmod(next_time, dt0) + dt0;
      }
      msg << next_time;
      AddParameter(pb, "next_time", msg.str().c_str(), "# Updated during run time");
    }
    pb = pb->pnext;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::ParameterDump(std::ostream& os)
//! \brief output entire InputBlock/InputLine hierarchy to specified stream

void ParameterInput::ParameterDump(std::ostream& os) {
  InputBlock *pb;
  InputLine *pl;
  std::string param_name,param_value;
  std::size_t len;

  os<< "#------------------------- PAR_DUMP -------------------------" << std::endl;

  for (pb = pfirst_block; pb != nullptr; pb = pb->pnext) { // loop over InputBlocks
    os<< "<" << pb->block_name << ">" << std::endl;     // write block name
    for (pl = pb->pline; pl != nullptr; pl = pl->pnext) {   // loop over InputLines
      param_name.assign(pl->param_name);
      param_value.assign(pl->param_value);

      len = pb->max_len_parname - param_name.length() + 1;
      param_name.append(len,' ');                      // pad name to align vertically
      len = pb->max_len_parvalue - param_value.length() + 1;
      param_value.append(len,' ');                     // pad value to align vertically

      os<< param_name << "= " << param_value << pl->param_comment <<  std::endl;
    }
  }

  os<< "#------------------------- PAR_DUMP -------------------------" << std::endl;
  os<< "<par_end>" << std::endl;    // finish with par-end (useful in restart files)
}

//----------------------------------------------------------------------------------------
//! \fn InputLine* InputBlock::GetPtrToLine(std::string name)
//! \brief return pointer to InputLine containing specified parameter if it exists

InputLine* InputBlock::GetPtrToLine(std::string name) {
  for (InputLine* pl = pline; pl != nullptr; pl = pl->pnext) {
    if (name.compare(pl->param_name) == 0) return pl;
  }
  return nullptr;
}


//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::Lock()
//! \brief Lock ParameterInput for reading and writing
void ParameterInput::Lock() {
#ifdef OPENMP_PARALLEL
  omp_set_lock(&lock_);
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ParameterInput::Unlock()
//! \brief Unlock ParameterInput for reading and writing
void ParameterInput::Unlock() {
#ifdef OPENMP_PARALLEL
  omp_unset_lock(&lock_);
#endif
  return;
}
